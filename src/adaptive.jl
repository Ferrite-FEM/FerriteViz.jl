# View-adaptive tessellation: the glue between FEData and the isubd core
# (src/isubd.jl), and the compute-graph wiring for `solutionplot(...;
# adaptive=true)` (#161).
#
# Graph topology (all nodes prefixed subd_ to avoid clashes):
#
#   camera (projectionview, eyeposition, resolution)  ─┐
#   px_target, max_depth                               ├─> subd_keys ──> subd_positions, subd_ξ, subd_faces
#   ds_u (into subd_keys only when the dataset warps) ─┘        │              │
#   ds_u ──────────────────────────────────────────────> subd_color <── subd_ξ, subd_keys
#
# The persistent LEB key buffer lives in plot-local state captured by the
# subd_keys computation; the node returns `nothing` when a camera move does
# not actually change the key set, so downstream decode and transfer are
# skipped. With a static camera the inputs never dirty and nothing at all
# recomputes.

# One dof field, evaluated at arbitrary (cell, ξ) — the continuous-geometry
# and per-vertex-color workhorse. PointValues are reinitialized per query, the
# same cost profile as transfer_solution (#153 notes the potential speedup).
struct FieldEvaluator{PV}
    pvs::Vector{PV}                   # one PointValues per subdofhandler with the field
    sdh_of_cell::Vector{Int}          # cell -> index into pvs, 0 when the field is absent
    celldofs_field::Vector{Vector{Int}}  # cell -> global dofs of the field, empty when absent
    ncomps::Int
end

function FieldEvaluator(dh::Ferrite.DofHandler, field::Symbol)
    sdhs = getsubdofhandlers(dh, field)
    isempty(sdhs) && error("field :$field not found in the DofHandler")
    ip_field = Ferrite.getfieldinterpolation(first(sdhs), field)
    ξ0 = Ferrite.Vec(ntuple(d -> 0.0, Ferrite.getrefdim(ip_field)))
    ncomps = length(Ferrite.reference_shape_value(ip_field, ξ0, 1))
    pvs = [Ferrite.PointValues(Ferrite.getfieldinterpolation(sdh, field),
                               Ferrite.geometric_interpolation(Ferrite.getcelltype(sdh));
                               update_gradients=false) for sdh in sdhs]
    ncells = Ferrite.getncells(Ferrite.get_grid(dh))
    sdh_of_cell = zeros(Int, ncells)
    celldofs_field = [Int[] for _ in 1:ncells]
    for (si, sdh) in enumerate(sdhs)
        rng = Ferrite.dof_range(sdh, field)
        for cell_idx in sdh.cellset
            sdh_of_cell[cell_idx] = si
            celldofs_field[cell_idx] = Ferrite.celldofs(dh, cell_idx)[rng]
        end
    end
    return FieldEvaluator(pvs, sdh_of_cell, celldofs_field, ncomps)
end

# Evaluate at one reference point of one cell; `nothing` outside the field's
# subdomain. `cellcoords` are the cell's node coordinates (caller-cached).
function evaluate_at(ev::FieldEvaluator, cell_idx::Int, cellcoords, ξ, u::AbstractVector)
    si = ev.sdh_of_cell[cell_idx]
    si == 0 && return nothing
    pv = ev.pvs[si]
    Ferrite.reinit!(pv, cellcoords, ξ)
    return Ferrite.function_value(pv, 1, @views(u[ev.celldofs_field[cell_idx]]))
end

# The isubd base domain of a dataset: the *unsubdivided* reference tessellation
# of every visible cell, one LEB-ordered corner triple per base triangle, and
# the triangle -> cell map. Static per plot, like the visibility mask.
function _isubd_base_triangles(ds::FEData)
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    refdim = Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(first(cells))))
    corners = NTuple{3,Ferrite.Vec{refdim,Float64}}[]
    cellmap = Int[]
    for (cell_id, cell) in enumerate(cells)
        ds.visible[cell_id] || continue
        Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(cell))) == refdim ||
            error("adaptive tessellation requires a single reference dimension across the grid")
        tess = reference_tessellation(getrefshape(cell))
        for tri in tess.triangles
            push!(corners, leb_order((tess.coords[tri[1]], tess.coords[tri[2]], tess.coords[tri[3]])))
            push!(cellmap, cell_id)
        end
    end
    isempty(corners) && error("adaptive tessellation found no visible cells to tessellate")
    return corners, cellmap
end

# Continuous geometry x(ξ) [+ Σ scaleᵢ · fieldᵢ(ξ) for upstream warps] of one
# cell, as the isubd mapping closure. Reads the current solution and warp
# scales non-reactively — reactivity is the graph nodes' concern.
function _isubd_mapping(ds::FEData{dim}, cellmap::Vector{Int}) where {dim}
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    gips = [Ferrite.geometric_interpolation(typeof(c)) for c in cells]
    cellcoords = [Ferrite.getcoordinates(grid, i) for i in 1:length(cells)]
    warps = [(FieldEvaluator(ds.dh, _resolve_name(ds, fname[])), scale)
             for (fname, scale) in ds.deformation]
    u_obs = ds.u
    function mapping(base_id::Int, ξ)
        cell_id = cellmap[base_id]
        x = geometric_map(gips[cell_id], cellcoords[cell_id], ξ)
        for (ev, scale) in warps
            d = evaluate_at(ev, cell_id, cellcoords[cell_id], ξ, u_obs[])
            d === nothing && continue
            x += Float64(scale[]) * d
        end
        return GeometryBasics.Point{dim,Float32}(x...)
    end
    return mapping
end

function _isubd_base(ds::FEData)
    corners, cellmap = _isubd_base_triangles(ds)
    return IsubdBase(corners, _isubd_mapping(ds, cellmap)), cellmap
end

# Per-vertex field values at the decoded reference coordinates: the adaptive
# counterpart of transfer_solution, evaluating at the sub-triangle vertices
# instead of the static tessellation vertices.
function _transfer_at!(out::Vector{Float32}, ev::FieldEvaluator, ds::FEData, cellmap::Vector{Int},
                       keys::Vector{UInt64}, ξs, u::AbstractVector; reduce::Bool)
    grid = Ferrite.get_grid(ds.dh)
    resize!(out, length(ξs))
    local_coords = Ferrite.getcoordinates(grid, 1)
    lastcell = 0
    for (i, k) in enumerate(keys)
        cell_id = cellmap[key_base(k)]
        if cell_id != lastcell
            Ferrite.getcoordinates!(local_coords, grid, cell_id)
            lastcell = cell_id
        end
        for j in (3i - 2):(3i)
            val = evaluate_at(ev, cell_id, local_coords, ξs[j], u)
            out[j] = val === nothing ? NaN32 :
                     reduce ? Float32(LinearAlgebra.norm(val)) : Float32(val[1])
        end
    end
    return out
end

# The adaptive branch of solutionplot's plot!. The dataset, its visibility and
# the color resolution are read eagerly (as everywhere in the recipes); the
# camera, px_target, max_depth and the solution drive the graph.
function _adaptive_solutionplot!(SP, ds::FEData{dim}) where {dim}
    graph = SP.attributes
    cam = Makie.camera(Makie.parent_scene(SP))
    ComputePipeline.add_input!(graph, :subd_projectionview, cam.projectionview)
    ComputePipeline.add_input!(graph, :subd_eyeposition, cam.eyeposition)
    ComputePipeline.add_input!(graph, :subd_resolution, cam.resolution)
    ComputePipeline.add_input!(graph, :subd_u, ds.u)

    base, cellmap = _isubd_base(ds)
    state = (keys=root_keys(base), scratch=UInt64[], prev=UInt64[])

    # geometry (and therefore the LoD) depends on the solution only through warps
    keyinputs = [:subd_projectionview, :subd_eyeposition, :subd_resolution, :px_target, :max_depth]
    isempty(ds.deformation) || push!(keyinputs, :subd_u)
    ComputePipeline.register_computation!(graph, keyinputs, [:subd_keys]) do inputs, changed, cached
        lod = ScreenSpaceLoD(inputs.subd_projectionview,
                             Tuple(inputs.subd_eyeposition),
                             (Float64(inputs.subd_resolution[1]), Float64(inputs.subd_resolution[2])),
                             Float64(inputs.px_target))
        refine_keys!(state.keys, state.scratch, base, lod; max_depth=Int(inputs.max_depth))
        # a camera move that does not change the key set leaves the outputs
        # clean, so decode and color transfer are skipped entirely
        cached !== nothing && state.keys == state.prev && return nothing
        copy!(state.prev, state.keys)
        return (state.keys,)
    end

    mesh_buf = IsubdMesh(base)
    ComputePipeline.register_computation!(graph, [:subd_keys], [:subd_positions, :subd_ξ, :subd_faces]) do inputs, changed, cached
        decode_keys!(mesh_buf, inputs.subd_keys, base)
        # positions/faces enter Makie's mesh internals, which may retain them:
        # hand over fresh arrays. The reference coordinates stay in-graph
        # (only the color node below reads them), so the buffer is shared.
        faces = [GeometryBasics.GLTriangleFace(f[1], f[2], f[3]) for f in mesh_buf.faces]
        return (copy(mesh_buf.positions), mesh_buf.refcoords, faces)
    end

    colorval = SP.color[]
    if colorval isa Symbol && (colorval === :default || _data_association(ds, colorval) !== :none)
        fname = _resolve_name(ds, colorval)
        fname in Ferrite.getfieldnames(ds.dh) ||
            error("adaptive solutionplot colors by evaluating a dof field at the refined vertices; " *
                  ":$fname is a registered data array on the static tessellation and cannot be resampled. " *
                  "Color by a dof field or a plain color, or use adaptive=false.")
        ev = FieldEvaluator(ds.dh, fname)
        reduce = colorval === :default && ev.ncomps > 1
        reduce || ev.ncomps == 1 ||
            error("field :$fname has $(ev.ncomps) components; adaptive coloring needs a scalar " *
                  "(or :default, which reduces to the magnitude)")
        colors = Float32[]
        ComputePipeline.register_computation!(graph, [:subd_keys, :subd_ξ, :subd_u], [:subd_color]) do inputs, changed, cached
            _transfer_at!(colors, ev, ds, cellmap, inputs.subd_keys, inputs.subd_ξ, inputs.subd_u; reduce)
            return (copy(colors),)
        end
        colornode = SP.subd_color
    else
        colornode = SP.color   # plain color, converted by the mesh child
    end

    # plain mesh child from graph nodes; the ShaderAbstractions.Buffer path of
    # `_mesh!` is tied to the static tessellation and does not apply here
    return Makie.mesh!(SP, SP.attributes, SP.subd_positions, SP.subd_faces, color=colornode)
end
