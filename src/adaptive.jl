# View-adaptive tessellation: the glue between FEData and the isubd core
# (src/isubd.jl), and the compute-graph wiring for `solutionplot(...;
# adaptive=true)` (#161).
#
# Refinement is driven by two interpolation-error estimators, both instances
# of DeviationLoD along a triangle's split edge, refining when *either* asks:
#   - geometry error: the exact geometry (dofhandler interpolation + warps)
#     vs the flat triangle, relative to the grid's bounding-box diagonal;
#   - solution error: the exact field polynomial vs the linear vertex-color
#     interpolation, relative to the field's value span.
# The camera is deliberately not an input: keys recompute only when the
# solution or a tolerance changes. An optional screen-space criterion
# (px_target) can be added on top, which wires the camera in.
#
# Graph topology (all nodes prefixed subd_ to avoid clashes):
#
#   ds_u, geometry_tol, solution_tol, max_depth ──> subd_keys ─┐
#   [camera, px_target — only when px_target is set] ──┘       ├─> subd_positions, subd_ξ,
#   ds_u ──────────────────────────────────────────────────────┘   subd_faces, subd_color
#
# The persistent LEB key buffer lives in plot-local state captured by the
# subd_keys computation; the node returns `nothing` when an update does not
# change the key set, so the decode is skipped. Geometry, faces and colors
# are emitted by ONE edge, so they can never be resolved against different
# key sets (mismatched buffers rendered as dropped triangles).

# One dof field, evaluated at arbitrary (cell, ξ) — the continuous-geometry
# and per-vertex-colour workhorse, and the hot path of everything adaptive
# (the error estimators query it far more often than the rendering does).
#
# The value is summed straight from the reference shape functions, the way
# `geometric_map` does, rather than through `PointValues`: for H1 (Lagrange)
# fields the value mapping is the identity, so `reinit!`'s machinery — which
# evaluates every shape function into a buffer for `function_value` to read
# back — is pure overhead, worth a measured 3x. Fields whose values are
# mapped (Piola, i.e. H(curl)/H(div)) would need the real thing; they are
# unsupported here as in the rest of the package (#151).
struct FieldEvaluator{IP,PF}
    ips::Vector{IP}                   # one interpolation per subdofhandler with the field
    sdh_of_cell::Vector{Int}          # cell -> index into ips, 0 when the field is absent
    celldofs_field::Vector{Vector{Int}}  # cell -> global dofs of the field, empty when absent
    ncomps::Int
    # Monomial coefficients per cell, when the interpolation admits them (see
    # polyeval.jl). `prepare!` fills them for the cells that will be sampled;
    # `nothing` means every evaluation sums shape functions instead.
    poly::PF
end

function FieldEvaluator(dh::Ferrite.DofHandler, field::Symbol)
    sdhs = getsubdofhandlers(dh, field)
    isempty(sdhs) && error("field :$field not found in the DofHandler")
    ip_field = Ferrite.getfieldinterpolation(first(sdhs), field)
    ξ0 = Ferrite.Vec(ntuple(d -> 0.0, Ferrite.getrefdim(ip_field)))
    ncomps = length(Ferrite.reference_shape_value(ip_field, ξ0, 1))
    ips = [Ferrite.getfieldinterpolation(sdh, field) for sdh in sdhs]
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
    # only a single interpolation can share one coefficient layout
    basis = length(ips) == 1 ? PolyBasis(only(ips)) : nothing
    T = ncomps == 1 ? Float64 : Tensors.Vec{ncomps,Float64}
    poly = basis === nothing ? nothing : PolyField(basis, ncells, T)
    return FieldEvaluator(ips, sdh_of_cell, celldofs_field, ncomps, poly)
end

# Fill the coefficients of the cells that will be evaluated. Cheap (one small
# matvec per cell) next to the sampling that follows, and skipped entirely
# when the field has no polynomial form.
prepare!(::FieldEvaluator{IP,Nothing}, cells, u::AbstractVector) where {IP} = nothing
function prepare!(ev::FieldEvaluator, cells, u::AbstractVector)
    pf = ev.poly
    nodal = _nodal_buffer(pf)
    for cell in cells
        ev.sdh_of_cell[cell] == 0 && continue
        dofs = ev.celldofs_field[cell]
        _gather_nodal!(nodal, dofs, u, ev.ncomps)
        refresh_cell!(pf, cell, nodal)
    end
    return nothing
end

_nodal_buffer(pf::PolyField{refdim,N,P,T}) where {refdim,N,P,T} = Vector{T}(undef, N)

function _gather_nodal!(nodal::Vector{Float64}, dofs, u, ::Int)
    @inbounds for k in eachindex(nodal)
        nodal[k] = u[dofs[k]]
    end
    return nodal
end
function _gather_nodal!(nodal::Vector{Tensors.Vec{vdim,Float64}}, dofs, u, ::Int) where {vdim}
    @inbounds for k in eachindex(nodal)
        o = (k - 1) * vdim
        nodal[k] = Tensors.Vec{vdim}(ntuple(c -> u[dofs[o + c]], vdim))
    end
    return nodal
end

# Evaluate at one reference point of one cell; `nothing` outside the field's
# subdomain.
@inline function evaluate_at(ev::FieldEvaluator, cell_idx::Int, ξ, u::AbstractVector)
    si = ev.sdh_of_cell[cell_idx]
    si == 0 && return nothing
    pf = ev.poly
    if pf !== nothing && pf.filled[cell_idx]
        return evaluate(pf, cell_idx, ξ)
    end
    return _sum_shape_values(ev.ips[si], ev.celldofs_field[cell_idx], ξ, u)
end

# function barrier: the interpolation is only abstractly typed in the vector
# above, and this loop must specialize on it
function _sum_shape_values(ip, dofs::Vector{Int}, ξ, u::AbstractVector)
    val = Ferrite.reference_shape_value(ip, ξ, 1) * u[dofs[1]]
    @inbounds for i in 2:length(dofs)
        val += Ferrite.reference_shape_value(ip, ξ, i) * u[dofs[i]]
    end
    return val
end

# A vectorized interpolation's basis functions are the scalar ones repeated
# once per component, and asking for each of them re-evaluates the scalar
# underneath — so a displacement field costs `vdim` times what it needs to.
# Evaluate the scalar basis once instead and gather the components (Ferrite
# numbers them consecutively per scalar base function).
function _sum_shape_values(ip::Ferrite.VectorizedInterpolation{vdim}, dofs::Vector{Int},
                           ξ, u::AbstractVector) where {vdim}
    sip = ip.ip
    val = zero(Tensors.Vec{vdim,Float64})
    @inbounds for i in 1:Ferrite.getnbasefunctions(sip)
        N = Ferrite.reference_shape_value(sip, ξ, i)
        o = (i - 1) * vdim
        val += N * Tensors.Vec{vdim}(ntuple(c -> u[dofs[o + c]], vdim))
    end
    return val
end

# Global identity of every tessellation vertex of a cell, so that shared base
# edges can be matched *exactly* — by integer id, never by comparing floating
# point coordinates. A vertex coinciding with one of the cell's geometric nodes
# takes that node's global id (which is how two cells recognize their shared
# edge, and how two facets of one cell recognize theirs); anything else is
# cell-interior (a fan centre) and gets a fresh negative id that nobody else
# can collide with.
function _vertex_gids(cell, coords, counter::Base.RefValue{Int})
    refcoords = Ferrite.reference_coordinates(Ferrite.geometric_interpolation(typeof(cell)))
    nodes = cell.nodes
    return map(coords) do ξ
        j = findfirst(rc -> isapprox(rc, ξ; atol=1e-12), refcoords)
        j === nothing ? (counter[] -= 1) : Int(nodes[j])
    end
end

# The isubd base domain of a dataset: one LEB-ordered corner triple per base
# triangle, the triangle -> cell map, and a global vertex id per corner (for
# the adjacency table). Static per plot, like the visibility mask.
#
# In 2D the base is built as a fan from each cell's centre over its *element
# edges*, so a base triangle's split edge is always an element edge. That is
# what makes the base compatible in the sense conforming refinement needs:
# element edges are shared by exactly two base triangles which both treat them
# as their split edge (a "diamond"), while the fan's interior edges are legs on
# both sides. For quadrilaterals this reproduces the existing centre fan; for
# triangles it replaces the single base triangle, whose split edge (the longest
# reference edge) would generally meet a neighbour's leg.
#
# In 3D the same construction applies per *surface facet*: a facet is on the
# surface when the cell across it is missing or not part of the body (see
# `FEData.solid`), and it is fanned from its centre, so a base triangle's split
# edge is again an element edge — now shared by exactly the two surface facets
# meeting there. Note this draws strictly less than the static path, which
# tessellates every facet of every visible cell including the ones buried
# inside the body.
# Is this facet of the cell part of the drawn surface? It is when nothing sits
# across it, or what sits there is not part of the body (an interior cell of a
# solid mesh is, a cell a CrinkleClip removed is not). Without topology every
# facet is taken, which is the 2D case and the pre-topology fallback.
function _is_surface_facet(ds::FEData, cell_id::Int, facet::Int)
    ds.topology === nothing && return true
    neighbours = ds.topology.face_face_neighbor[cell_id, facet]
    isempty(neighbours) && return true
    return !ds.solid[first(neighbours)[1]]
end

# One centre fan of a convex polygon given by its rim vertices in order:
# (v[i+1], centre, v[i]) keeps the rim's winding while making the rim edge the
# triangle's split edge, which is what the conforming refinement needs.
function _push_fan!(corners, cornergids, cellmap, cell_id, rim, rimgids, counter)
    centre = sum(rim) / length(rim)
    centre_gid = (counter[] -= 1)
    for i in eachindex(rim)
        j = mod1(i + 1, length(rim))
        push!(corners, (rim[j], centre, rim[i]))
        push!(cornergids, (rimgids[j], centre_gid, rimgids[i]))
        push!(cellmap, cell_id)
    end
    return nothing
end

function _isubd_base_triangles(ds::FEData)
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    refdim = Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(first(cells))))
    corners = NTuple{3,Ferrite.Vec{refdim,Float64}}[]
    cornergids = NTuple{3,Int}[]
    cellmap = Int[]
    counter = Ref(0)
    for (cell_id, cell) in enumerate(cells)
        ds.visible[cell_id] || continue
        Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(cell))) == refdim ||
            error("adaptive tessellation requires a single reference dimension across the grid")
        refshape = getrefshape(cell)
        tess = reference_tessellation(refshape)
        isempty(tess.triangles) && continue     # e.g. line cells carry no surface
        if refdim == 2 && !isempty(tess.edges)
            gids = _vertex_gids(cell, tess.coords, counter)
            # the element edges, in whatever order the tessellation lists them
            rim = unique(Iterators.flatten(tess.edges))
            centre = sum(tess.coords[i] for i in rim) / length(rim)
            centre_gid = (counter[] -= 1)
            for (a, b) in tess.edges
                # (b, centre, a) keeps the fan's winding while making the
                # element edge (b, a) the triangle's split edge
                push!(corners, (tess.coords[b], centre, tess.coords[a]))
                push!(cornergids, (gids[b], centre_gid, gids[a]))
                push!(cellmap, cell_id)
            end
        elseif refdim == 3
            vcoords = Ferrite.reference_coordinates(Ferrite.Lagrange{refshape,1}())
            vgids = _vertex_gids(cell, vcoords, counter)
            for (f, face) in enumerate(Ferrite.reference_faces(refshape))
                _is_surface_facet(ds, cell_id, f) || continue
                _push_fan!(corners, cornergids, cellmap, cell_id,
                           [vcoords[v] for v in face], [vgids[v] for v in face], counter)
            end
        else
            gids = _vertex_gids(cell, tess.coords, counter)
            for tri in tess.triangles
                c = leb_order((tess.coords[tri[1]], tess.coords[tri[2]], tess.coords[tri[3]]))
                # recover the permutation leb_order applied, to keep the ids aligned
                perm = map(x -> findfirst(i -> tess.coords[i] === x, tri), c)
                push!(corners, c)
                push!(cornergids, (gids[tri[perm[1]]], gids[tri[perm[2]]], gids[tri[perm[3]]]))
                push!(cellmap, cell_id)
            end
        end
    end
    isempty(corners) && error("adaptive tessellation found no visible cells to tessellate")
    return corners, cornergids, cellmap
end

# Pair up base triangles along shared edges. Only exact 2-triangle matches
# become neighbours, and only when both sides agree on the edge's role: a split
# edge may pair with a split edge and a leg with a leg, never across. Anything
# else stays a boundary, which costs conformity along that edge and nothing
# else (see `key_neighbour`).
function _base_adjacency(cornergids::Vector{NTuple{3,Int}})
    edges = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int,Bool}}}()
    for (t, g) in enumerate(cornergids)
        for (e, (i, j)) in enumerate(((1, 3), (1, 2), (2, 3)))   # EDGE_S, EDGE_L, EDGE_R
            a, b = g[i], g[j]
            forward = a < b
            push!(get!(Vector{Tuple{Int,Int,Bool}}, edges, forward ? (a, b) : (b, a)),
                  (t, e, forward))
        end
    end
    adjacency = [(NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR) for _ in cornergids]
    incompatible = 0
    for entries in values(edges)
        length(entries) == 2 || continue                     # boundary, or a non-manifold edge
        (t1, e1, f1), (t2, e2, f2) = entries
        if (e1 == EDGE_S) != (e2 == EDGE_S)
            incompatible += 1
            continue
        end
        reversed = f1 != f2
        adjacency[t1] = Base.setindex(adjacency[t1], (t2, e2, reversed), e1)
        adjacency[t2] = Base.setindex(adjacency[t2], (t1, e1, reversed), e2)
    end
    return adjacency, incompatible
end

# Continuous geometry x(ξ) [+ Σ scaleᵢ · fieldᵢ(ξ) for upstream warps] of one
# cell, as the isubd mapping closure. Reads the current solution and warp
# scales non-reactively — reactivity is the graph nodes' concern. Full Float64:
# the geometry-error estimator measures deviations far below Float32 noise
# (only the positions handed to Makie get truncated, in the decode node).
function _isubd_mapping(ds::FEData{dim}, cellmap::Vector{Int}, cellcoords) where {dim}
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    gips = [Ferrite.geometric_interpolation(typeof(c)) for c in cells]
    warps = [(FieldEvaluator(ds.dh, _resolve_name(ds, fname[])), scale)
             for (fname, scale) in ds.deformation]
    u_obs = ds.u
    # The geometry's own coefficients never change — the node coordinates are
    # fixed — so they are computed once here rather than per update.
    geo = _geometry_poly(gips, cellcoords, cellmap, Val(dim))
    function mapping(base_id::Int, ξ)
        cell_id = cellmap[base_id]
        x = (geo !== nothing && geo.filled[cell_id]) ? evaluate(geo, cell_id, ξ) :
            geometric_map(gips[cell_id], cellcoords[cell_id], ξ)
        for (ev, scale) in warps
            d = evaluate_at(ev, cell_id, ξ, u_obs[])
            d === nothing && continue
            x += Float64(scale[]) * d
        end
        return x
    end
    return mapping, warps
end

function _geometry_poly(gips, cellcoords, cellmap, ::Val{dim}) where {dim}
    basis = PolyBasis(first(gips))
    basis === nothing && return nothing
    poly = PolyField(basis, length(gips), Tensors.Vec{dim,Float64})
    for cell in unique(cellmap)
        # a differently interpolated cell keeps the shape-function path
        PolyBasis(gips[cell]) === nothing && continue
        length(cellcoords[cell]) == size(poly.coeffs, 1) || continue
        refresh_cell!(poly, cell, cellcoords[cell])
    end
    return poly
end

# Positions leave the graph as Float32 points for rendering; the pipeline
# itself stays in Float64 (the geometry estimator measures deviations far
# below Float32 noise). Filled into a buffer that is handed out as is, see
# `_adaptive_solutionplot!`.
function _render_positions!(out::Vector{GeometryBasics.Point{dim,Float32}}, ps) where {dim}
    resize!(out, length(ps))
    @inbounds for i in eachindex(ps)
        out[i] = GeometryBasics.Point{dim,Float32}(ps[i]...)
    end
    return out
end

function _isubd_base(ds::FEData)
    corners, cornergids, cellmap = _isubd_base_triangles(ds)
    adjacency, _ = _base_adjacency(cornergids)
    grid = Ferrite.get_grid(ds.dh)
    cellcoords = [Ferrite.getcoordinates(grid, i) for i in 1:Ferrite.getncells(grid)]
    mapping, warps = _isubd_mapping(ds, cellmap, cellcoords)
    base = IsubdBase(corners, mapping, adjacency)
    return base, cellmap, cellcoords, warps
end

# Per-vertex field values at the decoded reference coordinates: the adaptive
# counterpart of transfer_solution, evaluating at the sub-triangle vertices
# instead of the static tessellation vertices.
function _transfer_at!(out::Vector{Float32}, ev::FieldEvaluator, cellmap::Vector{Int},
                       mesh, u::AbstractVector; reduce::Bool)
    resize!(out, length(mesh.refcoords))
    @inbounds for v in eachindex(mesh.refcoords)
        val = evaluate_at(ev, cellmap[mesh.vertex_base[v]], mesh.refcoords[v], u)
        out[v] = val === nothing ? NaN32 :
                 reduce ? Float32(LinearAlgebra.norm(val)) : Float32(val[1])
    end
    return out
end

# Refresh the per-cell coefficients of everything sampled from the solution,
# before anything reads them. Coefficients that do not match the `u` in hand
# would evaluate silently wrong values, so this runs at the top of every graph
# node that samples, not once per update somewhere central.
function _prepare_fields!(warps, ev, cells, u::AbstractVector)
    for (wev, _) in warps
        prepare!(wev, cells, u)
    end
    ev === nothing || prepare!(ev, cells, u)
    return nothing
end

# The scalar drawn as color, as a (base_id, ξ) probe for the solution-error
# estimator; `u` is bound per key-node invocation.
function _scalar_probe(ev::FieldEvaluator, cellmap::Vector{Int}, cellcoords, u::AbstractVector; reduce::Bool)
    function probe(base_id::Int, ξ)
        cell_id = cellmap[base_id]
        val = evaluate_at(ev, cell_id, ξ, u)
        val === nothing && return 0.0
        return reduce ? Float64(LinearAlgebra.norm(val)) : Float64(val[1])
    end
    return probe
end

# Value span of the color over the base triangle corners — the reference scale
# for the relative solution tolerance.
function _base_span(probe, base::IsubdBase)
    lo, hi = Inf, -Inf
    for b in 1:length(base.corners), c in base.corners[b]
        v = probe(b, c)
        isfinite(v) || continue
        lo, hi = min(lo, v), max(hi, v)
    end
    return hi > lo ? hi - lo : 0.0
end

function _grid_diagonal(grid)
    nodes = Ferrite.getnodes(grid)
    lo = hi = Ferrite.get_node_coordinate(first(nodes))
    for n in nodes
        x = Ferrite.get_node_coordinate(n)
        lo, hi = min.(lo, x), max.(hi, x)
    end
    return LinearAlgebra.norm(hi - lo)
end

# The wireframe of an adaptively refined mesh: the element edges, subdivided
# exactly as the surface subdivided them.
#
# An element edge *is* a base triangle's split edge (that is what the fan
# construction guarantees), and a point lies on it precisely when its
# barycentric weight for the apex vanishes — which `key_xform` gives exactly,
# since the bisection weights are dyadic. So a leaf contributes a wireframe
# segment when two of its corners have zero apex weight, and conformity
# guarantees the two triangles sharing an element edge subdivide it
# identically. Each edge is therefore emitted once, by the base triangle with
# the smaller id of the pair.
function _element_edge_segments!(out::Vector{PT}, keys::Vector{UInt64}, base::IsubdBase) where {PT}
    empty!(out)
    for k in keys
        b = key_base(k)
        nb = base.adjacency[b][EDGE_S]
        (nb[1] != 0 && nb[1] < b) && continue   # the partner draws this one
        X = key_xform(k)
        on = (X[2, 1] == 0.0, X[2, 2] == 0.0, X[2, 3] == 0.0)
        (count(on) == 2) || continue            # a leaf touches the edge with at most one of its own
        c = key_corners(base, k)
        for (i, j) in ((1, 2), (2, 3), (3, 1))
            (on[i] && on[j]) || continue
            push!(out, PT(base.mapping(b, c[i])...), PT(base.mapping(b, c[j])...))
        end
    end
    return out
end

# The adaptive branch of meshplot's plot!: same refinement machinery as the
# surface, but only the geometry criterion — a wireframe has no field to
# resolve, only a curve to follow.
function _adaptive_wireframe!(WF, ds::FEData{dim}) where {dim}
    graph = WF.attributes
    ComputePipeline.add_input!(graph, :subd_u, ds.u)
    base, cellmap, _, warps = _isubd_base(ds)
    conforming = WF.conforming[] && is_conformable(base)
    used_cells = unique(cellmap)
    state = (keys=root_keys(base), prev=UInt64[], scratch=UInt64[])
    diag = _grid_diagonal(Ferrite.get_grid(ds.dh))
    ComputePipeline.register_computation!(graph, [:subd_u, :geometry_tol, :max_depth],
                                          [:subd_keys]) do inputs, changed, cached
        _prepare_fields!(warps, nothing, used_cells, inputs.subd_u)
        lod = DeviationLoD(base.mapping, Float64(inputs.geometry_tol) * diag)
        refine_keys!(state.keys, state.scratch, base, lod;
                     max_depth=Int(inputs.max_depth), conforming)
        cached !== nothing && state.keys == state.prev && return nothing
        copy!(state.prev, state.keys)
        return (state.keys,)
    end
    # Makie draws 2D/3D points; pad 1D grids with a zero y-coordinate
    segments = GeometryBasics.Point{max(dim, 2),Float32}[]
    Makie.map!(graph, [:subd_keys, :subd_u], :edge_lines) do keys, u
        _prepare_fields!(warps, nothing, used_cells, u)
        return _element_edge_segments!(segments, keys, base)
    end
    return nothing
end

# The adaptive branch of solutionplot's plot!. The dataset, its visibility and
# the color resolution are read eagerly (as everywhere in the recipes); the
# solution and the tolerances drive the graph.
function _adaptive_solutionplot!(SP, ds::FEData{dim}) where {dim}
    graph = SP.attributes
    ComputePipeline.add_input!(graph, :subd_u, ds.u)

    base, cellmap, cellcoords, warps = _isubd_base(ds)
    conforming = SP.conforming[] && is_conformable(base)
    # cells the plot samples; their coefficients are refreshed per update
    used_cells = unique(cellmap)
    state = (keys=root_keys(base), scratch=UInt64[], prev=UInt64[])

    # color resolution (eager): a dof field evaluated per vertex, or a plain color
    colorval = SP.color[]
    ev = nothing
    reduce = false
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
    end

    # optional screen-space criterion: only then does the camera enter the graph
    px_target = SP.px_target[]
    keyinputs = [:subd_u, :geometry_tol, :solution_tol, :max_depth]
    if px_target !== nothing
        cam = Makie.camera(Makie.parent_scene(SP))
        ComputePipeline.add_input!(graph, :subd_projectionview, cam.projectionview)
        ComputePipeline.add_input!(graph, :subd_eyeposition, cam.eyeposition)
        ComputePipeline.add_input!(graph, :subd_resolution, cam.resolution)
        append!(keyinputs, [:subd_projectionview, :subd_eyeposition, :subd_resolution, :px_target])
    end

    diag = _grid_diagonal(Ferrite.get_grid(ds.dh))
    span = Ref(NaN)
    ComputePipeline.register_computation!(graph, keyinputs, [:subd_keys]) do inputs, changed, cached
        _prepare_fields!(warps, ev, used_cells, inputs.subd_u)
        geo = DeviationLoD(base.mapping, Float64(inputs.geometry_tol) * diag)
        lods = (geo,)
        if ev !== nothing
            probe = _scalar_probe(ev, cellmap, cellcoords, inputs.subd_u; reduce)
            if isnan(span[]) || changed.subd_u
                span[] = _base_span(probe, base)
            end
            # a (near-)constant field never asks for refinement
            tol = max(Float64(inputs.solution_tol) * span[], 1e-12)
            lods = (lods..., DeviationLoD(probe, tol))
        end
        if px_target !== nothing
            lods = (lods..., ScreenSpaceLoD(inputs.subd_projectionview,
                                            Tuple(inputs.subd_eyeposition),
                                            (Float64(inputs.subd_resolution[1]), Float64(inputs.subd_resolution[2])),
                                            Float64(inputs.px_target)))
        end
        refine_keys!(state.keys, state.scratch, base, CombinedLoD(lods);
                     max_depth=Int(inputs.max_depth), conforming)
        # an update that does not change the key set leaves the output clean,
        # so the decode below is skipped entirely
        cached !== nothing && state.keys == state.prev && return nothing
        copy!(state.prev, state.keys)
        return (state.keys,)
    end

    # geometry, connectivity and colors from ONE edge: they can never be
    # resolved against different key sets. positions/faces enter Makie's mesh
    # internals, which may retain them: hand over fresh arrays. The reference
    # coordinates stay in-graph, so the buffer is shared.
    # Connectivity and values are separate nodes: the vertex layout depends
    # only on the key set, so a solution change on an unchanged mesh skips the
    # vertex-sharing lookups and re-evaluates one point per vertex — the same
    # work the static path does. Values still cannot disagree with the mesh
    # they belong to, because they are computed *from* the layout node's
    # output rather than beside it.
    #
    # Vertices are shared within a cell (where the drawn field is continuous)
    # and duplicated across cells, so element-boundary jumps survive.
    # Every node hands out the same buffer it filled, rather than a copy: the
    # pipeline compares an output against its previous value to decide what is
    # stale, and an array that *is* the previous one (same pointer) counts as
    # changed — so in-place reuse propagates exactly like a fresh array, while
    # a mesh that keeps its size stops allocating per update altogether.
    mesh_buf = IsubdMesh(base)
    faces_out = GeometryBasics.GLTriangleFace[]
    positions_out = GeometryBasics.Point{dim,Float32}[]
    colors_out = Float32[]
    ComputePipeline.register_computation!(graph, [:subd_keys], [:subd_ξ, :subd_faces]) do inputs, changed, cached
        decode_topology!(mesh_buf, inputs.subd_keys, base; groups=cellmap)
        resize!(faces_out, length(mesh_buf.faces))
        @inbounds for i in eachindex(mesh_buf.faces)
            f = mesh_buf.faces[i]
            faces_out[i] = GeometryBasics.GLTriangleFace(f[1], f[2], f[3])
        end
        return (mesh_buf.refcoords, faces_out)
    end
    Makie.map!(graph, [:subd_ξ, :subd_u], :subd_positions) do _ξ, uu
        _prepare_fields!(warps, ev, used_cells, uu)
        decode_positions!(mesh_buf, base)
        return _render_positions!(positions_out, mesh_buf.positions)
    end
    if ev === nothing
        colornode = SP.color   # plain color, converted by the mesh child
    else
        Makie.map!(graph, [:subd_ξ, :subd_u], :subd_color) do _ξ, u
            _prepare_fields!(warps, ev, used_cells, u)
            return _transfer_at!(colors_out, ev, cellmap, mesh_buf, u; reduce)
        end
        colornode = SP.subd_color
    end

    # plain mesh child from graph nodes; the ShaderAbstractions.Buffer path of
    # `_mesh!` is tied to the static tessellation and does not apply here
    return Makie.mesh!(SP, SP.attributes, SP.subd_positions, SP.subd_faces, color=colornode)
end
