# Layer 2b: the reactive filter graph.
#
# A filter maps FEData -> FEData. Filters are callable, so pipelines compose
# with |>:  ds |> WarpByVector(:u, 2.0) |> Gradient(:u) |> VonMises().
# Every apply returns a new FEData sharing the source solution observable, so
# pipelines stay live under FerriteViz.update! and fork without clobbering each
# other. Geometry-rebuilding filters (Refine, FirstOrderRefinement) rebuild
# from the base geometry — apply WarpByVector after them (CrinkleClip and
# Gradient share their input's coordinates, so warps survive them).
#
# Derived datasets are wired to the source observable and are kept alive by it;
# there is no explicit disposal — drop the source to release a pipeline.

"""
    apply(f::AbstractFilter, ds::FEData) -> FEData

Apply a filter to a dataset. Filters are callable, so `ds |> f` is equivalent.
"""
function apply end

(f::AbstractFilter)(ds::FEData) = apply(f, ds)

# Rebind dh/u on identical geometry (same grid ⇒ identical tessellation).
function _rebind(ds::FEData{dim}, dh, u::Makie.Observable;
                 point_data=Dict{Symbol,Makie.Observable}(), cell_data=Dict{Symbol,Makie.Observable}()) where {dim}
    return FEData{dim,typeof(dh),eltype(u[]),typeof(ds.topology),typeof(ds.source_u),typeof(ds.mesh),eltype(ds.all_triangles)}(
        dh, u, ds.source_u, ds.topology, ds.visible, ds.gridnodes, ds.coords, ds.coords_buffer,
        ds.all_triangles, ds.vis_triangles, ds.triangle_cell_map, ds.cell_triangle_offsets,
        ds.cell_vertex_offsets, ds.reference_coords, ds.mesh, point_data, cell_data)
end

################
# WarpByVector #
################

"""
    WarpByVector(field=:default, scale=1.0)

Filter displacing the geometry (tessellation vertices and grid nodes) by
`scale` times the vector-valued `field`. `scale` may be a number or an
`Observable` (e.g. driven by a slider). Warps compose: the displacement is
added to the input dataset's current coordinates.

`field` may name any vector-valued point-data array; the original grid nodes
(drawn by `meshplot`) can only be displaced when it is a dof-backed field of
the dof handler and stay put otherwise.
"""
struct WarpByVector{F<:Union{Symbol,Makie.Observable{Symbol}},S} <: AbstractFilter
    field::F
    scale::S
end
WarpByVector(field=:default) = WarpByVector(field, 1.0)

function apply(w::WarpByVector, ds::FEData{dim}) where {dim}
    scale = w.scale isa Makie.Observable ? w.scale : Makie.Observable(w.scale)
    fname = w.field isa Makie.Observable ? w.field : Makie.Observable(_resolve_name(ds, w.field))
    disp = _switching_point_data(ds, fname)
    size(disp[], 2) == dim || error("deformation field :$(fname[]) has $(size(disp[], 2)) components, expected $dim")
    coords = Makie.lift(ds.coords, disp, scale) do base, d, s
        size(d, 2) == dim || error("deformation field has $(size(d, 2)) components, expected $dim")
        [base[i] .+ Float32(s) .* GeometryBasics.Point{dim,Float32}(view(d, i, :)...) for i in eachindex(base)]
    end
    gridnodes = Makie.lift(ds.gridnodes, ds.u, scale, fname) do nodes, u, s, fn
        fn = _resolve_name(ds, fn)
        # named (non-dof) arrays live on the tessellation vertices only
        fn in Ferrite.getfieldnames(ds.dh) || return nodes
        vals = Ferrite.evaluate_at_grid_nodes(ds.dh, u, fn)
        [nodes[i] .+ Float32(s) .* GeometryBasics.Point{dim,Float32}(vals[i]...) for i in eachindex(nodes)]
    end
    coords_buffer = ShaderAbstractions.Buffer(coords)
    mesh = GeometryBasics.Mesh(coords_buffer, ds.vis_triangles)
    return FEData{dim,typeof(ds.dh),eltype(ds.u[]),typeof(ds.topology),typeof(ds.source_u),typeof(mesh),eltype(ds.all_triangles)}(
        ds.dh, ds.u, ds.source_u, ds.topology, ds.visible, gridnodes, coords, coords_buffer,
        ds.all_triangles, ds.vis_triangles, ds.triangle_cell_map, ds.cell_triangle_offsets,
        ds.cell_vertex_offsets, ds.reference_coords, mesh, copy(ds.point_data), copy(ds.cell_data))
end

# Register a listener for cleanup when `owner` (a plot) is deleted; without an
# owner the listener lives as long as the observed observable.
_register_listener!(::Nothing, obsfunc) = obsfunc
_register_listener!(owner, obsfunc) = (push!(owner.deregister_callbacks, obsfunc); obsfunc)

# Flatten an observable name -> Observable of the named point-data array,
# rewiring the inner listener when the name changes.
function _switching_point_data(ds::FEData, name_obs; owner=nothing)
    inner = point_data(ds, name_obs[])
    out = Makie.Observable(inner[])
    listener = Ref(_register_listener!(owner, Makie.on(v -> out[] = v, inner)))
    _register_listener!(owner, Makie.on(name_obs) do name
        Makie.Observables.off(listener[]) # double-off at plot deletion is harmless
        new_inner = point_data(ds, name)
        listener[] = _register_listener!(owner, Makie.on(v -> out[] = v, new_inner))
        out[] = new_inner[]
    end)
    return out
end

###############
# CrinkleClip #
###############

"""
    ClipPlane{T}(normal, distance_to_origin)

Clip plane described by its normal and distance to the coordinate origin, for
use as the decision function of [`CrinkleClip`](@ref): callable as
`plane(grid, cellid)`, returning whether the cell lies on the non-clipped side.
"""
struct ClipPlane{T}
    normal::Tensors.Vec{3,T}
    distance::T
end

function (plane::ClipPlane)(grid::Ferrite.AbstractGrid, cellid::Int)
    cell = getcells(grid, cellid)
    coords = Ferrite.get_node_coordinate.(Ferrite.getnodes(grid)[[cell.nodes...]])
    for coord ∈ coords
        if coord ⋅ plane.normal > plane.distance
            return false
        end
    end
    return true
end

"""
    CrinkleClip(decision)

Filter hiding the cells for which `decision(grid, cellid)` is false, revealing
the (crinkled) interior along the clip surface. `decision` is typically a
[`ClipPlane`](@ref).
"""
struct CrinkleClip{F} <: AbstractFilter
    decision::F
end

function apply(c::CrinkleClip, ds::FEData{3})
    ds.topology === nothing && error("CrinkleClip needs the dataset's topology; construct FEData with one")
    grid = Ferrite.get_grid(ds.dh)
    visible = copy(ds.visible)
    for cell_id in 1:Ferrite.getncells(grid)
        if c.decision(grid, cell_id)
            cell_neighbors = Ferrite.getneighborhood(ds.topology, grid, Ferrite.CellIndex(cell_id))
            visible[cell_id] = !all(c.decision.((grid,), cell_neighbors)) || ds.visible[cell_id]
        else
            visible[cell_id] = false
        end
    end
    vis_triangles = ShaderAbstractions.Buffer(Makie.Observable(_visibility_triangles(ds.all_triangles, visible, ds.triangle_cell_map)))
    mesh = GeometryBasics.Mesh(ds.coords_buffer, vis_triangles)
    return FEData{3,typeof(ds.dh),eltype(ds.u[]),typeof(ds.topology),typeof(ds.source_u),typeof(mesh),eltype(ds.all_triangles)}(
        ds.dh, ds.u, ds.source_u, ds.topology, visible, ds.gridnodes, ds.coords, ds.coords_buffer,
        ds.all_triangles, vis_triangles, ds.triangle_cell_map, ds.cell_triangle_offsets,
        ds.cell_vertex_offsets, ds.reference_coords, mesh,
        Dict{Symbol,Makie.Observable}(), copy(ds.cell_data))
end

##########
# Refine #
##########

"""
    Refine(n=1)

Filter subdividing every triangle into 4 (in reference space, orientation
preserving), `n` times. New vertices are mapped through the cell's geometric
interpolation, so both solution resolution and curved geometry improve.

!!! danger
    This filter has high RAM usage!
"""
struct Refine <: AbstractFilter
    n::Int
end
Refine() = Refine(1)

function apply(r::Refine, ds::FEData)
    out = ds
    for _ in 1:r.n
        out = _refine_once(out)
    end
    return out
end

# The reference dimension is the interpolation's, not the spatial one (padded
# reference-coordinate rows may be wider, e.g. surface cells in a 3D grid).
function _map_refined_vertex(ip_geo::Ferrite.ScalarInterpolation, node_coords, refc)
    ξ = Tensors.Vec{Ferrite.getrefdim(ip_geo)}(d -> refc[d])
    return geometric_map(ip_geo, node_coords, ξ)
end

function _refine_once(ds::FEData{dim}) where {dim}
    grid = Ferrite.get_grid(ds.dh)
    total_triangles = length(ds.all_triangles)
    refined_reference_coords = Matrix{Float64}(undef, 4 * 3 * total_triangles, dim)
    refined_physical_coords = Vector{GeometryBasics.Point{dim,Float32}}(undef, 4 * 3 * total_triangles)
    refined_triangle_cell_map = Vector{Int}(undef, 4 * total_triangles)
    refined_triangles = Matrix{Int}(undef, 4 * total_triangles, 3)

    for cell_id in 1:Ferrite.getncells(grid)
        ip_geo = Ferrite.geometric_interpolation(typeof(Ferrite.getcells(grid, cell_id)))
        node_coords = Ferrite.getcoordinates(grid, cell_id)
        for triangle_index in triangles_on_cell(ds, cell_id)
            tri = ds.all_triangles[triangle_index]
            v1 = ds.reference_coords[tri[1], :]
            v2 = ds.reference_coords[tri[2], :]
            v3 = ds.reference_coords[tri[3], :]
            m12 = (v1 + v2) / 2.0
            m23 = (v2 + v3) / 2.0
            m31 = (v3 + v1) / 2.0
            voff = 4 * 3 * (triangle_index - 1)
            # 4 sub-triangles, orientation preserving
            for (k, refc) in enumerate((v1, m12, m31, v2, m23, m12, v3, m31, m23, m12, m23, m31))
                refined_reference_coords[voff+k, :] = refc
                x = _map_refined_vertex(ip_geo, node_coords, refc)
                refined_physical_coords[voff+k] = GeometryBasics.Point{dim,Float32}(x...)
            end
            toff = 4 * (triangle_index - 1)
            refined_triangle_cell_map[(toff+1):(toff+4)] .= cell_id
            for t in 1:4
                refined_triangles[toff+t, :] = voff .+ (3 * (t - 1)) .+ (1:3)
            end
        end
    end

    all_triangles = convert(Vector{GeometryBasics.GLTriangleFace}, Makie.to_triangles(refined_triangles))
    vis_triangles = ShaderAbstractions.Buffer(Makie.Observable(_visibility_triangles(all_triangles, ds.visible, refined_triangle_cell_map)))
    coords = Makie.Observable(refined_physical_coords)
    coords_buffer = ShaderAbstractions.Buffer(coords)
    mesh = GeometryBasics.Mesh(coords_buffer, vis_triangles)
    return FEData{dim,typeof(ds.dh),eltype(ds.u[]),typeof(ds.topology),typeof(ds.source_u),typeof(mesh),eltype(all_triangles)}(
        ds.dh, ds.u, ds.source_u, ds.topology, ds.visible, ds.gridnodes, coords, coords_buffer,
        all_triangles, vis_triangles, refined_triangle_cell_map, ds.cell_triangle_offsets .* 4,
        ds.cell_triangle_offsets .* 12, refined_reference_coords, mesh,
        Dict{Symbol,Makie.Observable}(), copy(ds.cell_data))
end

########################
# FirstOrderRefinement #
########################

"""
    FirstOrderRefinement()

Filter replacing a high-order discretization by the first-order one spanned by
its nodes (see [`first_order_subcells`](@ref)), transferring the solution.
Requires a single subdofhandler with a single (scalar or vector) Lagrange field.
"""
struct FirstOrderRefinement <: AbstractFilter end

function apply(::FirstOrderRefinement, ds::FEData)
    dh = ds.dh
    _check_full_domain(dh, "FirstOrderRefinement")
    length(Ferrite.getfieldnames(dh)) == 1 || error("FirstOrderRefinement supports only a single field")
    dh_new, transfer = _first_order_discretization(dh)
    u_new = Makie.lift(transfer, ds.u)
    return FEData(dh_new, u_new; source_u=ds.source_u)
end

function _first_order_discretization(dh)
    sdh = dh.subdofhandlers[1]
    field_name = first(Ferrite.getfieldnames(sdh))
    grid = Ferrite.get_grid(dh)
    ip = Ferrite.getfieldinterpolation(sdh, field_name)
    vdim = Ferrite.n_components(dh, field_name)
    sip = ip isa VectorizedInterpolation ? ip.ip : ip
    subcells = first_order_subcells(sip)
    CT = linear_celltype(Ferrite.getrefshape(sip))

    ref_coords = Ferrite.reference_coordinates(sip)
    nodes_per_cell = length(ref_coords)
    ip_geo = Ferrite.geometric_interpolation(Ferrite.getcelltype(sdh))
    qr = Ferrite.QuadratureRule{Ferrite.getrefshape(sip)}(zeros(nodes_per_cell), ref_coords)
    cv = Ferrite.CellValues(qr, sip, ip_geo)

    # One new grid node per scalar basis function; with a single field the dofs
    # of node j are vdim*(j-1)+1 : vdim*j, which is what makes the node
    # identification below work.
    nnodes_new = Ferrite.ndofs(dh) ÷ vdim
    nodes = Vector{Ferrite.Node{Ferrite.getspatialdim(grid),Float64}}(undef, nnodes_new)
    cells = Vector{CT}(undef, Ferrite.getncells(grid) * length(subcells))
    nodeid(dofs_f, q) = div(dofs_f[vdim*(q-1)+1] - 1, vdim) + 1
    cellidx = 1
    for cell in Ferrite.CellIterator(sdh)
        Ferrite.reinit!(cv, cell)
        coords = Ferrite.getcoordinates(cell)
        dofs_f = Ferrite.celldofs(cell)[Ferrite.dof_range(sdh, field_name)]
        for q in 1:nodes_per_cell
            nodes[nodeid(dofs_f, q)] = Ferrite.Node(Ferrite.spatial_coordinate(cv, q, coords))
        end
        for sub in subcells
            cells[cellidx] = CT(map(k -> nodeid(dofs_f, k), sub))
            cellidx += 1
        end
    end

    grid_new = Ferrite.Grid(cells, nodes)
    dh_new = Ferrite.DofHandler(grid_new)
    lip = Ferrite.Lagrange{Ferrite.getrefshape(sip),1}()
    add!(dh_new, field_name, vdim > 1 ? lip^vdim : lip)
    close!(dh_new)

    rng = Ferrite.dof_range(dh_new.subdofhandlers[1], field_name)
    function transfer(u)
        u_new = zeros(eltype(u), Ferrite.ndofs(dh_new))
        cdofs = zeros(Int, Ferrite.ndofs_per_cell(dh_new))
        for cell_idx in 1:Ferrite.getncells(grid_new)
            Ferrite.celldofs!(cdofs, dh_new, cell_idx)
            dofs = @view cdofs[rng]
            for (k, node) in enumerate(Ferrite.getcells(grid_new, cell_idx).nodes), c in 1:vdim
                u_new[dofs[vdim*(k-1)+c]] = u[vdim*(node-1)+c]
            end
        end
        return u_new
    end
    return dh_new, transfer
end

############
# Gradient #
############

"""
    Gradient(field=:default; copy_fields=Symbol[])

Filter computing the piecewise discontinuous gradient of `field` (via
[`interpolate_gradient_field`](@ref)). The output dataset's field is named
`:gradient`; fields listed in `copy_fields` are carried along. Geometry
(including an upstream warp) is shared with the input.
"""
struct Gradient <: AbstractFilter
    field::Symbol
    copy_fields::Vector{Symbol}
end
Gradient(field::Symbol=:default; copy_fields::Vector{Symbol}=Symbol[]) = Gradient(field, copy_fields)

function apply(g::Gradient, ds::FEData)
    fname = _resolve_name(ds, g.field)
    dh = ds.dh
    _check_full_domain(dh, "Gradient")
    dh_grad = _gradient_dofhandler(dh, fname, g.copy_fields)
    u_grad = Makie.lift(u -> _compute_gradient_values(dh, dh_grad, u, fname, g.copy_fields), ds.u)
    return _rebind(ds, dh_grad, u_grad)
end

####################
# Data derivations #
####################

"""
    vonmises(σ)

Von Mises equivalent stress `√(3/2 dev(σ) ⊡ dev(σ))` of a second-order tensor.
"""
function vonmises(σ::Union{Tensors.Tensor{2},Tensors.SymmetricTensor{2}})
    s = Tensors.dev(σ)
    return sqrt(3.0 / 2.0 * (s ⊡ s))
end

"""
    Component(i; input=:default, output=Symbol("x", i))

Filter extracting component `i` of a data array into a named scalar array.
"""
struct Component <: AbstractFilter
    i::Int
    input::Symbol
    output::Symbol
end
Component(i::Int; input::Symbol=:default, output::Symbol=Symbol("x", i)) = Component(i, input, output)
_valfun(f::Component) = x -> x[f.i]

"""
    Magnitude(; input=:default, output=:magnitude)

Filter computing the euclidean norm of a data array into a named scalar array.
"""
struct Magnitude <: AbstractFilter
    input::Symbol
    output::Symbol
end
Magnitude(; input::Symbol=:default, output::Symbol=:magnitude) = Magnitude(input, output)
_valfun(::Magnitude) = LinearAlgebra.norm

"""
    Norm1(; input=:default, output=:norm1)

Filter computing the 1-norm of a data array into a named scalar array.
"""
struct Norm1 <: AbstractFilter
    input::Symbol
    output::Symbol
end
Norm1(; input::Symbol=:default, output::Symbol=:norm1) = Norm1(input, output)
_valfun(::Norm1) = x -> sum(abs, x)

"""
    VonMises(; input=:default, output=:vonMises)

Filter computing the von Mises invariant ([`vonmises`](@ref)) of a
tensor-valued data array into a named scalar array. Note that this is only a
stress if the input array holds stresses — apply a constitutive law with
[`Derive`](@ref) first when starting from a displacement gradient.
"""
struct VonMises <: AbstractFilter
    input::Symbol
    output::Symbol
end
VonMises(; input::Symbol=:default, output::Symbol=:vonMises) = VonMises(input, output)
_valfun(::VonMises) = vonmises

"""
    Deviator(; input=:default, output=:deviator)

Filter computing the deviatoric part of a tensor-valued data array.
"""
struct Deviator <: AbstractFilter
    input::Symbol
    output::Symbol
end
Deviator(; input::Symbol=:default, output::Symbol=:deviator) = Deviator(input, output)
_valfun(::Deviator) = Tensors.dev

"""
    Derive(f; input=:default, output=:derived)

Generic derivation filter mapping each entry of a data array through `f`.
Point-data rows are passed to `f` as a scalar (1 component), `Vec` (spatial-dim
components) or `Tensor{2}` (spatial-dim² components); cell-data entries are
passed as-is. `f` may return a scalar, `Vec`, `Tensor` or `Tuple`.

# Example
```julia
σ(∇u) = 2G*dev(symmetric(∇u)) + K*tr(∇u)*one(∇u)
ds |> Gradient(:u) |> Derive(∇u -> vonmises(σ(∇u)); output=:σvM)
```
"""
struct Derive{F} <: AbstractFilter
    f::F
    input::Symbol
    output::Symbol
end
Derive(f; input::Symbol=:default, output::Symbol=:derived) = Derive(f, input, output)
_valfun(d::Derive) = d.f

for T in (:Component, :Magnitude, :Norm1, :VonMises, :Deviator, :Derive)
    @eval apply(f::$T, ds::FEData) = _apply_derivation(f, ds)
end

# Wrap a point-data row into the value it represents (see Derive docstring).
function _wrap_row(row, sdim::Int)
    n = length(row)
    n == 1 && return row[1]
    n == sdim && return Tensors.Vec{sdim}(NTuple{sdim,Float64}(row))
    n == sdim * sdim && return Tensors.Tensor{2,sdim}(NTuple{sdim * sdim,Float64}(row))
    return Tensors.Vec{n}(NTuple{n,Float64}(row))
end
_components(v::Number) = (v,)
_components(v) = Tuple(v)

function _rows_to_matrix(vf, A::AbstractMatrix, sdim::Int)
    n = size(A, 1)
    n == 0 && return Matrix{Float64}(undef, 0, 1)
    first_val = _components(vf(_wrap_row(view(A, 1, :), sdim)))
    out = Matrix{Float64}(undef, n, length(first_val))
    out[1, :] .= first_val
    for i in 2:n
        out[i, :] .= _components(vf(_wrap_row(view(A, i, :), sdim)))
    end
    return out
end

function _apply_derivation(f, ds::FEData{dim}) where {dim}
    input = _resolve_name(ds, f.input)
    assoc = _data_association(ds, input)
    assoc === :none && error("no data named :$input; available: $(_available_data(ds))")
    pd = copy(ds.point_data)
    cd = copy(ds.cell_data)
    vf = _valfun(f)
    if assoc === :point
        pd[f.output] = Makie.lift(A -> _rows_to_matrix(vf, A, dim), point_data(ds, input))
    else
        cd[f.output] = Makie.lift(v -> map(vf, v), cell_data(ds, input))
    end
    return _rebind(ds, ds.dh, ds.u; point_data=pd, cell_data=cd)
end

"""
    Threshold(; input=:default, output=:threshold, min=-Inf, max=Inf)

Filter copying a data array with values outside `[min, max]` replaced by `NaN`
(rendered as `nan_color`).
"""
struct Threshold <: AbstractFilter
    input::Symbol
    output::Symbol
    min::Float64
    max::Float64
end
Threshold(; input::Symbol=:default, output::Symbol=:threshold, min::Real=-Inf, max::Real=Inf) =
    Threshold(input, output, Float64(min), Float64(max))

function apply(t::Threshold, ds::FEData)
    input = _resolve_name(ds, t.input)
    assoc = _data_association(ds, input)
    assoc === :none && error("no data named :$input; available: $(_available_data(ds))")
    pd = copy(ds.point_data)
    cd = copy(ds.cell_data)
    clampnan(x) = t.min <= x <= t.max ? Float64(x) : NaN
    if assoc === :point
        pd[t.output] = Makie.lift(A -> map(clampnan, A), point_data(ds, input))
    else
        cd[t.output] = Makie.lift(v -> map(clampnan, v), cell_data(ds, input))
    end
    return _rebind(ds, ds.dh, ds.u; point_data=pd, cell_data=cd)
end
