# Layer 2b: the reactive filter graph.
#
# A filter maps FEData -> FEData. Filters are callable, so pipelines compose
# with |>:  ds |> WarpByVector(:u, 2.0) |> Gradient(:u) |> VonMises().
# Every apply returns a new FEData sharing the source solution observable, so
# pipelines stay live under FerriteViz.update! and fork without clobbering each
# other. Geometry-rebuilding filters (Refine, AddQuadraturePointData) rebuild
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

# Filter arguments may be given as a plain value or as an Observable the caller
# drives (a slider, a menu, ...); wrap the plain ones so the code downstream
# only ever deals with Observables.
make_observable(x) = Makie.Observable(x)
make_observable(x::Makie.Observable) = x

# Rebind dh/u on identical geometry (same grid ⇒ identical tessellation).
function _rebind(ds::FEData{dim}, dh, u::Makie.Observable;
                 point_data=Dict{Symbol,Makie.Observable}(), cell_data=Dict{Symbol,Makie.Observable}()) where {dim}
    return FEData{dim,typeof(dh),eltype(u[]),typeof(ds.topology),typeof(ds.source_u),typeof(ds.mesh),eltype(ds.all_triangles)}(
        dh, u, ds.source_u, ds.topology, ds.visible, ds.gridnodes, ds.coords, ds.coords_buffer,
        ds.all_triangles, ds.vis_triangles, ds.triangle_cell_map, ds.cell_triangle_offsets,
        ds.cell_vertex_offsets, ds.all_edges, ds.edge_cell_map, ds.cell_edge_offsets,
        ds.reference_coords, ds.mesh, point_data, cell_data)
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

`field` may name any vector-valued point-data array; the tessellation vertices
(surfaces and the `meshplot` wireframe) always follow. The original grid nodes
(`meshplot`'s node markers and labels) can only be displaced when the field is
a dof-backed field of the dof handler and stay put otherwise.
"""
struct WarpByVector{F<:Union{Symbol,Makie.Observable{Symbol}},S} <: AbstractFilter
    field::F
    scale::S
end
WarpByVector(field=:default) = WarpByVector(field, 1.0)

# Row `i` of a displacement container: a point-data Matrix or a Vector of
# per-node values (e.g. from `evaluate_at_grid_nodes`).
@inline _disp_component(d::AbstractMatrix, i::Int, j::Int) = d[i, j]
@inline _disp_component(d::AbstractVector, i::Int, j::Int) = d[i][j]

# Per-vertex displacement of WarpByVector, per frame. Function barrier: in the
# lift closures `dim` is captured as a plain `Int`, which would make
# `Point{dim,Float32}` a dynamic type application on every vertex — here it is
# a static parameter recovered from the points' element type.
function _displaced(points::Vector{GeometryBasics.Point{dim,Float32}}, d, scale::Real) where {dim}
    s = Float32(scale)
    return [points[i] .+ s .* GeometryBasics.Point{dim,Float32}(ntuple(j -> Float32(_disp_component(d, i, j)), Val(dim))) for i in eachindex(points)]
end

function apply(w::WarpByVector, ds::FEData{dim}) where {dim}
    scale = make_observable(w.scale)
    fname = w.field isa Makie.Observable ? w.field : make_observable(_resolve_name(ds, w.field))
    disp = _switching_point_data(ds, fname)
    size(disp[], 2) == dim || error("deformation field :$(fname[]) has $(size(disp[], 2)) components, expected $dim")
    coords = Makie.lift(ds.coords, disp, scale) do base, d, s
        size(d, 2) == dim || error("deformation field has $(size(d, 2)) components, expected $dim")
        _displaced(base, d, s)
    end
    gridnodes = Makie.lift(ds.gridnodes, ds.u, scale, fname) do nodes, u, s, fn
        fn = _resolve_name(ds, fn)
        # named (non-dof) arrays live on the tessellation vertices only
        fn in Ferrite.getfieldnames(ds.dh) || return nodes
        vals = Ferrite.evaluate_at_grid_nodes(ds.dh, u, fn)
        _displaced(nodes, vals, s)
    end
    coords_buffer = ShaderAbstractions.Buffer(coords)
    mesh = GeometryBasics.Mesh(coords_buffer, ds.vis_triangles)
    return FEData{dim,typeof(ds.dh),eltype(ds.u[]),typeof(ds.topology),typeof(ds.source_u),typeof(mesh),eltype(ds.all_triangles)}(
        ds.dh, ds.u, ds.source_u, ds.topology, ds.visible, gridnodes, coords, coords_buffer,
        ds.all_triangles, ds.vis_triangles, ds.triangle_cell_map, ds.cell_triangle_offsets,
        ds.cell_vertex_offsets, ds.all_edges, ds.edge_cell_map, ds.cell_edge_offsets,
        ds.reference_coords, mesh, copy(ds.point_data), copy(ds.cell_data))
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
        ds.cell_vertex_offsets, ds.all_edges, ds.edge_cell_map, ds.cell_edge_offsets,
        ds.reference_coords, mesh,
        Dict{Symbol,Makie.Observable}(), copy(ds.cell_data))
end

##########
# Refine #
##########

"""
    Refine()                       # automatic, what FEData applies by default
    Refine(n::Int; edges=n)
    Refine(; surface=nothing, edges=nothing)

Filter re-tessellating every cell from its reference shape with a subdivided
reference tessellation (see [`FerriteViz.subdivide`](@ref)): `surface` rounds
for the rendered triangles (each round quadruples them, refining the rendered
solution), `edges` rounds for the wireframe segments drawn by
[`meshplot`](@ref) (each round doubles them, refining the rendered geometry
edges). The subdivided reference vertices are mapped through the cell's
geometric interpolation, so curved (high-order) geometry and high-order
deformation render curved instead of as flat facets and straight chords.

The counts are absolute, not relative to the input's tessellation: `Refine(2)`
yields 2 subdivision rounds regardless of how the dataset was tessellated
before.

A count given as `nothing` is chosen per cell type: no subdivision when the
geometry and every dof field are (multi-)linear, otherwise 1 surface and
3 edge rounds. [`FEData`](@ref) applies this automatic mode by default —
construct with `adaptive=false` to opt out, e.g. to subdivide only one branch
of a pipeline:

```julia
ds = FEData(dh, u; adaptive=false)
meshplot(ds)                          # flat, cheap
solutionplot(ds |> Refine(2))         # this plot resolved finer
```

!!! note "Memory usage"
    Every surface round quadruples the rendered triangles and roughly triples
    the tessellation vertices (each of which carries solution values per
    field). The automatic mode therefore costs high-order cell types about 4×
    the memory of the flat tessellation; edge rounds are comparatively cheap
    (segments only double). On large high-order grids opt out with
    `FEData(dh, u; adaptive=false)`.

!!! note
    The choice made for a `nothing` count may change in a future release; such a
    change is breaking. Explicit counts are stable.

Rebuilds the geometry from the grid (a quadrature-point partition of
[`AddQuadraturePointData`](@ref) does not survive — nor would it gain anything
from refinement, its data being piecewise constant), so apply
[`WarpByVector`](@ref) *after* it; registered point data is dropped, cell data
survives.
"""
struct Refine <: AbstractFilter
    surface::Union{Nothing,Int}
    edges::Union{Nothing,Int}
end
Refine(n::Int; edges::Int=n) = Refine(n, edges)
Refine(; surface::Union{Nothing,Int}=nothing, edges::Union{Nothing,Int}=nothing) = Refine(surface, edges)

# Per-cell tessellation choice of a Refine filter, shared with the FEData
# constructor (which builds through this directly so the default application
# costs nothing over constructing flat and filtering afterwards).
function _tessellation_provider(f::Refine, dh::Ferrite.AbstractDofHandler)
    cache = Dict{Type,ReferenceTessellation}()
    return function (cell)
        return get!(cache, typeof(cell)) do
            degree = _render_degree(typeof(cell), dh)
            _cell_tessellation(reference_tessellation(Ferrite.getrefshape(cell)),
                               something(f.surface, _auto_surface_resolution(degree)),
                               something(f.edges, _auto_edge_resolution(degree)))
        end
    end
end

function apply(f::Refine, ds::FEData)
    out = _build_dataset(ds.dh, ds.u, ds.source_u, ds.topology, ds.visible,
                         _tessellation_provider(f, ds.dh))
    merge!(out.cell_data, ds.cell_data) # cell data is layout independent
    return out
end

##########################
# AddQuadraturePointData #
##########################

"""
    AddQuadraturePointData(qr, values; output=:qpdata, extract=identity)

Filter for internal variables, i.e. quantities that carry a value only at the
quadrature points and have no interpolation defining them anywhere else (in FEM
terms: L2 data). Every cell is partitioned into the Voronoi regions of its
quadrature points (see [`FerriteViz.qp_voronoi_tessellation`](@ref)) and each
region is filled with its quadrature point's value, giving a piecewise constant
("flat") rendering that neither averages over the cell nor smooths the data onto
a nodal field.

`qr` is a `Ferrite.QuadratureRule`, or a `Dict` mapping reference shapes to rules
for grids with mixed cell types.

`values` may be

  * a `Vector` of per-cell vectors (`values[cell][qp]`; `nqp` may differ per cell),
  * a `Matrix` (`values[cell, qp]`, requiring a uniform `nqp`), or
  * an `Observable` of either — updating it refreshes all open plots.

`extract` maps a stored entry to the plotted value, so Ferrite material states can
be handed over directly (`extract = s -> s.εₚ`). Scalars, `Vec`s and (symmetric)
second order tensors are supported. The result is registered as point data named
`output` and can be reduced further with [`VonMises`](@ref), [`Derive`](@ref), ...

Rebuilds the geometry, so apply [`WarpByVector`](@ref) *after* this filter.

# Example
```julia
FEData(dh, u) |> AddQuadraturePointData(qr, states; extract = s -> s.σ) |> VonMises(input = :qpdata)
```
"""
struct AddQuadraturePointData{Q,V,F} <: AbstractFilter
    qr::Q
    values::V
    output::Symbol
    extract::F
end

function AddQuadraturePointData(qr, values; output::Symbol=:qpdata, extract=identity)
    return AddQuadraturePointData(qr, make_observable(values), output, extract)
end

_qr_for(qr::Ferrite.QuadratureRule, ::Type) = qr
function _qr_for(qrs::AbstractDict, ::Type{RS}) where {RS}
    haskey(qrs, RS) ||
        error("no quadrature rule for reference shape $RS; add it to the `qr` mapping (have $(collect(keys(qrs))))")
    return qrs[RS]
end

# values[cell][qp] (ragged) or values[cell, qp] (uniform nqp)
_qp_ncells(values::AbstractMatrix) = size(values, 1)
_qp_ncells(values::AbstractVector) = length(values)
_qp_nqp(values::AbstractMatrix, ::Int) = size(values, 2)
_qp_nqp(values::AbstractVector, cell::Int) = length(values[cell])
_qp_at(values::AbstractMatrix, cell::Int, qp::Int) = values[cell, qp]
_qp_at(values::AbstractVector, cell::Int, qp::Int) = values[cell][qp]

# Voronoi assignment of an off-quadrature-point vertex (the appended wireframe
# edge vertices): the nearest quadrature point in reference space.
function _nearest_qp(ξ::Ferrite.Vec, points)
    best, bestdist = 1, Inf
    for (i, p) in enumerate(points)
        dist = sum(abs2, ξ - p)
        if dist < bestdist
            best, bestdist = i, dist
        end
    end
    return best
end

# Symmetric tensors are expanded to the full component order, so that the stored
# row has spatial-dim² entries and `_wrap_row` hands a `Tensor{2}` back to the
# derivation filters (VonMises, Deviator, ...).
_qp_components(v) = _components(v)
_qp_components(v::Tensors.SymmetricTensor{2,dim}) where {dim} = _components(convert(Tensors.Tensor{2,dim}, v))

function apply(f::AddQuadraturePointData, ds::FEData{dim}) where {dim}
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    ncells = length(cells)
    values = f.values[]
    _qp_ncells(values) == ncells ||
        error("quadrature point data must have one entry per cell ($ncells), got $(_qp_ncells(values))")

    tess_cache = Dict{Type,QPTessellation}()
    function tess_for(cell)
        RS = Ferrite.getrefshape(cell)
        return get!(() -> qp_voronoi_tessellation(RS, _qr_for(f.qr, RS)), tess_cache, RS)
    end
    # The rebuilt layout still carries the FE cell edges, so meshplot keeps
    # working downstream. Edge vertices are appended after the Voronoi vertices
    # and valued by their nearest quadrature point (the Voronoi assignment).
    edge_cache = Dict{Type,Any}()
    function edges_for(cell)
        return get!(edge_cache, typeof(cell)) do
            edge_res = _auto_edge_resolution(_render_degree(typeof(cell), ds.dh))
            edge_geometry(_cell_tessellation(reference_tessellation(Ferrite.getrefshape(cell)), 0, edge_res))
        end
    end

    cell_triangle_offsets = Vector{Int}(undef, ncells + 1)
    cell_vertex_offsets = Vector{Int}(undef, ncells + 1)
    cell_edge_offsets = Vector{Int}(undef, ncells + 1)
    cell_triangle_offsets[1] = 0
    cell_vertex_offsets[1] = 0
    cell_edge_offsets[1] = 0
    for (cell_id, cell) in enumerate(cells)
        tess = tess_for(cell)
        ecoords, eedges = edges_for(cell)
        nqp = length(Ferrite.getpoints(_qr_for(f.qr, Ferrite.getrefshape(cell))))
        _qp_nqp(values, cell_id) == nqp ||
            error("cell $cell_id carries $(_qp_nqp(values, cell_id)) quadrature values, but its rule has $nqp points")
        cell_triangle_offsets[cell_id+1] = cell_triangle_offsets[cell_id] + ntriangles(tess)
        cell_vertex_offsets[cell_id+1] = cell_vertex_offsets[cell_id] + nvertices(tess) + length(ecoords)
        cell_edge_offsets[cell_id+1] = cell_edge_offsets[cell_id] + length(eedges)
    end
    num_triangles = cell_triangle_offsets[end]
    num_verts = cell_vertex_offsets[end]
    num_edges = cell_edge_offsets[end]

    triangles = Matrix{Int}(undef, num_triangles, 3)
    triangle_cell_map = Vector{Int}(undef, num_triangles)
    all_edges = Vector{NTuple{2,Int}}(undef, num_edges)
    edge_cell_map = Vector{Int}(undef, num_edges)
    physical_coords = Vector{GeometryBasics.Point{dim,Float32}}(undef, num_verts)
    reference_coords = zeros(Float64, num_verts, dim)
    # static vertex -> (cell, quadrature point) map; the value lift is a gather
    vertex_cell = Vector{Int}(undef, num_verts)
    vertex_qp = Vector{Int}(undef, num_verts)

    for (cell_id, cell) in enumerate(cells)
        tess = tess_for(cell)
        ecoords, eedges = edges_for(cell)
        qpoints = Ferrite.getpoints(_qr_for(f.qr, Ferrite.getrefshape(cell)))
        ip_geo = Ferrite.geometric_interpolation(typeof(cell))
        node_coords = Ferrite.getcoordinates(grid, cell_id)
        coff = cell_vertex_offsets[cell_id]
        for (k, ξ) in enumerate(tess.coords)
            x = geometric_map(ip_geo, node_coords, ξ)
            physical_coords[coff+k] = GeometryBasics.Point{dim,Float32}(x...)
            for d in 1:length(ξ)
                reference_coords[coff+k, d] = ξ[d]
            end
            vertex_cell[coff+k] = cell_id
            vertex_qp[coff+k] = tess.vertex_qp[k]
        end
        evoff = coff + nvertices(tess)
        for (k, ξ) in enumerate(ecoords)
            x = geometric_map(ip_geo, node_coords, ξ)
            physical_coords[evoff+k] = GeometryBasics.Point{dim,Float32}(x...)
            for d in 1:length(ξ)
                reference_coords[evoff+k, d] = ξ[d]
            end
            vertex_cell[evoff+k] = cell_id
            vertex_qp[evoff+k] = _nearest_qp(ξ, qpoints)
        end
        toff = cell_triangle_offsets[cell_id]
        for (t, tri) in enumerate(tess.triangles)
            for j in 1:3
                triangles[toff+t, j] = tri[j] + coff
            end
            triangle_cell_map[toff+t] = cell_id
        end
        eoff = cell_edge_offsets[cell_id]
        for (e, edge) in enumerate(eedges)
            all_edges[eoff+e] = (edge[1] + evoff, edge[2] + evoff)
            edge_cell_map[eoff+e] = cell_id
        end
    end

    all_triangles = convert(Vector{GeometryBasics.GLTriangleFace}, Makie.to_triangles(triangles))
    vis_triangles = ShaderAbstractions.Buffer(Makie.Observable(_visibility_triangles(all_triangles, ds.visible, triangle_cell_map)))
    coords = Makie.Observable(physical_coords)
    coords_buffer = ShaderAbstractions.Buffer(coords)
    mesh = GeometryBasics.Mesh(coords_buffer, vis_triangles)
    # Rebuilding the geometry normally invalidates the upstream point data. A
    # second AddQuadraturePointData with the same rule, however, lays out exactly
    # the same vertices, so those arrays stay valid — which is what lets several
    # quadrature point quantities be combined (e.g. σ and εᵖ in one Derive).
    same_layout = size(ds.reference_coords) == size(reference_coords) &&
                  ds.cell_vertex_offsets == cell_vertex_offsets &&
                  ds.reference_coords == reference_coords
    out = FEData{dim,typeof(ds.dh),eltype(ds.u[]),typeof(ds.topology),typeof(ds.source_u),typeof(mesh),eltype(all_triangles)}(
        ds.dh, ds.u, ds.source_u, ds.topology, ds.visible, ds.gridnodes, coords, coords_buffer,
        all_triangles, vis_triangles, triangle_cell_map, cell_triangle_offsets,
        cell_vertex_offsets, all_edges, edge_cell_map, cell_edge_offsets,
        reference_coords, mesh,
        same_layout ? copy(ds.point_data) : Dict{Symbol,Makie.Observable}(), copy(ds.cell_data))

    ncomponents = length(_qp_components(f.extract(_qp_at(values, 1, 1))))
    data = Makie.lift(f.values) do vals
        _qp_ncells(vals) == ncells ||
            error("quadrature point data must have one entry per cell ($ncells), got $(_qp_ncells(vals))")
        A = Matrix{Float64}(undef, num_verts, ncomponents)
        for v in 1:num_verts
            components = _qp_components(f.extract(_qp_at(vals, vertex_cell[v], vertex_qp[v])))
            for d in 1:ncomponents
                A[v, d] = components[d]
            end
        end
        return A
    end
    set_point_data!(out, f.output, data)
    return out
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
    ExtractComponent(i; input=:default, output=Symbol("x", i))

Filter extracting component `i` of a data array into a named scalar array.
"""
struct ExtractComponent <: AbstractFilter
    i::Int
    input::Symbol
    output::Symbol
end
ExtractComponent(i::Int; input::Symbol=:default, output::Symbol=Symbol("x", i)) = ExtractComponent(i, input, output)
_valfun(f::ExtractComponent) = x -> x[f.i]

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

Generic derivation filter mapping the entries of one or more data arrays through `f`.

`input` is a single name or a vector of names. `f` receives one argument per input,
taken from the same tessellation vertex (point data) or the same cell (cell data),
so `Derive(g; input=[:a, :b])` calls `g(a_i, b_i)`. All inputs must be of the same
kind, either all point data or all cell data.

Point-data rows are passed to `f` as a scalar (1 component), `Vec` (spatial-dim
components) or `Tensor{2}` (spatial-dim² components). Rows holding a tensor of a
different dimension than the grid — a shell or plane-strain problem carrying 3D
stresses on a 2D grid — are recognised by their component count as
`Tensor{2,2}` (4), `SymmetricTensor{2,3}` (6) or `Tensor{2,3}` (9). Cell-data
entries are passed as-is. `f` may return a scalar, `Vec`, `Tensor` or `Tuple`.

# Examples
```julia
σ(∇u) = 2G*dev(symmetric(∇u)) + K*tr(∇u)*one(∇u)
ds |> Gradient(:u) |> Derive(∇u -> vonmises(σ(∇u)); output=:σvM)

# several inputs -> one argument each
ds |> Derive((σ, εᵖ) -> σ ⊡ εᵖ; input=[:σ, :εᵖ], output=:dissipation)
```
"""
struct Derive{F} <: AbstractFilter
    f::F
    input::Vector{Symbol}
    output::Symbol
end
Derive(f; input=:default, output::Symbol=:derived) = Derive(f, _derive_inputs(input), output)
_derive_inputs(name::Symbol) = [name]
_derive_inputs(names) = collect(Symbol, names)
_valfun(d::Derive) = d.f

for T in (:ExtractComponent, :Magnitude, :Norm1, :VonMises, :Deviator, :Derive)
    @eval apply(f::$T, ds::FEData) = _apply_derivation(f, ds)
end

# Wrap a point-data row into the value it represents (see Derive docstring).
#
# The component count is matched against the spatial dimension first, so an
# ordinary vector or tensor field on the grid keeps its natural meaning. A row
# may however carry a tensor of a *different* dimension than the grid: shells
# and plane-strain problems store full 3D stresses on a grid whose spatial
# dimension is 2, and those counts (4, 6, 9) are unambiguous once the
# dimension-matched cases above have been ruled out.
#
# Anything left over cannot be interpreted — Tensors.jl has no `Vec` beyond
# dimension 3, so the old fallback `Vec{n}` raised a MethodError from deep
# inside Tensors for every such row. Failing here with the component count
# named is considerably more useful.
function _wrap_row(row, sdim::Int)
    n = length(row)
    n == 1 && return row[1]
    n == sdim && return Tensors.Vec{sdim}(NTuple{sdim,Float64}(row))
    n == sdim * sdim && return Tensors.Tensor{2,sdim}(NTuple{sdim * sdim,Float64}(row))
    n == 9 && return Tensors.Tensor{2,3}(NTuple{9,Float64}(row))
    n == 6 && return Tensors.SymmetricTensor{2,3}(NTuple{6,Float64}(row))
    n == 4 && return Tensors.Tensor{2,2}(NTuple{4,Float64}(row))
    n <= 3 && return Tensors.Vec{n}(NTuple{n,Float64}(row))
    error("cannot interpret a $n-component data row on a $(sdim)D grid as a scalar, vector or " *
          "second order tensor; reduce it first, e.g. with ExtractComponent(i) or Derive")
end
_components(v::Number) = (v,)
_components(v) = Tuple(v)

_rows_to_matrix(vf, A::AbstractMatrix, sdim::Int) = _rows_to_matrix(vf, (A,), sdim)

# `As` holds one point-data array per input; row `i` of each is wrapped and
# handed to `vf` as a separate argument.
function _rows_to_matrix(vf, As::Tuple, sdim::Int)
    n = size(first(As), 1)
    all(A -> size(A, 1) == n, As) ||
        error("derivation inputs must have the same number of rows, got $(map(A -> size(A, 1), As))")
    n == 0 && return Matrix{Float64}(undef, 0, 1)
    wrap(i) = map(A -> _wrap_row(view(A, i, :), sdim), As)
    first_val = _components(vf(wrap(1)...))
    out = Matrix{Float64}(undef, n, length(first_val))
    out[1, :] .= first_val
    for i in 2:n
        out[i, :] .= _components(vf(wrap(i)...))
    end
    return out
end

# The unary reductions carry a single input name, Derive a vector of them.
_inputs(f) = (f.input,)
_inputs(d::Derive) = Tuple(d.input)

function _apply_derivation(f, ds::FEData{dim}) where {dim}
    names = map(n -> _resolve_name(ds, n), _inputs(f))
    assocs = map(n -> _data_association(ds, n), names)
    for (name, assoc) in zip(names, assocs)
        assoc === :none && error("no data named :$name; available: $(_available_data(ds))")
    end
    all(a -> a === first(assocs), assocs) ||
        error("cannot derive from a mix of point and cell data: " *
              join(("$n => $a" for (n, a) in zip(names, assocs)), ", "))
    pd = copy(ds.point_data)
    cd = copy(ds.cell_data)
    vf = _valfun(f)
    if first(assocs) === :point
        pd[f.output] = Makie.lift((As...) -> _rows_to_matrix(vf, As, dim),
                                  map(n -> point_data(ds, n), names)...)
    else
        cd[f.output] = Makie.lift((vs...) -> map(vf, vs...),
                                  map(n -> cell_data(ds, n), names)...)
    end
    return _rebind(ds, ds.dh, ds.u; point_data=pd, cell_data=cd)
end

"""
    Threshold(; input=:default, output=:threshold, min=-Inf, max=Inf)

Filter copying a data array with values outside `[min, max]` replaced by `NaN`.

!!! note
    This masks values, it does not remove geometry. The tessellation is passed
    through unchanged and the `NaN`s reach Makie as *colors*, so the affected
    triangles are still drawn — in `nan_color` (`:red` by default), and blended
    across a triangle whose other vertices are inside the range. Set
    `nan_color=:transparent` on the representation to hide them. Removing cells
    from the mesh is what [`CrinkleClip`](@ref) does.
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
