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

# Filters that rebuild whole cells from the grid cannot run on a dataset whose
# cells were partially cut away — they would resurrect the removed geometry.
function _check_whole_cells(ds::FEData, what::String)
    ds.cells_intact ||
        error("$what rebuilds whole cells from the grid, but this dataset's cells were cut " *
              "(by Clip or ExtractIsosurfaces); apply $what before cutting")
    return nothing
end

# Filter arguments may be given as a plain value or as an Observable the caller
# drives (a slider, a menu, ...); wrap the plain ones so the code downstream
# only ever deals with Observables.
make_observable(x) = Makie.Observable(x)
make_observable(x::Makie.Observable) = x

# Rebind dh/u on identical geometry (same grid ⇒ identical tessellation).
function _rebind(ds::FEData, dh, u::Makie.Observable;
                 point_data=Dict{Symbol,Makie.Observable}(), cell_data=Dict{Symbol,Makie.Observable}(),
                 derivations=copy(ds.point_derivations))
    return _derive(ds; dh, u, point_data, cell_data, point_derivations=derivations)
end

##############
# Adaptivity #
##############

# Applying an `Adaptivity` (defined next to FEData, see its docstring)
# re-configures how the dataset's plots refine: the result is the same
# dataset — geometry, data and provenance shared — with the given settings
# and a fresh substrate cache. This is also what the constructor's
# `adaptivity` keyword stores; `FEData(dh, u; adaptivity=a)` and
# `FEData(dh, u; adaptivity=false) |> a` are equivalent.
apply(a::Adaptivity, ds::FEData) = _derive(ds; adaptivity=a)

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
    return _derive(ds; coords, gridnodes,
                   deformation=vcat(ds.deformation, [Deformation(ds.dh, ds.u, fname, scale)]))
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
    solid = copy(ds.solid)
    visible = copy(ds.visible)
    # a neighbor counts as kept only if it survived every previous filter too,
    # so boundaries created by an earlier (crinkle or exact) clip are revealed
    kept(cell_id) = ds.solid[cell_id] && c.decision(grid, cell_id)
    for cell_id in 1:Ferrite.getncells(grid)
        if kept(cell_id)
            cell_neighbors = Ferrite.getneighborhood(ds.topology, grid, Ferrite.CellIndex(cell_id))
            visible[cell_id] = !all(kept, cell_neighbors) || ds.visible[cell_id]
        else
            solid[cell_id] = false
            visible[cell_id] = false
        end
    end
    # registered point data survives (the vertex layout is unchanged); cached
    # dof-field transfers are dropped and re-resolve against the new masks
    return _derive(ds; solid, visible, point_data=_registered_point_data(ds))
end

########
# Clip #
########

"""
    Clip(plane::ClipPlane)

Filter cutting a 3D dataset exactly at `plane`, keeping the side
`normal ⋅ x ≤ distance`. Unlike [`CrinkleClip`](@ref) (which hides whole
cells), the finite elements themselves are cut: surface triangles and
wireframe edges are clipped at the plane and the cross-section is capped with
triangles showing the interior field values. Cap triangles belong to the cell
they cut through, cut cells become visible (previously hidden interior cells
included), and the remaining per-cell volume is carried along — so a second
`Clip` cuts the already-clipped volume, quadrature-point Voronoi regions are
cut exactly along their walls, and [`ExtractIsosurfaces`](@ref) of a clipped
dataset stays inside the kept volume.

The cut *topology* (which edges cross the plane, and the interpolation
weights) is fixed when the filter is applied: positions and data stay reactive
under [`FerriteViz.update!`](@ref) — cut vertices follow their parent edges —
but a deformation that moves vertices across the plane needs the filter to be
re-applied to re-cut. Dof-backed fields evaluate exactly at the cut positions;
registered point data is interpolated linearly along the cut edges.

Nonlinear geometry is treated as linear (the plane cuts the tessellation's
straight edges); apply [`Refine`](@ref) *before* clipping to resolve
curvature. Plots of a cut dataset always draw the static tessellation (the
error-adaptive path refines whole cells and cannot represent cut ones).
"""
struct Clip{T} <: AbstractFilter
    plane::ClipPlane{T}
end

apply(::Clip, ::FEData{dim}) where {dim} =
    error("Clip supports only 3D datasets (got spatial dimension $dim); in 2D consider Threshold or ExtractIsosurfaces")

function apply(c::Clip, ds::FEData{3})
    nn = LinearAlgebra.norm(c.plane.normal)
    (isfinite(nn) && nn > 0 && isfinite(c.plane.distance)) ||
        error("Clip plane must have a finite nonzero normal and a finite distance")
    n = c.plane.normal / nn
    d = c.plane.distance / nn

    pc = ds.coords[]
    g = Vector{Float64}(undef, length(pc))
    lo = Tensors.Vec(Inf, Inf, Inf)
    hi = -lo
    maxabs = 0.0
    for (i, p) in enumerate(pc)
        x = Tensors.Vec{3,Float64}(NTuple{3,Float64}(p))
        g[i] = x ⋅ n - d
        if isfinite(g[i])
            lo = min.(lo, x)
            hi = max.(hi, x)
            maxabs = max(maxabs, maximum(abs, x))
        end
    end
    # one dataset-global classification tolerance (cell-local tolerances could
    # classify the duplicated copies of a shared vertex differently and crack
    # the surface); the eps term floors it at Float32 roundoff of the data
    tol = 1e-6 * LinearAlgebra.norm(hi - lo) + 4 * Float64(eps(Float32(maxabs)))
    snap!(g, tol)

    grid = Ferrite.get_grid(ds.dh)
    ncells = Ferrite.getncells(grid)
    pool = CutVertexPool(pc)
    out_tris = NTuple{3,Int}[]
    out_tri_cells = Int[]
    out_tets = NTuple{4,Int}[]
    out_tet_cells = Int[]
    out_edges = NTuple{2,Int}[]
    out_edge_cells = Int[]
    solid = copy(ds.solid)
    visible = copy(ds.visible)
    cell_triangle_offsets = zeros(Int, ncells + 1)
    cell_vertex_offsets = zeros(Int, ncells + 1)
    cell_simplex_offsets = zeros(Int, ncells + 1)
    cell_edge_offsets = zeros(Int, ncells + 1)
    any_cut = false

    for cell in 1:ncells
        if ds.solid[cell]
            verts = vertices_on_cell(ds, cell)
            has_volume = ds.cell_simplex_offsets[cell+1] > ds.cell_simplex_offsets[cell]
            # degeneracy pruning is cell-local (a global threshold could
            # discard valid geometry in the small cells of a graded mesh);
            # only the classification tolerance above is global
            tol_len = 1e-9 * _cell_diag(pc, verts)
            area_tol = tol_len^2
            vol_tol = tol_len^3
            nneg = count(v -> g[v] < 0, verts)
            nzero = count(v -> g[v] == 0, verts)
            nin = nneg + nzero
            # a volume cell whose kept part has no interior (nothing strictly
            # inside) is removed — keeping it would only duplicate on-plane
            # faces the inside neighbor already draws; measure-zero shells on
            # the plane are kept
            kept = has_volume ? nneg > 0 : nin > 0
            if !kept
                solid[cell] = false
                visible[cell] = false
            elseif nin == length(verts)
                # fully on the kept side: carry the cell over unchanged
                for v in verts
                    out_vertex!(pool, v)
                end
                for t in triangles_on_cell(ds, cell)
                    tri = ds.all_triangles[t]
                    push!(out_tris, (pool.remap[convert(Int, tri[1])], pool.remap[convert(Int, tri[2])], pool.remap[convert(Int, tri[3])]))
                    push!(out_tri_cells, cell)
                end
                for s in simplices_on_cell(ds, cell)
                    tet = ds.simplices[s]
                    push!(out_tets, (pool.remap[tet[1]], pool.remap[tet[2]], pool.remap[tet[3]], pool.remap[tet[4]]))
                    push!(out_tet_cells, cell)
                end
                for e in edges_on_cell(ds, cell)
                    edge = ds.all_edges[e]
                    push!(out_edges, (pool.remap[edge[1]], pool.remap[edge[2]]))
                    push!(out_edge_cells, cell)
                end
                nzero > 0 && (visible[cell] = true)   # touches the plane: reveal
            else
                any_cut = true
                visible[cell] = true
                for t in triangles_on_cell(ds, cell)
                    tri = ds.all_triangles[t]
                    clip_triangle!(pool, out_tris, out_tri_cells,
                                   (convert(Int, tri[1]), convert(Int, tri[2]), convert(Int, tri[3])),
                                   cell, g, area_tol)
                end
                # caps go into the same triangle list, keeping the cell's
                # triangles contiguous
                for s in simplices_on_cell(ds, cell)
                    clip_tet!(pool, out_tets, out_tet_cells, out_tris, out_tri_cells,
                              ds.simplices[s], cell, g, vol_tol, area_tol)
                end
                for e in edges_on_cell(ds, cell)
                    clip_edge!(pool, out_edges, out_edge_cells, ds.all_edges[e], cell, g)
                end
            end
        end
        cell_vertex_offsets[cell+1] = nvertices(pool)
        cell_triangle_offsets[cell+1] = length(out_tris)
        cell_simplex_offsets[cell+1] = length(out_tets)
        cell_edge_offsets[cell+1] = length(out_edges)
    end

    combos = pool.combos
    tri_matrix = Matrix{Int}(undef, length(out_tris), 3)
    for (t, tri) in enumerate(out_tris), j in 1:3
        tri_matrix[t, j] = tri[j]
    end
    all_triangles = convert(Vector{GeometryBasics.GLTriangleFace}, Makie.to_triangles(tri_matrix))
    coords = Makie.lift(p -> combine_points(combos, p), ds.coords)
    reference_coords = combine_rows(combos, ds.reference_coords)
    point_data = Dict{Symbol,Makie.Observable}(
        k => Makie.lift(A -> combine_rows(combos, A), v) for (k, v) in _registered_point_data(ds))
    return _derive(ds; solid, visible, cells_intact=ds.cells_intact && !any_cut,
                   coords, all_triangles, triangle_cell_map=out_tri_cells,
                   cell_triangle_offsets, cell_vertex_offsets,
                   simplices=out_tets, simplex_cell_map=out_tet_cells, cell_simplex_offsets,
                   all_edges=out_edges, edge_cell_map=out_edge_cells, cell_edge_offsets,
                   reference_coords, point_data)
end

function _bbox_diag(coords)
    isempty(coords) && return 0.0
    lo = hi = Float64.(coords[1])
    for p in coords
        x = Float64.(p)
        any(!isfinite, x) && continue
        lo = min.(lo, x)
        hi = max.(hi, x)
    end
    return LinearAlgebra.norm(hi .- lo)
end

_cell_diag(coords, verts) = isempty(verts) ? 0.0 : _bbox_diag(view(coords, verts))

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
3 edge rounds. [`FEData`](@ref) itself always builds the flat base
tessellation (curved rendering comes from the adaptive path); uniform static
subdivision is an explicit composition, e.g. for one branch of a pipeline:

```julia
ds = FEData(dh, u; adaptivity=false)
meshplot(ds)                          # flat, cheap
solutionplot(ds |> Refine(2))         # this plot resolved finer
```

!!! note "Memory usage"
    Every surface round quadruples the rendered triangles and roughly triples
    the tessellation vertices (each of which carries solution values per
    field). The automatic mode therefore costs high-order cell types about 4×
    the memory of the flat tessellation; edge rounds are comparatively cheap
    (segments only double).

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

# How a cell type's tessellation is chosen, top to bottom:
#
#   _tessellation_provider          one cached tessellation per cell type
#   └─ _build_cell_tessellation     base shape -> subdivided tessellation
#      ├─ _pick_subdivision_rounds  explicit Refine counts win, `nothing`
#      │  │                         falls back to the automatic defaults
#      │  ├─ _max_render_order      highest order of geometry and dof fields
#      │  └─ _default_*_rounds      flat for order 1, subdivided above
#      ├─ reference_tessellation    flat tessellation of the reference shape
#      └─ _subdivided_tessellation  applies the subdivision rounds

# Highest polynomial order among the dof fields. A quadratic displacement on a
# linear grid bends edges just like curved geometry does, so the fields count
# toward the render order alongside the geometric interpolation.
function _max_field_order(dh::Ferrite.DofHandler)
    order = 1
    for sdh in dh.subdofhandlers, name in sdh.field_names
        order = max(order, Ferrite.getorder(Ferrite.getfieldinterpolation(sdh, name)))
    end
    return order
end
_max_field_order(::Ferrite.AbstractDofHandler) = 1

# Highest polynomial order a cell of this type may have to render: its
# geometric interpolation's order or any dof field's, whichever is larger.
function _max_render_order(celltype::Type{<:Ferrite.AbstractCell}, dh::Ferrite.AbstractDofHandler)
    return max(Ferrite.getorder(Ferrite.geometric_interpolation(celltype)), _max_field_order(dh))
end

# Defaults of Refine's automatic mode: cells that render (multi-)linearly are
# exact on the flat base tessellation and get no subdivision; higher-order
# cells get 1 surface round (4× the triangles) and 3 edge rounds (each cell
# edge drawn as 8 segments).
_default_surface_rounds(render_order::Int) = render_order > 1 ? 1 : 0
_default_edge_rounds(render_order::Int) = render_order > 1 ? 3 : 0

# The subdivision rounds a Refine filter applies to one cell type: counts set
# explicitly on the filter are used as given, counts left as `nothing` fall
# back to the automatic defaults for the cell type's render order.
function _pick_subdivision_rounds(f::Refine, celltype::Type{<:Ferrite.AbstractCell}, dh::Ferrite.AbstractDofHandler)
    render_order = _max_render_order(celltype, dh)
    surface_rounds = f.surface === nothing ? _default_surface_rounds(render_order) : f.surface
    edge_rounds = f.edges === nothing ? _default_edge_rounds(render_order) : f.edges
    return surface_rounds, edge_rounds
end

function _build_cell_tessellation(f::Refine, celltype::Type{<:Ferrite.AbstractCell}, dh::Ferrite.AbstractDofHandler)
    surface_rounds, edge_rounds = _pick_subdivision_rounds(f, celltype, dh)
    base = reference_tessellation(Ferrite.getrefshape(celltype))
    return _subdivided_tessellation(base, surface_rounds, edge_rounds)
end

# Per-cell tessellation lookup handed to `_build_dataset`, shared with the
# FEData constructor (which builds through this directly so the default
# application costs nothing over constructing flat and filtering afterwards).
# All cells of one type share the same reference tessellation, so the build
# runs once per cell type and is cached.
function _tessellation_provider(f::Refine, dh::Ferrite.AbstractDofHandler)
    cache = Dict{Type,ReferenceTessellation}()
    tessellation_for_cell(cell) = get!(() -> _build_cell_tessellation(f, typeof(cell), dh), cache, typeof(cell))
    return tessellation_for_cell
end

function apply(f::Refine, ds::FEData)
    _check_whole_cells(ds, "Refine")
    out = _build_dataset(ds.dh, ds.u, ds.source_u, ds.topology, ds.visible,
                         _tessellation_provider(f, ds.dh); adaptivity=ds.adaptivity, solid=ds.solid)
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
    _check_whole_cells(ds, "AddQuadraturePointData")
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
            edge_rounds = _default_edge_rounds(_max_render_order(typeof(cell), ds.dh))
            base = reference_tessellation(Ferrite.getrefshape(cell))
            edge_geometry(_subdivided_tessellation(base, 0, edge_rounds))
        end
    end

    cell_triangle_offsets = Vector{Int}(undef, ncells + 1)
    cell_vertex_offsets = Vector{Int}(undef, ncells + 1)
    cell_edge_offsets = Vector{Int}(undef, ncells + 1)
    cell_simplex_offsets = Vector{Int}(undef, ncells + 1)
    cell_triangle_offsets[1] = 0
    cell_vertex_offsets[1] = 0
    cell_edge_offsets[1] = 0
    cell_simplex_offsets[1] = 0
    for (cell_id, cell) in enumerate(cells)
        # removed (e.g. crinkle-clipped) cells contribute no geometry
        if !ds.solid[cell_id]
            cell_triangle_offsets[cell_id+1] = cell_triangle_offsets[cell_id]
            cell_vertex_offsets[cell_id+1] = cell_vertex_offsets[cell_id]
            cell_edge_offsets[cell_id+1] = cell_edge_offsets[cell_id]
            cell_simplex_offsets[cell_id+1] = cell_simplex_offsets[cell_id]
            continue
        end
        tess = tess_for(cell)
        ecoords, eedges = edges_for(cell)
        nqp = length(Ferrite.getpoints(_qr_for(f.qr, Ferrite.getrefshape(cell))))
        _qp_nqp(values, cell_id) == nqp ||
            error("cell $cell_id carries $(_qp_nqp(values, cell_id)) quadrature values, but its rule has $nqp points")
        cell_triangle_offsets[cell_id+1] = cell_triangle_offsets[cell_id] + ntriangles(tess)
        cell_vertex_offsets[cell_id+1] = cell_vertex_offsets[cell_id] + nvertices(tess) + length(ecoords)
        cell_edge_offsets[cell_id+1] = cell_edge_offsets[cell_id] + length(eedges)
        cell_simplex_offsets[cell_id+1] = cell_simplex_offsets[cell_id] + nsimplices(tess)
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
    S = NTuple{dim + 1,Int}
    simplices = Vector{S}(undef, cell_simplex_offsets[end])
    simplex_cell_map = Vector{Int}(undef, cell_simplex_offsets[end])
    # static vertex -> (cell, quadrature point) map; the value lift is a gather
    vertex_cell = Vector{Int}(undef, num_verts)
    vertex_qp = Vector{Int}(undef, num_verts)

    for (cell_id, cell) in enumerate(cells)
        ds.solid[cell_id] || continue
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
        soff = cell_simplex_offsets[cell_id]
        for (s, simplex) in enumerate(tess.simplices)
            simplices[soff+s] = simplex .+ coff
            simplex_cell_map[soff+s] = cell_id
        end
    end

    all_triangles = convert(Vector{GeometryBasics.GLTriangleFace}, Makie.to_triangles(triangles))
    coords = Makie.Observable(physical_coords)
    # Rebuilding the geometry normally invalidates the upstream point data. A
    # second AddQuadraturePointData with the same rule, however, lays out exactly
    # the same vertices, so those arrays stay valid — which is what lets several
    # quadrature point quantities be combined (e.g. σ and εᵖ in one Derive).
    same_layout = size(ds.reference_coords) == size(reference_coords) &&
                  ds.cell_vertex_offsets == cell_vertex_offsets &&
                  ds.reference_coords == reference_coords
    out = _derive(ds; coords, all_triangles, triangle_cell_map,
                  cell_triangle_offsets, cell_vertex_offsets,
                  simplices, simplex_cell_map, cell_simplex_offsets,
                  all_edges, edge_cell_map, cell_edge_offsets,
                  reference_coords,
                  point_data=same_layout ? copy(ds.point_data) : Dict{Symbol,Makie.Observable}(),
                  point_derivations=same_layout ? copy(ds.point_derivations) : Dict{Symbol,DerivedPointData}(),
                  qp_partition=QPPartition(f.qr, f.values, f.extract, f.output),
                  # the geometry was rebuilt from the grid: upstream warps are
                  # baked into nothing here — apply WarpByVector after this
                  # filter (the static coords and the adaptive substrate then
                  # agree on the deformation)
                  deformation=Deformation[])

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

# The symbolic side of a pointwise derivation: how to recompute `vf` at an
# arbitrary (cell, ξ) — the array registered by `_apply_derivation` only knows
# the static tessellation's vertices. Inputs resolve to dof fields (captured
# with their handler and solution, like a warp's `Deformation`) or to earlier
# records; an input with no pointwise meaning — a raw registered array,
# quadrature-point data — returns `nothing`, and the derived quantity stays
# bound to the static tessellation (the adaptive path keeps refusing it).
function _derivation_record(ds::FEData, vf, names)
    inputs = Any[]
    for name in names
        if name in Ferrite.getfieldnames(ds.dh)
            push!(inputs, FieldSource(ds.dh, ds.u, name))
        elseif haskey(ds.point_derivations, name)
            push!(inputs, ds.point_derivations[name])
        else
            return nothing
        end
    end
    return DerivedPointData(vf, inputs)
end

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
    dv = copy(ds.point_derivations)
    vf = _valfun(f)
    if first(assocs) === :point
        pd[f.output] = Makie.lift((As...) -> _rows_to_matrix(vf, As, dim),
                                  map(n -> point_data(ds, n), names)...)
        rec = _derivation_record(ds, vf, names)
        rec === nothing ? delete!(dv, f.output) : (dv[f.output] = rec)
    else
        cd[f.output] = Makie.lift((vs...) -> map(vf, vs...),
                                  map(n -> cell_data(ds, n), names)...)
    end
    return _rebind(ds, ds.dh, ds.u; point_data=pd, cell_data=cd, derivations=dv)
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
    dv = copy(ds.point_derivations)
    clampnan(x) = t.min <= x <= t.max ? Float64(x) : NaN
    if assoc === :point
        pd[t.output] = Makie.lift(A -> map(clampnan, A), point_data(ds, input))
        rec = _derivation_record(ds, clampnan, (input,))
        rec === nothing ? delete!(dv, t.output) : (dv[t.output] = rec)
    else
        cd[t.output] = Makie.lift(v -> map(clampnan, v), cell_data(ds, input))
    end
    return _rebind(ds, ds.dh, ds.u; point_data=pd, cell_data=cd, derivations=dv)
end
