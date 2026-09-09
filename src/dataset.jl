# Layer 2a: the FEData source.
#
# FEData wraps a dof handler and a solution Observable together with the static
# triangulation built from the tessellation interface.
#
# Two Makie mechanisms are used directly here rather than left to the recipes,
# and both are load bearing:
#
#  * Observables. Every quantity that can change after a plot exists (the dof
#    vector, the coordinates, each named data array) is an Observable, and
#    filters build new datasets by `lift`ing from their input's observables
#    instead of copying values. That is what makes a pipeline reactive: one
#    FerriteViz.update! on the root solution propagates through every filter
#    down to the open plots, and two pipelines forked off the same source stay
#    independent because each `lift` has its own output. Storing plain arrays
#    would mean rebuilding the pipeline (and the plots) on every time step.
#
#  * ShaderAbstractions.Buffers. The vertex coordinates and the triangle index
#    list are wrapped in Buffers, which are then shared into the
#    GeometryBasics.Mesh handed to Makie. A Buffer is the CPU-side handle of a
#    GPU buffer, so writing into it uploads in place: an updated solution moves
#    the existing vertices instead of allocating a new mesh and forcing Makie to
#    tear down and re-upload the plot. Meshes derived by filters that do not
#    change the geometry (Gradient, the derivation filters) deliberately reuse
#    the *same* buffer objects, so several plots of one pipeline share a single
#    GPU upload. Filters that do rebuild the geometry (Refine,
#    AddQuadraturePointData) allocate fresh buffers, which is why upstream point
#    data cannot survive them unless the vertex layout is reproduced exactly.
#
# The consequence for anyone touching this file: never replace an Observable's
# or a Buffer's *content* by assigning a new object to the field — set the
# observable (`obs[] = ...`) so the listeners downstream fire.

# One WarpByVector application, recorded for consumers that need the
# deformation as a continuous function rather than baked into coordinates.
# `dh` and `u` are those of the dataset the warp was *applied to*: a later
# filter may rebind both (Gradient replaces the dof handler entirely), and the
# warp field need not exist in the rebound handler — evaluating against the
# captured pair is what keeps `warp |> Gradient |> adaptive plot` working.
struct Deformation{DH<:Ferrite.AbstractDofHandler,UO<:Makie.Observable,SO<:Makie.Observable}
    dh::DH
    u::UO                             # the warped stage's dof vector
    field::Makie.Observable{Symbol}   # resolved at apply time; may be switched
    scale::SO
end

# A dof field captured together with the handler and solution it lives on (a
# later filter may rebind both), as the leaf of a derivation chain.
struct FieldSource{DH<:Ferrite.AbstractDofHandler,UO<:Makie.Observable}
    dh::DH
    u::UO
    name::Symbol
end

# The symbolic side of a pointwise derivation (VonMises, Derive, ...): the
# closure and its inputs, recorded alongside the sampled array the filter
# registers as point data. The array knows the static tessellation's vertices
# and nothing else; the record is what lets the adaptive path re-evaluate the
# quantity at arbitrary (cell, ξ) — and, following the regularity principle
# (a derived quantity is at most as regular as its sources), what names the
# dof fields the refinement criterion should sample. Inputs are `FieldSource`
# leaves or nested `DerivedPointData` (untyped vector: the chain is walked
# once at plot creation, never per vertex).
struct DerivedPointData{F}
    f::F
    inputs::Vector{Any}
end

# The quadrature-point partition of an `AddQuadraturePointData` stage: the
# rule (or refshape → rule mapping), the live values, and how a stored entry
# maps to the plotted value. The filter bakes the partition's *static*
# tessellation into the dataset's geometry; this record is what lets the
# adaptive path rebuild the partition as its base domain — each Voronoi
# region fanned so the piecewise-constant regions survive refinement — and
# look values up at arbitrary refined vertices.
struct QPPartition{Q,V<:Makie.Observable,F}
    qr::Q
    values::V
    extract::F
    output::Symbol
end

"""
    Adaptivity(; geometry_tol=1e-3, solution_tol=5e-3, max_depth=10, sample_type=Float32)

How a dataset's plots refine — the error-adaptive tessellation settings,
owned by the [`FEData`](@ref) (every plot of a dataset follows the same
settings; there are no per-plot overrides). The tolerances and the depth cap
are Observables: assigning them (`ds.adaptivity.solution_tol[] = 1e-4`)
re-refines every open plot of the dataset.

`Adaptivity` is a filter: `ds |> Adaptivity(...)` returns the same dataset
(geometry and data shared) with these settings — configuring a dataset that
had adaptivity disabled, or replacing another configuration. The `FEData`
constructor applies it automatically: `adaptivity=true` (the default) with
these defaults, `adaptivity=Adaptivity(...)` with yours, `adaptivity=false`
for the static tessellation.

- `geometry_tol`: geometry-error tolerance, as a fraction of the grid's
  bounding-box diagonal. The drawn triangles approximate the exact geometry
  (dofhandler interpolation, including warps) to within it.
- `solution_tol`: solution-error tolerance, as a fraction of the color
  field's value span. The linear vertex-color interpolation approximates the
  exact field polynomial to within it.
- `max_depth`: maximum bisection depth per base triangle.
- `sample_type`: the number type the pipeline samples geometry and fields in,
  by default what the renderer draws (GLMakie uploads Float32). Tolerances
  are floored at its resolution; pass `Float64` to sample at full precision.

Construct `FEData(dh, u; adaptivity=Adaptivity(...))` to tweak,
`adaptivity=false` to disable and always draw the static tessellation.
"""
struct Adaptivity <: AbstractFilter
    geometry_tol::Makie.Observable{Float64}
    solution_tol::Makie.Observable{Float64}
    max_depth::Makie.Observable{Int}
    sample_type::DataType
end

function Adaptivity(; geometry_tol::Real=1e-3, solution_tol::Real=5e-3, max_depth::Integer=10,
                    sample_type::Type{<:AbstractFloat}=Float32)
    return Adaptivity(Makie.Observable(Float64(geometry_tol)),
                      Makie.Observable(Float64(solution_tol)),
                      Makie.Observable(Int(max_depth)), sample_type)
end

# The constructor's `adaptivity` keyword: a config, `true` for the defaults,
# or `false`/`nothing` to disable.
_adaptivity_config(a::Adaptivity) = a
_adaptivity_config(a::Bool) = a ? Adaptivity() : nothing
_adaptivity_config(::Nothing) = nothing

"""
    FEData(dh::Ferrite.AbstractDofHandler, u::Vector; topology, adaptivity=true)

Source node of the visualization pipeline: builds the static "L2" triangulation
of `Ferrite.get_grid(dh)` (nodes shared between cells are duplicated per cell so
discontinuous fields render with their jumps) and holds `u` as an Observable
for live updating via [`FerriteViz.update!`](@ref).

Named data arrays are resolved with [`point_data`](@ref)/[`cell_data`](@ref)
and registered with [`set_point_data!`](@ref)/[`set_cell_data!`](@ref).
Transformations are applied by piping into filters:
`ds |> WarpByVector(:u, 2.0) |> Gradient(:u) |> VonMises()`.

For large 3D grids, pass a precomputed `topology::Ferrite.ExclusiveTopology`
to avoid rebuilding it.

`adaptivity` controls how this dataset's plots refine: `adaptivity=true`
(the default) applies the [`Adaptivity`](@ref) filter with its default
tolerances, so `solutionplot` and `meshplot` re-tessellate the visible cells
by longest-edge bisection until the drawn triangles resolve both the exact
geometry and the color field — watertight, camera-independent, following
[`FerriteViz.update!`](@ref). Pass an `Adaptivity(...)` of your own to
tweak the tolerances (equivalent to `FEData(...; adaptivity=false) |>
Adaptivity(...)`), or `adaptivity=false` to always draw the static
tessellation. The setting is a dataset property, shared by every plot of it
and carried through filters. A color the adaptive path cannot re-evaluate
at refined vertices (a raw point-data array) falls back to the static
tessellation for that plot.

The *static* tessellation the constructor builds is the flat base for every
cell — with adaptivity on, curved rendering comes from the adaptive path.
For a uniformly subdivided static tessellation instead, compose explicitly:
`FEData(dh, u; adaptivity=false) |> Refine(2)` (or [`Refine`](@ref)`()` for
its automatic per-cell-type choice).
"""
struct FEData{dim,DH<:Ferrite.AbstractDofHandler,T1,TOP<:Union{Nothing,Ferrite.AbstractTopology},SU<:Makie.Observable,M,TRI,S} <: AbstractPlotter
    dh::DH
    u::Makie.Observable{Vector{T1}}   # this dataset's dof vector (possibly lifted from source_u)
    source_u::SU                      # the root solution observable; update! target
    topology::TOP
    visible::Vector{Bool}             # per-cell visibility; immutable after construction
    gridnodes::Makie.Observable{Vector{GeometryBasics.Point{dim,Float32}}}  # the grid's nodes (meshplot)
    # Coordinates of the tessellation vertices, i.e. the vertices of the
    # rendered triangulation, not the finite element cells' vertices. The
    # Observable is what a warp lifts from; the Buffer is what the GPU sees.
    coords::Makie.Observable{Vector{GeometryBasics.Point{dim,Float32}}}
    coords_buffer::ShaderAbstractions.Buffer{GeometryBasics.Point{dim,Float32},Vector{GeometryBasics.Point{dim,Float32}}}
    all_triangles::Vector{TRI}                        # every triangle of the tessellation
    vis_triangles::ShaderAbstractions.Buffer{TRI,Vector{TRI}}  # the subset actually drawn (see CrinkleClip)
    triangle_cell_map::Vector{Int}    # triangle -> owning cell
    cell_triangle_offsets::Vector{Int}  # cell -> range in all_triangles (see triangles_on_cell)
    cell_vertex_offsets::Vector{Int}    # cell -> range in coords (see vertices_on_cell)
    # Decomposition of each cell's volume into simplices (tets in 3D, triangles
    # in 2D; S is NTuple{dim+1,Int}), indexing into the same vertex array as the
    # triangles, so point data covers the simplex vertices too. Only cells whose
    # reference dimension equals the spatial dimension carry simplices (embedded
    # shells/lines have no volume). This is what Clip cuts and
    # ExtractIsosurfaces marches; datasets without volume (e.g. an extracted
    # isosurface) have it empty.
    simplices::Vector{S}
    simplex_cell_map::Vector{Int}       # simplex -> owning cell
    cell_simplex_offsets::Vector{Int}   # cell -> range in simplices (see simplices_on_cell)
    # Wireframe segments along the FE cell edges, as pairs of indices into
    # coords. Their endpoints are ordinary tessellation vertices, which is what
    # makes the meshplot wireframe follow warps, solutions and clips for free.
    all_edges::Vector{NTuple{2,Int}}
    edge_cell_map::Vector{Int}          # edge segment -> owning cell
    cell_edge_offsets::Vector{Int}      # cell -> range in all_edges (see edges_on_cell)
    reference_coords::Matrix{Float64}   # per tessellation vertex, in the cell's reference coordinates
    mesh::M                             # coords_buffer + vis_triangles, handed to Makie as is
    point_data::Dict{Symbol,Makie.Observable}  # arrays on the tessellation vertices
    cell_data::Dict{Symbol,Makie.Observable}   # arrays on the cells
    # Derivation provenance for point-data arrays that have one (see
    # `DerivedPointData`); copied and cleared exactly alongside `point_data`.
    point_derivations::Dict{Symbol,DerivedPointData}
    # Quadrature-point partition provenance (see `QPPartition`); `nothing`
    # unless an `AddQuadraturePointData` produced this dataset. Survives
    # geometry-preserving filters (a CrinkleClip of internal variables is the
    # main consumer), dropped by geometry-rebuilding ones.
    qp_partition::Union{Nothing,QPPartition}
    # Deformation provenance: one record per WarpByVector applied upstream, in
    # application order. The displaced coordinates are baked into `coords`, but
    # consumers that need the deformation as a *continuous* function of the
    # reference coordinate (adaptive tessellation evaluates geometry at
    # arbitrary ξ) reconstruct it from this record.
    deformation::Vector{Deformation}
    # Which cells make up the body, as opposed to `visible`, which marks the
    # cells contributing surface. In 3D an interior cell is solid but not
    # visible; a cell removed by CrinkleClip is neither (visible implies
    # solid). The distinction is what identifies the surface facets — those
    # whose neighbour is missing or not solid. Consumed only by the adaptive
    # path, through `_is_surface_facet`, which both base builds
    # (`_isubd_base_cells`, `_isubd_base_qp`) gate their 3D facet loops on;
    # the static path deliberately ignores it and draws every facet of every
    # visible cell. Clip filters remove cells by clearing it, and
    # `transfer_solution` gates on it, so hidden (but solid) interior cells
    # carry real values — which is what lets the volume-based filters (Clip,
    # ExtractIsosurfaces) and warps see the field everywhere.
    solid::Vector{Bool}
    # True while every solid cell's geometry is the full tessellation of the
    # cell. Exact cuts (Clip, ExtractIsosurfaces) clear it, which blocks filters
    # that rebuild whole cells from the grid (Refine, AddQuadraturePointData)
    # from resurrecting cut-away geometry, and makes plots of the dataset fall
    # back to the static path (the adaptive base is fanned from whole cells).
    cells_intact::Bool
    # How this dataset's plots refine (`Adaptivity`), or `nothing` for the
    # static tessellation. A dataset property — the substrate below is shared
    # by every adaptive plot, so there are no per-plot overrides — carried
    # through filters *by reference*: a filtered dataset shares the config,
    # so one knob steers the whole pipeline family. Untyped concerns
    # (sample_type as a plain DataType) stay behind the `_substrate` function
    # barrier.
    adaptivity::Union{Nothing,Adaptivity}
    # Lazily built adaptive-tessellation substrate (`IsubdSubstrate`), shared
    # by every adaptive plot of this dataset; `nothing` until the first one
    # asks. Everything in it is a pure function of the fields above, so it is
    # never invalidated — filters return new FEData instances, each with a
    # fresh (empty) cache. Untyped on purpose: plots retrieve it through the
    # `_substrate` function barrier, keeping this struct's parameters stable.
    subd_cache::Base.RefValue{Any}
end

# Only volumetric 3D grids get a default topology (used to hide interior
# cells): ExclusiveTopology does not support embedded cells, and in 2D
# everything is visible anyway.
function _default_topology(grid)
    Ferrite.getspatialdim(grid) > 2 || return nothing
    all(c -> Ferrite.getrefdim(Ferrite.getrefshape(c)) == 3, Ferrite.getcells(grid)) || return nothing
    return Ferrite.ExclusiveTopology(grid)
end

# `:default` is the sentinel every filter and representation resolves to the
# *first* field of the dof handler (see `_resolve_name`). A dof field actually
# named `default` could therefore never be addressed: naming it would silently
# select the first field instead, which is the wrong array whenever `default`
# is not itself in first position. Reject it where it enters the package.
function _check_reserved_fieldnames(dh::Ferrite.AbstractDofHandler)
    :default in Ferrite.getfieldnames(dh) &&
        error("`:default` is reserved by FerriteViz: it names the first field of the dof handler " *
              "wherever a field, point-data or cell-data name is expected. Rename the dof field " *
              "`:default` to address it explicitly.")
    return nothing
end

function FEData(dh::Ferrite.AbstractDofHandler, u::AbstractVector;
                topology=_default_topology(Ferrite.get_grid(dh)), adaptivity=true)
    # copy: update! writes into this array and must not mutate the caller's u
    return FEData(dh, Makie.Observable(collect(u)); topology, adaptivity)
end

function FEData(dh::Ferrite.AbstractDofHandler, u::Makie.Observable;
                topology=_default_topology(Ferrite.get_grid(dh)), source_u::Makie.Observable=u,
                adaptivity=true)
    _check_reserved_fieldnames(dh)
    grid = Ferrite.get_grid(dh)
    sdim = Ferrite.getspatialdim(grid)
    ncells = Ferrite.getncells(grid)

    visible = zeros(Bool, ncells)
    if sdim > 2 && topology !== nothing
        boundaryfaces = findall(isempty, topology.face_face_neighbor)
        visible[Ferrite.getindex.(boundaryfaces, 1)] .= true
    else
        # without a topology (2D, or embedded cells) everything is drawn
        visible .= true
    end

    # The static tessellation is the flat base for every cell: with
    # adaptivity on (the default) the curved rendering comes from the
    # error-adaptive path, and uniform static subdivision is an explicit
    # composition — `FEData(dh, u; adaptivity=false) |> Refine(n)` (or
    # `Refine()` for the automatic per-cell-type choice).
    return _build_dataset(dh, u, source_u, topology, visible,
                          _tessellation_provider(Refine(0), dh);
                          adaptivity=_adaptivity_config(adaptivity))
end

# Volume simplices exist only for cells whose reference dimension matches the
# spatial one (embedded shells/lines have no volume to decompose). In 2D the
# tessellation triangles already tile the cell; in 3D the surface triangles are
# fanned into tets from one extra centroid vertex per cell (valid because the
# reference shapes are convex, hence star-shaped).
_cell_has_volume(cell, sdim::Int) = sdim >= 2 && Ferrite.getrefdim(Ferrite.getrefshape(cell)) == sdim

# Shared tessellation-instantiation core of the FEData constructor and the
# Refine filter: lay out `tess_for(cell)` per cell with duplicated vertices.
function _build_dataset(dh::Ferrite.AbstractDofHandler, u::Makie.Observable, source_u::Makie.Observable,
                        topology, visible::Vector{Bool}, tess_for;
                        adaptivity::Union{Nothing,Adaptivity}=Adaptivity(),
                        solid::Vector{Bool}=fill(true, Ferrite.getncells(Ferrite.get_grid(dh))))
    grid = Ferrite.get_grid(dh)
    cells = Ferrite.getcells(grid)
    sdim = Ferrite.getspatialdim(grid)
    ncells = length(cells)

    cell_triangle_offsets = Vector{Int}(undef, ncells + 1)
    cell_vertex_offsets = Vector{Int}(undef, ncells + 1)
    cell_edge_offsets = Vector{Int}(undef, ncells + 1)
    cell_simplex_offsets = Vector{Int}(undef, ncells + 1)
    cell_triangle_offsets[1] = 0
    cell_vertex_offsets[1] = 0
    cell_edge_offsets[1] = 0
    cell_simplex_offsets[1] = 0
    for (i, cell) in enumerate(cells)
        tess = tess_for(cell)
        hasvol = _cell_has_volume(cell, sdim)
        cell_triangle_offsets[i+1] = cell_triangle_offsets[i] + ntriangles(tess)
        # 3D volume cells carry one extra vertex, the fan centroid
        cell_vertex_offsets[i+1] = cell_vertex_offsets[i] + nvertices(tess) + (hasvol && sdim == 3 ? 1 : 0)
        cell_edge_offsets[i+1] = cell_edge_offsets[i] + nedges(tess)
        cell_simplex_offsets[i+1] = cell_simplex_offsets[i] + (hasvol ? ntriangles(tess) : 0)
    end
    num_triangles = cell_triangle_offsets[end]
    num_verts = cell_vertex_offsets[end]
    num_edges = cell_edge_offsets[end]

    triangles = Matrix{Int}(undef, num_triangles, 3)
    triangle_cell_map = Vector{Int}(undef, num_triangles)
    all_edges = Vector{NTuple{2,Int}}(undef, num_edges)
    edge_cell_map = Vector{Int}(undef, num_edges)
    physical_coords = Vector{GeometryBasics.Point{sdim,Float32}}(undef, num_verts)
    reference_coords = zeros(Float64, num_verts, sdim)
    S = NTuple{sdim + 1,Int}
    simplices = Vector{S}(undef, cell_simplex_offsets[end])
    simplex_cell_map = Vector{Int}(undef, cell_simplex_offsets[end])

    for (cell_id, cell) in enumerate(cells)
        # Function barrier: `tess_for` and `geometric_interpolation` are only
        # abstractly inferred here (the tessellation cache is heterogeneous),
        # so instantiate through a call specialized on the concrete types —
        # one dynamic dispatch per cell instead of per tessellation vertex.
        _instantiate_cell!(physical_coords, reference_coords, triangles, triangle_cell_map,
                           all_edges, edge_cell_map, simplices, simplex_cell_map,
                           _cell_has_volume(cell, sdim), tess_for(cell),
                           Ferrite.geometric_interpolation(typeof(cell)),
                           Ferrite.getcoordinates(grid, cell_id),
                           cell_vertex_offsets[cell_id], cell_triangle_offsets[cell_id],
                           cell_edge_offsets[cell_id], cell_simplex_offsets[cell_id], cell_id)
    end

    # convert: to_triangles yields an untyped empty vector for 0 triangles
    all_triangles = convert(Vector{GeometryBasics.GLTriangleFace}, Makie.to_triangles(triangles))
    vis_triangles = ShaderAbstractions.Buffer(Makie.Observable(_visibility_triangles(all_triangles, visible, triangle_cell_map)))
    coords = Makie.Observable(physical_coords)
    coords_buffer = ShaderAbstractions.Buffer(coords)
    mesh = GeometryBasics.Mesh(coords_buffer, vis_triangles)
    gridnodes = Makie.Observable([GeometryBasics.Point{sdim,Float32}(Ferrite.get_node_coordinate(n)...) for n in Ferrite.getnodes(grid)])

    return FEData{sdim,typeof(dh),eltype(u[]),typeof(topology),typeof(source_u),typeof(mesh),eltype(all_triangles),S}(
        dh, u, source_u, topology, visible, gridnodes, coords, coords_buffer,
        all_triangles, vis_triangles, triangle_cell_map, cell_triangle_offsets,
        cell_vertex_offsets, simplices, simplex_cell_map, cell_simplex_offsets,
        all_edges, edge_cell_map, cell_edge_offsets,
        reference_coords, mesh,
        Dict{Symbol,Makie.Observable}(), Dict{Symbol,Makie.Observable}(),
        Dict{Symbol,DerivedPointData}(), nothing,
        Deformation[], solid, true, adaptivity, Ref{Any}(nothing))
end

function _instantiate_cell!(physical_coords::Vector{GeometryBasics.Point{sdim,Float32}}, reference_coords,
                            triangles, triangle_cell_map, all_edges, edge_cell_map,
                            simplices, simplex_cell_map, hasvol::Bool,
                            tess::ReferenceTessellation, ip_geo::Ferrite.ScalarInterpolation,
                            node_coords::AbstractVector, coff::Int, toff::Int, eoff::Int,
                            soff::Int, cell_id::Int) where {sdim}
    for (k, ξ) in enumerate(tess.coords)
        x = geometric_map(ip_geo, node_coords, ξ)
        physical_coords[coff+k] = GeometryBasics.Point{sdim,Float32}(x...)
        for d in 1:length(ξ)
            reference_coords[coff+k, d] = ξ[d]
        end
    end
    for (t, tri) in enumerate(tess.triangles)
        for j in 1:3
            triangles[toff+t, j] = tri[j] + coff
        end
        triangle_cell_map[toff+t] = cell_id
    end
    for (e, edge) in enumerate(tess.edges)
        all_edges[eoff+e] = (edge[1] + coff, edge[2] + coff)
        edge_cell_map[eoff+e] = cell_id
    end
    if hasvol
        if sdim == 2
            for (t, tri) in enumerate(tess.triangles)
                simplices[soff+t] = (tri[1] + coff, tri[2] + coff, tri[3] + coff)
                simplex_cell_map[soff+t] = cell_id
            end
        elseif sdim == 3
            center = coff + nvertices(tess) + 1
            # The fan apex is the affine mean of the cell's tessellation
            # vertices in both reference and physical space (the same equal
            # weights), not the geometric map of the reference centroid: the
            # fan tiles the *linearized* cell (the flat surface triangles are
            # what is rendered and cut), and the mapped reference centroid of
            # a strongly curved cell can fall outside that polyhedron, which
            # would invert part of the fan. The mean of the mapped vertices is
            # inside their convex hull, so the fan of a convex cell is valid.
            ξc = sum(tess.coords) / nvertices(tess)
            xc = zero(Tensors.Vec{3,Float64})
            for k in 1:nvertices(tess)
                xc += Tensors.Vec{3,Float64}(NTuple{3,Float64}(physical_coords[coff+k]))
            end
            physical_coords[center] = GeometryBasics.Point{sdim,Float32}((xc / nvertices(tess))...)
            for d in 1:length(ξc)
                reference_coords[center, d] = ξc[d]
            end
            # tets are stored positively oriented; the surface triangles are
            # consistently oriented per cell, so one shared flip decides the
            # whole fan — taken from the *total* signed fan volume (the
            # enclosed volume up to sign), which stays reliable when
            # individual triangles are degenerate
            cpos = Tensors.Vec{3,Float64}(NTuple{3,Float64}(physical_coords[center]))
            signed = 0.0
            for tri in tess.triangles
                a = Tensors.Vec{3,Float64}(NTuple{3,Float64}(physical_coords[tri[1] + coff]))
                b = Tensors.Vec{3,Float64}(NTuple{3,Float64}(physical_coords[tri[2] + coff]))
                c = Tensors.Vec{3,Float64}(NTuple{3,Float64}(physical_coords[tri[3] + coff]))
                signed += _signed_tet_volume(a, b, c, cpos)
            end
            flip = signed < 0
            for (t, tri) in enumerate(tess.triangles)
                simplices[soff+t] = flip ? (tri[2] + coff, tri[1] + coff, tri[3] + coff, center) :
                                           (tri[1] + coff, tri[2] + coff, tri[3] + coff, center)
                simplex_cell_map[soff+t] = cell_id
            end
        end
    end
    return nothing
end

function _visibility_triangles(all_triangles, visible, triangle_cell_map)
    vis_triangles = copy(all_triangles)
    for (i, cell_id) in enumerate(triangle_cell_map)
        if !visible[cell_id]
            vis_triangles[i] = GeometryBasics.GLTriangleFace(1, 1, 1)
        end
    end
    return vis_triangles
end

# Derive a new FEData from `ds`, sharing every field that is not overridden.
# This helper owns the dependent-field invariants so no filter can leave the
# GPU-facing objects pointing at stale geometry: passing `coords` rebuilds the
# coordinate buffer, passing any of `all_triangles`/`visible`/
# `triangle_cell_map` rebuilds the visibility-filtered triangle buffer, and a
# rebuild of either rebuilds the mesh. The substrate cache is always fresh
# (filters return new FEData instances; the cache is a pure function of them).
function _derive(ds::FEData{dim};
                 dh=ds.dh,
                 u::Makie.Observable=ds.u,
                 visible::Union{Nothing,Vector{Bool}}=nothing,
                 gridnodes::Makie.Observable=ds.gridnodes,
                 coords::Union{Nothing,Makie.Observable}=nothing,
                 all_triangles::Union{Nothing,Vector}=nothing,
                 triangle_cell_map::Vector{Int}=ds.triangle_cell_map,
                 cell_triangle_offsets::Vector{Int}=ds.cell_triangle_offsets,
                 cell_vertex_offsets::Vector{Int}=ds.cell_vertex_offsets,
                 simplices::Vector=ds.simplices,
                 simplex_cell_map::Vector{Int}=ds.simplex_cell_map,
                 cell_simplex_offsets::Vector{Int}=ds.cell_simplex_offsets,
                 all_edges::Vector{NTuple{2,Int}}=ds.all_edges,
                 edge_cell_map::Vector{Int}=ds.edge_cell_map,
                 cell_edge_offsets::Vector{Int}=ds.cell_edge_offsets,
                 reference_coords::Matrix{Float64}=ds.reference_coords,
                 point_data::Dict{Symbol,Makie.Observable}=copy(ds.point_data),
                 cell_data::Dict{Symbol,Makie.Observable}=copy(ds.cell_data),
                 point_derivations::Dict{Symbol,DerivedPointData}=copy(ds.point_derivations),
                 qp_partition::Union{Nothing,QPPartition}=ds.qp_partition,
                 deformation::Vector{Deformation}=ds.deformation,
                 solid::Vector{Bool}=ds.solid,
                 cells_intact::Bool=ds.cells_intact,
                 adaptivity::Union{Nothing,Adaptivity}=ds.adaptivity) where {dim}
    new_coords = coords !== nothing
    new_coords || (coords = ds.coords)
    coords_buffer = new_coords ? ShaderAbstractions.Buffer(coords) : ds.coords_buffer
    new_triangles = all_triangles !== nothing || visible !== nothing
    all_triangles === nothing && (all_triangles = ds.all_triangles)
    visible === nothing && (visible = ds.visible)
    vis_triangles = new_triangles ?
        ShaderAbstractions.Buffer(Makie.Observable(_visibility_triangles(all_triangles, visible, triangle_cell_map))) :
        ds.vis_triangles
    mesh = (new_coords || new_triangles) ? GeometryBasics.Mesh(coords_buffer, vis_triangles) : ds.mesh
    return FEData{dim,typeof(dh),eltype(u[]),typeof(ds.topology),typeof(ds.source_u),typeof(mesh),eltype(all_triangles),eltype(simplices)}(
        dh, u, ds.source_u, ds.topology, visible, gridnodes, coords, coords_buffer,
        all_triangles, vis_triangles, triangle_cell_map, cell_triangle_offsets, cell_vertex_offsets,
        simplices, simplex_cell_map, cell_simplex_offsets, all_edges, edge_cell_map, cell_edge_offsets,
        reference_coords, mesh, point_data, cell_data, point_derivations, qp_partition,
        deformation, solid, cells_intact, adaptivity, Ref{Any}(nothing))
end

# The point-data arrays that survive a geometry rebuild or visibility change:
# registered arrays are plain per-vertex values, but keys naming dof fields are
# cached lazy transfers whose values depend on the dataset they were computed
# on — carrying them over would make results depend on whether somebody
# accessed the field upstream. They are dropped and re-resolve lazily.
function _registered_point_data(ds::FEData)
    fieldnames = Ferrite.getfieldnames(ds.dh)
    return Dict{Symbol,Makie.Observable}(k => v for (k, v) in ds.point_data if !(k in fieldnames))
end

"""
Total number of tessellation vertices, i.e. vertices of the rendered
triangulation. These are not the vertices of the finite element cells: cells do
not share them (they are duplicated per cell, so discontinuities render), and a
tessellated cell generally carries more of them than it has corners.
"""
num_vertices(ds::FEData) = length(ds.coords[])

vertices_on_cell(ds::FEData, cell_idx::Int) = (ds.cell_vertex_offsets[cell_idx]+1):ds.cell_vertex_offsets[cell_idx+1]
triangles_on_cell(ds::FEData, cell_idx::Int) = (ds.cell_triangle_offsets[cell_idx]+1):ds.cell_triangle_offsets[cell_idx+1]
edges_on_cell(ds::FEData, cell_idx::Int) = (ds.cell_edge_offsets[cell_idx]+1):ds.cell_edge_offsets[cell_idx+1]
simplices_on_cell(ds::FEData, cell_idx::Int) = (ds.cell_simplex_offsets[cell_idx]+1):ds.cell_simplex_offsets[cell_idx+1]

# Flat vertex-index list (2 entries per segment) of the wireframe of the
# visible cells. Static per dataset (visibility is immutable after
# construction), so meshplot computes it once and per-frame work is only the
# coordinate gather.
function _visible_edge_indices(ds::FEData)
    indices = Int[]
    for (e, cell_id) in enumerate(ds.edge_cell_map)
        ds.visible[cell_id] || continue
        edge = ds.all_edges[e]
        push!(indices, edge[1], edge[2])
    end
    return indices
end

# Grid nodes belonging to at least one visible cell (meshplot's node markers
# and labels follow clips like the surface does).
function _visible_node_ids(ds::FEData)
    grid = Ferrite.get_grid(ds.dh)
    mask = falses(Ferrite.getnnodes(grid))
    for (cell_id, cell) in enumerate(Ferrite.getcells(grid))
        ds.visible[cell_id] || continue
        for n in cell.nodes
            mask[n] = true
        end
    end
    return findall(mask)
end

"""
    FerriteViz.update!(ds::FEData, u::Vector)

Update the source solution observable, propagating through all filters and open
plots down to the GPU buffers. Can be called on any dataset of a pipeline; it
always updates the root solution (`u` must match its length).

Not exported: both Makie and Ferrite export distinct functions named `update!`,
so call this one qualified.
"""
function update!(ds::FEData, u::Vector)
    length(ds.source_u[]) == length(u) || error("length mismatch: source solution has $(length(ds.source_u[])) dofs, got $(length(u))")
    ds.source_u[] .= u
    Makie.notify(ds.source_u)
    return nothing
end

##############
# Named data #
##############

function default_field(ds::FEData)
    fieldnames = Ferrite.getfieldnames(ds.dh)
    isempty(fieldnames) && error("the DofHandler has no fields, specify data by name")
    return first(fieldnames)
end

_resolve_name(ds::FEData, name::Symbol) = name === :default ? default_field(ds) : name

function _data_association(ds::FEData, name::Symbol)
    name in Ferrite.getfieldnames(ds.dh) && return :point
    haskey(ds.point_data, name) && return :point
    haskey(ds.cell_data, name) && return :cell
    return :none
end

_available_data(ds::FEData) =
    "fields $(collect(Ferrite.getfieldnames(ds.dh))), point data $(collect(keys(ds.point_data))), cell data $(collect(keys(ds.cell_data)))"

"""
    point_data(ds::FEData, name::Symbol) -> Observable{Matrix{Float64}}

The named data array on the tessellation vertices (nvertices × ncomponents;
tensor components in Tensors.jl linear order). These are the vertices of the
triangulation the dataset renders, not the vertices of the finite element cells
— see the [architecture overview](@ref "Architecture"). Fields of the dof handler
are transferred to the tessellation lazily and cached; `:default` resolves to
the first field.

!!! note
    `:default` is reserved for this purpose wherever a name is expected, so a
    dof handler carrying a field named `:default` is rejected by [`FEData`](@ref).
"""
function point_data(ds::FEData, name::Symbol)
    name = _resolve_name(ds, name)
    haskey(ds.point_data, name) && return ds.point_data[name]
    if name in Ferrite.getfieldnames(ds.dh)
        obs = Makie.lift(u -> transfer_solution(ds, u; field_name=name), ds.u)
        ds.point_data[name] = obs
        return obs
    end
    error("no point data named :$name; available: $(_available_data(ds))")
end

"""
    cell_data(ds::FEData, name::Symbol) -> Observable

The named per-cell data array (length ncells), registered via [`set_cell_data!`](@ref).
"""
function cell_data(ds::FEData, name::Symbol)
    haskey(ds.cell_data, name) && return ds.cell_data[name]
    error("no cell data named :$name; available: $(_available_data(ds))")
end

_canonical_point_data(A::AbstractVector) = reshape(convert(Vector{Float64}, A), :, 1)
_canonical_point_data(A::AbstractMatrix) = convert(Matrix{Float64}, A)

# Named data may not shadow a dof field: field names always resolve to the
# dof-backed arrays (e.g. WarpByVector's grid-node path relies on this).
function _check_no_field_shadow(ds::FEData, name::Symbol)
    name in Ferrite.getfieldnames(ds.dh) &&
        error("cannot register data named :$name, it would shadow the dof field of the same name; pick another name")
    return nothing
end

"""
    set_point_data!(ds::FEData, name::Symbol, data)

Register a named data array on the tessellation vertices (see
[`num_vertices`](@ref FerriteViz.num_vertices) — one row per vertex of the
rendered triangulation, *not* per grid node). `data` may be a `Vector`/`Matrix`
or an `Observable` of one — updates to a registered Observable propagate into
plots. Existing names are overwritten; dof field names cannot be shadowed.
"""
function set_point_data!(ds::FEData, name::Symbol, data::AbstractVecOrMat)
    _check_no_field_shadow(ds, name)
    size(data, 1) == num_vertices(ds) || error("point data must have $(num_vertices(ds)) rows, got $(size(data, 1))")
    ds.point_data[name] = Makie.Observable(_canonical_point_data(data))
    return ds
end
function set_point_data!(ds::FEData, name::Symbol, data::Makie.Observable)
    _check_no_field_shadow(ds, name)
    ds.point_data[name] = Makie.lift(data) do A
        size(A, 1) == num_vertices(ds) || error("point data must have $(num_vertices(ds)) rows, got $(size(A, 1))")
        _canonical_point_data(A)
    end
    return ds
end

"""
    set_cell_data!(ds::FEData, name::Symbol, data)

Register a named per-cell data array (length ncells, any element type — e.g.
stress tensors, to be reduced by filters like [`VonMises`](@ref)). `data` may
be a `Vector` or an `Observable` of one. Existing names are overwritten; dof
field names cannot be shadowed.
"""
function set_cell_data!(ds::FEData, name::Symbol, data::AbstractVector)
    _check_no_field_shadow(ds, name)
    ncells = Ferrite.getncells(Ferrite.get_grid(ds.dh))
    length(data) == ncells || error("cell data must have $ncells entries, got $(length(data))")
    ds.cell_data[name] = Makie.Observable(collect(data))
    return ds
end
function set_cell_data!(ds::FEData, name::Symbol, data::Makie.Observable)
    _check_no_field_shadow(ds, name)
    ncells = Ferrite.getncells(Ferrite.get_grid(ds.dh))
    ds.cell_data[name] = Makie.lift(data) do v
        length(v) == ncells || error("cell data must have $ncells entries, got $(length(v))")
        v
    end
    return ds
end

# NaN out the rows of invisible cells (mutating `x`, which must not alias a
# cached array): interior vertices carry real values since transfer gates on
# `solid` (the volume filters need them), but colors — and Makie's automatic
# colorrange — follow the *drawn* geometry, as they did before the volume
# decomposition existed.
function _mask_invisible!(ds::FEData, x::Vector{Float64})
    for cell in 1:length(ds.visible)
        ds.visible[cell] && continue
        for v in vertices_on_cell(ds, cell)
            x[v] = NaN
        end
    end
    return x
end

# Scalar per-vertex Observable for coloring: point data must have one component
# (except when resolving :default, where a vector field falls back to its
# magnitude), cell data (scalar) is expanded to the tessellation vertices.
function _scalar_data(ds::FEData, name::Symbol; reduce_default::Bool=false)
    name = _resolve_name(ds, name)
    assoc = _data_association(ds, name)
    assoc === :none && error("no data named :$name; available: $(_available_data(ds))")
    if assoc === :cell
        return Makie.lift(v -> transfer_scalar_celldata(ds, v), cell_data(ds, name))
    end
    all_visible = all(ds.visible)
    return Makie.lift(point_data(ds, name)) do A
        if size(A, 2) == 1
            all_visible ? vec(A) : _mask_invisible!(ds, Float64.(vec(A)))
        elseif reduce_default
            mags = [LinearAlgebra.norm(view(A, i, :)) for i in 1:size(A, 1)]
            all_visible ? mags : _mask_invisible!(ds, mags)
        else
            error("point data :$name has $(size(A, 2)) components; reduce it to a scalar first, e.g. with Magnitude(input=:$name) or ExtractComponent(i; input=:$name)")
        end
    end
end

#################
# Data transfer #
#################

# TODO upstream to Ferrite.jl
function getsubdofhandlers(dh::Ferrite.DofHandler, field_name::Symbol)
    sdhs = SubDofHandler[]
    for sdh in dh.subdofhandlers
        if field_name ∈ sdh.field_names
            push!(sdhs, sdh)
        end
    end
    return sdhs
end

"""
    transfer_solution(ds::FEData, u::Vector; field_name=:u) -> Matrix{Float64}

Evaluate the field at every tessellation vertex of every solid cell from the
owning element's dofs (preserving inter-element discontinuities). Vertices of
removed cells or cells outside the field's subdomain stay `NaN` — hidden (but
solid) interior cells are evaluated, so volume-based filters (`Clip`,
`ExtractIsosurfaces`) and warps see real values everywhere.
"""
function transfer_solution(ds::FEData, u::Vector; field_name::Symbol=:u)
    dh = ds.dh
    sdhs = getsubdofhandlers(dh, field_name)
    isempty(sdhs) && error("field :$field_name not found in the DofHandler")
    ip_field = Ferrite.getfieldinterpolation(first(sdhs), field_name)
    ref_dim = Ferrite.getrefdim(ip_field)
    ξ0 = Ferrite.Vec(ntuple(d -> 0.0, ref_dim))
    # NOTE this does not work for ansatz spaces where derivatives are mixed in (e.g. Hermite)
    ncomps = length(Ferrite.reference_shape_value(ip_field, ξ0, 1))
    data = fill(NaN, num_vertices(ds), ncomps)
    for sdh in sdhs
        ip_field = Ferrite.getfieldinterpolation(sdh, field_name)
        ip_geo = Ferrite.geometric_interpolation(Ferrite.getcelltype(sdh))
        pv = Ferrite.PointValues(ip_field, ip_geo; update_gradients=false)
        # function barrier for pv and the reference dimension
        _transfer_solution!(data, pv, sdh, field_name, ds, u, Val(Ferrite.getrefdim(ip_field)))
    end
    return data
end

function _transfer_solution!(data, pv, sdh, field_name::Symbol, ds::FEData, u::Vector, ::Val{refdim}) where {refdim}
    dh = ds.dh
    grid = Ferrite.get_grid(dh)
    cellset = collect(sdh.cellset)
    local_dof_range = Ferrite.dof_range(sdh, field_name)
    ncomps = size(data, 2)

    local_coords = Ferrite.getcoordinates(grid, first(cellset))
    local_celldofs = Ferrite.celldofs(dh, first(cellset))
    for cell_idx in cellset
        ds.solid[cell_idx] || continue
        Ferrite.getcoordinates!(local_coords, grid, cell_idx)
        Ferrite.celldofs!(local_celldofs, dh, cell_idx)
        celldofs_field = @view(local_celldofs[local_dof_range])
        for v in vertices_on_cell(ds, cell_idx)
            ξ = Tensors.Vec{refdim}(d -> ds.reference_coords[v, d])
            Ferrite.reinit!(pv, local_coords, ξ)
            val = Ferrite.function_value(pv, 1, @views(u[celldofs_field]))
            for d in 1:ncomps
                data[v, d] = val[d]
            end
        end
    end
    return data
end

"""
    transfer_scalar_celldata(ds::FEData, values::AbstractVector) -> Vector{Float64}

Expand one scalar per cell to the tessellation vertices.
"""
function transfer_scalar_celldata(ds::FEData, values::AbstractVector)
    ncells = Ferrite.getncells(Ferrite.get_grid(ds.dh))
    length(values) == ncells || error("expected one value per cell ($ncells), got $(length(values))")
    eltype(values) <: Real || error("cell data must be scalar for plotting; reduce it first, e.g. with VonMises() or Magnitude()")
    data = Vector{Float64}(undef, num_vertices(ds))
    for cell_idx in 1:ncells
        data[vertices_on_cell(ds, cell_idx)] .= Float64(values[cell_idx])
    end
    return data
end
