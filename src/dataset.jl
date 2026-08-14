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

"""
    FEData(dh::Ferrite.AbstractDofHandler, u::Vector; topology, adaptive=true)

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

`adaptive=true` (the default) tessellates with the [`Refine`](@ref) filter's
automatic choice: cell types whose geometry and fields are all (multi-)linear
keep the flat base tessellation, everything else is subdivided so curved and
high-order-deformed cells render curved — at the price of more triangles (see
[`Refine`](@ref) for the numbers). Opt out with `adaptive=false` (flat base
tessellation for every cell); custom levels are a filter application:
`FEData(dh, u; adaptive=false) |> Refine(2)`.

!!! note
    The tessellation `adaptive=true` picks may change in a future release; such
    a change is breaking. `adaptive=false` and explicit `Refine(n)` counts are
    stable.
"""
struct FEData{dim,DH<:Ferrite.AbstractDofHandler,T1,TOP<:Union{Nothing,Ferrite.AbstractTopology},SU<:Makie.Observable,M,TRI} <: AbstractPlotter
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
    # Deformation provenance: one record per WarpByVector applied upstream, in
    # application order. The displaced coordinates are baked into `coords`, but
    # consumers that need the deformation as a *continuous* function of the
    # reference coordinate (adaptive tessellation evaluates geometry at
    # arbitrary ξ) reconstruct it from this record.
    deformation::Vector{Deformation}
    # Which cells make up the body, as opposed to `visible`, which marks the
    # cells contributing surface. In 3D an interior cell is solid but not
    # visible; a cell removed by CrinkleClip is neither. The distinction is
    # what identifies the surface facets — those whose neighbour is missing or
    # not solid — for consumers that extract the boundary surface themselves.
    solid::Vector{Bool}
    # Lazily built adaptive-tessellation substrate (`IsubdSubstrate`), shared
    # by every adaptive plot of this dataset; `nothing` until the first one
    # asks. Everything in it is a pure function of the fields above, so it is
    # never invalidated — filters return new FEData instances, each with a
    # fresh (empty) cache. Untyped on purpose: plots retrieve it through the
    # `_substrate` function barrier, keeping this struct's parameters stable.
    subd_cache::Base.RefValue{Any}
end

function _default_topology(grid)
    return Ferrite.getspatialdim(grid) > 2 ? Ferrite.ExclusiveTopology(grid) : nothing
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
                topology=_default_topology(Ferrite.get_grid(dh)), adaptive::Bool=true)
    # copy: update! writes into this array and must not mutate the caller's u
    return FEData(dh, Makie.Observable(collect(u)); topology, adaptive)
end

function FEData(dh::Ferrite.AbstractDofHandler, u::Makie.Observable;
                topology=_default_topology(Ferrite.get_grid(dh)), source_u::Makie.Observable=u,
                adaptive::Bool=true)
    _check_reserved_fieldnames(dh)
    grid = Ferrite.get_grid(dh)
    sdim = Ferrite.getspatialdim(grid)
    ncells = Ferrite.getncells(grid)

    visible = zeros(Bool, ncells)
    if sdim > 2
        boundaryfaces = findall(isempty, topology.face_face_neighbor)
        visible[Ferrite.getindex.(boundaryfaces, 1)] .= true
    else
        visible .= true
    end

    # The tessellation choice is the Refine filter's; the constructor merely
    # applies its automatic mode by default — Refine() picks the subdivision
    # per cell type (see _pick_subdivision_rounds), Refine(0) pins every cell
    # to the flat base tessellation. Building through the provider directly
    # means the default costs nothing over constructing flat and filtering
    # afterwards.
    refinement = adaptive ? Refine() : Refine(0)
    return _build_dataset(dh, u, source_u, topology, visible, _tessellation_provider(refinement, dh))
end

# Shared tessellation-instantiation core of the FEData constructor and the
# Refine filter: lay out `tess_for(cell)` per cell with duplicated vertices.
function _build_dataset(dh::Ferrite.AbstractDofHandler, u::Makie.Observable, source_u::Makie.Observable,
                        topology, visible::Vector{Bool}, tess_for)
    grid = Ferrite.get_grid(dh)
    cells = Ferrite.getcells(grid)
    sdim = Ferrite.getspatialdim(grid)
    ncells = length(cells)

    cell_triangle_offsets = Vector{Int}(undef, ncells + 1)
    cell_vertex_offsets = Vector{Int}(undef, ncells + 1)
    cell_edge_offsets = Vector{Int}(undef, ncells + 1)
    cell_triangle_offsets[1] = 0
    cell_vertex_offsets[1] = 0
    cell_edge_offsets[1] = 0
    for (i, cell) in enumerate(cells)
        tess = tess_for(cell)
        cell_triangle_offsets[i+1] = cell_triangle_offsets[i] + ntriangles(tess)
        cell_vertex_offsets[i+1] = cell_vertex_offsets[i] + nvertices(tess)
        cell_edge_offsets[i+1] = cell_edge_offsets[i] + nedges(tess)
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

    for (cell_id, cell) in enumerate(cells)
        # Function barrier: `tess_for` and `geometric_interpolation` are only
        # abstractly inferred here (the tessellation cache is heterogeneous),
        # so instantiate through a call specialized on the concrete types —
        # one dynamic dispatch per cell instead of per tessellation vertex.
        _instantiate_cell!(physical_coords, reference_coords, triangles, triangle_cell_map,
                           all_edges, edge_cell_map, tess_for(cell),
                           Ferrite.geometric_interpolation(typeof(cell)),
                           Ferrite.getcoordinates(grid, cell_id),
                           cell_vertex_offsets[cell_id], cell_triangle_offsets[cell_id],
                           cell_edge_offsets[cell_id], cell_id)
    end

    # convert: to_triangles yields an untyped empty vector for 0 triangles
    all_triangles = convert(Vector{GeometryBasics.GLTriangleFace}, Makie.to_triangles(triangles))
    vis_triangles = ShaderAbstractions.Buffer(Makie.Observable(_visibility_triangles(all_triangles, visible, triangle_cell_map)))
    coords = Makie.Observable(physical_coords)
    coords_buffer = ShaderAbstractions.Buffer(coords)
    mesh = GeometryBasics.Mesh(coords_buffer, vis_triangles)
    gridnodes = Makie.Observable([GeometryBasics.Point{sdim,Float32}(Ferrite.get_node_coordinate(n)...) for n in Ferrite.getnodes(grid)])

    return FEData{sdim,typeof(dh),eltype(u[]),typeof(topology),typeof(source_u),typeof(mesh),eltype(all_triangles)}(
        dh, u, source_u, topology, visible, gridnodes, coords, coords_buffer,
        all_triangles, vis_triangles, triangle_cell_map, cell_triangle_offsets,
        cell_vertex_offsets, all_edges, edge_cell_map, cell_edge_offsets,
        reference_coords, mesh,
        Dict{Symbol,Makie.Observable}(), Dict{Symbol,Makie.Observable}(),
        Dict{Symbol,DerivedPointData}(),
        Deformation[], fill(true, ncells), Ref{Any}(nothing))
end

function _instantiate_cell!(physical_coords::Vector{GeometryBasics.Point{sdim,Float32}}, reference_coords,
                            triangles, triangle_cell_map, all_edges, edge_cell_map,
                            tess::ReferenceTessellation, ip_geo::Ferrite.ScalarInterpolation,
                            node_coords::AbstractVector, coff::Int, toff::Int, eoff::Int,
                            cell_id::Int) where {sdim}
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
    return Makie.lift(point_data(ds, name)) do A
        if size(A, 2) == 1
            vec(A)
        elseif reduce_default
            [LinearAlgebra.norm(view(A, i, :)) for i in 1:size(A, 1)]
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

Evaluate the field at every tessellation vertex of every visible cell from the
owning element's dofs (preserving inter-element discontinuities). Vertices of
invisible cells or cells outside the field's subdomain stay `NaN`.
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
        ds.visible[cell_idx] || continue
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
