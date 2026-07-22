# Layer 2a: the FEData source.
#
# FEData wraps a dof handler and a solution Observable together with the static
# triangulation built from the tessellation interface. Coordinates and triangle
# indices live in ShaderAbstractions.Buffers shared into a GeometryBasics.Mesh,
# so downstream Observable updates mutate the GPU data in place without
# rebuilding the plot.

"""
    FEData(dh::Ferrite.AbstractDofHandler, u::Vector; topology)

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
"""
struct FEData{dim,DH<:Ferrite.AbstractDofHandler,T1,TOP<:Union{Nothing,Ferrite.AbstractTopology},SU<:Makie.Observable,M,TRI} <: AbstractPlotter
    dh::DH
    u::Makie.Observable{Vector{T1}}   # this dataset's dof vector (possibly lifted from source_u)
    source_u::SU                      # the root solution observable; update! target
    topology::TOP
    visible::Vector{Bool}             # per-cell visibility; immutable after construction
    gridnodes::Makie.Observable{Vector{GeometryBasics.Point{dim,Float32}}}
    coords::Makie.Observable{Vector{GeometryBasics.Point{dim,Float32}}}  # tessellation vertex coords
    coords_buffer::ShaderAbstractions.Buffer{GeometryBasics.Point{dim,Float32},Vector{GeometryBasics.Point{dim,Float32}}}
    all_triangles::Vector{TRI}
    vis_triangles::ShaderAbstractions.Buffer{TRI,Vector{TRI}}
    triangle_cell_map::Vector{Int}
    cell_triangle_offsets::Vector{Int}
    cell_vertex_offsets::Vector{Int}
    reference_coords::Matrix{Float64} # per-vertex reference coordinates (padded to spatial dim)
    mesh::M
    point_data::Dict{Symbol,Makie.Observable}
    cell_data::Dict{Symbol,Makie.Observable}
end

function _default_topology(grid)
    return Ferrite.getspatialdim(grid) > 2 ? Ferrite.ExclusiveTopology(grid) : nothing
end

function FEData(dh::Ferrite.AbstractDofHandler, u::AbstractVector;
                topology=_default_topology(Ferrite.get_grid(dh)))
    # copy: update! writes into this array and must not mutate the caller's u
    return FEData(dh, Makie.Observable(collect(u)); topology)
end

function FEData(dh::Ferrite.AbstractDofHandler, u::Makie.Observable;
                topology=_default_topology(Ferrite.get_grid(dh)), source_u::Makie.Observable=u)
    grid = Ferrite.get_grid(dh)
    cells = Ferrite.getcells(grid)
    sdim = Ferrite.getspatialdim(grid)
    ncells = length(cells)

    visible = zeros(Bool, ncells)
    if sdim > 2
        boundaryfaces = findall(isempty, topology.face_face_neighbor)
        visible[Ferrite.getindex.(boundaryfaces, 1)] .= true
    else
        visible .= true
    end

    tess_cache = Dict{Type,ReferenceTessellation}()
    tess_for(cell) = get!(() -> reference_tessellation(Ferrite.getrefshape(cell)), tess_cache, Ferrite.getrefshape(cell))

    cell_triangle_offsets = Vector{Int}(undef, ncells + 1)
    cell_vertex_offsets = Vector{Int}(undef, ncells + 1)
    cell_triangle_offsets[1] = 0
    cell_vertex_offsets[1] = 0
    for (i, cell) in enumerate(cells)
        tess = tess_for(cell)
        cell_triangle_offsets[i+1] = cell_triangle_offsets[i] + ntriangles(tess)
        cell_vertex_offsets[i+1] = cell_vertex_offsets[i] + nvertices(tess)
    end
    num_triangles = cell_triangle_offsets[end]
    num_verts = cell_vertex_offsets[end]

    triangles = Matrix{Int}(undef, num_triangles, 3)
    triangle_cell_map = Vector{Int}(undef, num_triangles)
    physical_coords = Vector{GeometryBasics.Point{sdim,Float32}}(undef, num_verts)
    reference_coords = zeros(Float64, num_verts, sdim)

    for (cell_id, cell) in enumerate(cells)
        tess = tess_for(cell)
        ip_geo = Ferrite.geometric_interpolation(typeof(cell))
        node_coords = Ferrite.getcoordinates(grid, cell_id)
        coff = cell_vertex_offsets[cell_id]
        for (k, ξ) in enumerate(tess.coords)
            x = geometric_map(ip_geo, node_coords, ξ)
            physical_coords[coff+k] = GeometryBasics.Point{sdim,Float32}(x...)
            for d in 1:length(ξ)
                reference_coords[coff+k, d] = ξ[d]
            end
        end
        toff = cell_triangle_offsets[cell_id]
        for (t, tri) in enumerate(tess.triangles)
            for j in 1:3
                triangles[toff+t, j] = tri[j] + coff
            end
            triangle_cell_map[toff+t] = cell_id
        end
    end

    all_triangles = Makie.to_triangles(triangles)
    vis_triangles = ShaderAbstractions.Buffer(Makie.Observable(_visibility_triangles(all_triangles, visible, triangle_cell_map)))
    coords = Makie.Observable(physical_coords)
    coords_buffer = ShaderAbstractions.Buffer(coords)
    mesh = GeometryBasics.Mesh(coords_buffer, vis_triangles)
    gridnodes = Makie.Observable([GeometryBasics.Point{sdim,Float32}(Ferrite.get_node_coordinate(n)...) for n in Ferrite.getnodes(grid)])

    return FEData{sdim,typeof(dh),eltype(u[]),typeof(topology),typeof(source_u),typeof(mesh),eltype(all_triangles)}(
        dh, u, source_u, topology, visible, gridnodes, coords, coords_buffer,
        all_triangles, vis_triangles, triangle_cell_map, cell_triangle_offsets,
        cell_vertex_offsets, reference_coords, mesh,
        Dict{Symbol,Makie.Observable}(), Dict{Symbol,Makie.Observable}())
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
Total number of (duplicated) tessellation vertices.
"""
num_vertices(ds::FEData) = length(ds.coords[])

vertices_on_cell(ds::FEData, cell_idx::Int) = (ds.cell_vertex_offsets[cell_idx]+1):ds.cell_vertex_offsets[cell_idx+1]
triangles_on_cell(ds::FEData, cell_idx::Int) = (ds.cell_triangle_offsets[cell_idx]+1):ds.cell_triangle_offsets[cell_idx+1]

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

The named per-vertex data array (nvertices × ncomponents; tensor components in
Tensors.jl linear order). Fields of the dof handler are transferred to the
tessellation lazily and cached; `:default` resolves to the first field.
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

"""
    set_point_data!(ds::FEData, name::Symbol, data)

Register a named per-vertex data array. `data` may be a `Vector`/`Matrix`
(nvertices rows) or an `Observable` of one — updates to a registered Observable
propagate into plots. Existing names are overwritten.
"""
function set_point_data!(ds::FEData, name::Symbol, data::AbstractVecOrMat)
    size(data, 1) == num_vertices(ds) || error("point data must have $(num_vertices(ds)) rows, got $(size(data, 1))")
    ds.point_data[name] = Makie.Observable(_canonical_point_data(data))
    return ds
end
function set_point_data!(ds::FEData, name::Symbol, data::Makie.Observable)
    size(data[], 1) == num_vertices(ds) || error("point data must have $(num_vertices(ds)) rows, got $(size(data[], 1))")
    ds.point_data[name] = Makie.lift(_canonical_point_data, data)
    return ds
end

"""
    set_cell_data!(ds::FEData, name::Symbol, data)

Register a named per-cell data array (length ncells, any element type — e.g.
stress tensors, to be reduced by filters like [`VonMises`](@ref)). `data` may
be a `Vector` or an `Observable` of one. Existing names are overwritten.
"""
function set_cell_data!(ds::FEData, name::Symbol, data::AbstractVector)
    ncells = Ferrite.getncells(Ferrite.get_grid(ds.dh))
    length(data) == ncells || error("cell data must have $ncells entries, got $(length(data))")
    ds.cell_data[name] = Makie.Observable(collect(data))
    return ds
end
function set_cell_data!(ds::FEData, name::Symbol, data::Makie.Observable)
    ncells = Ferrite.getncells(Ferrite.get_grid(ds.dh))
    length(data[]) == ncells || error("cell data must have $ncells entries, got $(length(data[]))")
    ds.cell_data[name] = data
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
            error("point data :$name has $(size(A, 2)) components; reduce it to a scalar first, e.g. with Magnitude(input=:$name) or Component(i; input=:$name)")
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
