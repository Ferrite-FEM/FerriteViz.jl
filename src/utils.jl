# Helper... Refactoring needed.
function getfieldinterpolation(dh::Ferrite.DofHandler, field_name::Symbol)
    field_idx = indexin(dh.field_names, [:b])
    field_idx == nothing && error("did not find field $field_name")
    dh.field_interpolations[field_idx]
end

function field_offset(dh::Ferrite.DofHandler, field_idx::Int)
    offset = 0
    for i in 1:(field_idx-1)
        offset += Ferrite.getnbasefunctions(dh.field_interpolations[i])::Int * dh.field_dims[i]
    end
    return offset
end

function facedofs(cell, local_face_idx::Int, field_idx::Int)
    celldofs_ = Ferrite.celldofs(cell)
    offset_ = field_offset(cell.dh, field_idx)
    fip = Ferrite.getfieldinterpolation(cell.dh, field_idx)
    faces_ = Ferrite.faces(fip)
    isempty(faces_) && return () # no face dofs :)
    indices = faces_[local_face_idx] .+ offset_
    celldofs_[[indices...]]
end

Ferrite.vertices(cell::Ferrite.Cell{3,3,1}) = cell.nodes

Ferrite.default_interpolation(::Type{Ferrite.Cell{3,3,1}}) = Ferrite.Lagrange{2,Ferrite.RefTetrahedron,1}()

# Note: This extracts the face spanned by the vertices, not the actual face!
linear_face_cell(cell::Ferrite.Cell{3,N,4}, local_face_idx::Int) where N =  Ferrite.Cell{3,3,1}(Ferrite.faces(cell)[local_face_idx])
linear_face_cell(cell::Ferrite.Cell{3,N,6}, local_face_idx::Int) where N = Ferrite.Quadrilateral3D(Ferrite.faces(cell)[local_face_idx])

# Obtain the face interpolation on regular geometries.
getfaceip(ip::Ferrite.Interpolation{dim, shape, order}, local_face_idx::Int) where {dim, shape <: Union{Ferrite.RefTetrahedron, Ferrite.RefCube}, order} = Ferrite.getlowerdim(ip)

struct MakiePlotter{dim,DH<:Ferrite.AbstractDofHandler,T} <: AbstractPlotter
    dh::DH
    u::Makie.Observable{Vector{T}} # original solution on the original mesh (i.e. dh.mesh)

    gridnodes::Matrix{T} # coordinates of grid nodes in matrix form
    physical_coords::Matrix{T} # coordinates in physical space of a vertex
    triangles::Matrix{Int} # each row carries a triple with the indices into the coords matrix defining the triangle
    triangle_cell_map::Vector{Int}
    reference_coords::Matrix{T} # coordinates on the associated reference cell for the corresponding triangle vertex
end

"""
Build a static triangulation of the grid to render out the solution via Makie.
"""
function MakiePlotter(dh::Ferrite.AbstractDofHandler, u::Vector)
    cells = Ferrite.getcells(dh.grid)
    dim = Ferrite.getdim(dh.grid)

    # We do not take any assumptions on the mesh, so first we have to loopkup
    num_triangles = 0
    for cell in cells
        num_triangles += ntriangles(cell)
    end

    # Preallocate the matrices carrying the triangulation information
    triangles = Matrix{Int}(undef, num_triangles, 3)
    triangle_cell_map = Vector{Int}(undef, num_triangles)
    physical_coords = Matrix{Float64}(undef, num_triangles*3, dim)
    gridnodes = [node[i] for node in Ferrite.getcoordinates.(Ferrite.getnodes(dh.grid)), i in 1:dim]
    reference_coords = Matrix{Float64}(undef, num_triangles*3, dim) # NOTE this should have the dimension of the actual reference element

    # Decompose does the heavy lifting for us
    coord_offset = 1
    triangle_offset = 1
    for (cell_id,cell) ∈ enumerate(cells)
        triangle_offset_begin = triangle_offset
        (coord_offset, triangle_offset) = decompose!(coord_offset, physical_coords, reference_coords, triangle_offset, triangles, dh.grid, cell)
        triangle_cell_map[triangle_offset_begin:(triangle_offset-1)] .= cell_id
    end
    return MakiePlotter{dim,typeof(dh),eltype(u)}(dh,Observable(u),gridnodes,physical_coords,triangles,triangle_cell_map,reference_coords);
end

"""
Clip plane described by the normal and its distance to the coordinate origin.
"""
struct ClipPlane
    normal::Tensors.Vec{3}
    distance::Real
end

"""
Binary decision function to clip a cell with a plane for the crincle clip.
"""
function (plane::ClipPlane)(grid, cellid)
    cell = grid.cells[cellid]
    coords = Ferrite.getcoordinates.(Ferrite.getnodes(grid)[[cell.nodes...]])
    for coord ∈ coords
        if coord ⋅ plane.normal > plane.distance
            return false
        end
    end
    return true
end

"""
Crincle clip generates a new plotter that deletes some of the triangles, based on an
implicit description of the clipping surface. Here `decision_fun` takes the grid and
a cell index as input and returns whether the cell is visible or not.
"""
function crincle_clip(plotter::MakiePlotter{3,DH,T}, decision_fun) where {DH,T}
    dh = plotter.dh
    u = plotter.u
    grid = dh.grid

    # We iterate over all triangles and check if the corresponding cell is visible.
    visible_triangles = Vector{Bool}(undef, size(plotter.triangles, 1))
    visible_coords = Vector{Bool}(undef, 3*size(plotter.triangles, 1))
    for (i, triangle) ∈ enumerate(eachrow(plotter.triangles))
        cell_id = plotter.triangle_cell_map[i]
        visible_triangles[i] = decision_fun(grid, cell_id)
        visible_coords[3*(i-1)+1] = visible_coords[3*(i-1)+2] = visible_coords[3*(i-1)+3] = visible_triangles[i]
    end

    # Create a plotter with views on the data.
    return MakiePlotter{3,DH,T}(dh, u, plotter.gridnodes,
        plotter.physical_coords,
        plotter.triangles[visible_triangles, :],
        plotter.triangle_cell_map[visible_triangles],
        plotter.reference_coords);
end

"""
Total number of vertices
"""
num_vertices(p::MakiePlotter) = size(p.physical_coords,1)

#Visualization is just fancy triangle plotting, every element needs to be translatable to a triangle *sadnoises*
to_triangle(cell::Ferrite.AbstractCell{2,N,3}) where N = [Ferrite.vertices(cell)[1:3]]
to_triangle(cell::Ferrite.AbstractCell{3,N,4}) where N = [Ferrite.vertices(cell)[[1,2,3]], Ferrite.vertices(cell)[[1,2,4]], Ferrite.vertices(cell)[[2,3,4]], Ferrite.vertices(cell)[[1,4,3]]]
to_triangle(cell::Ferrite.AbstractCell{2,N,4}) where N = [Ferrite.vertices(cell)[[1,2,3]], Ferrite.vertices(cell)[[3,4,1]]]
to_triangle(cell::Ferrite.AbstractCell{3,N,6}) where N = [Ferrite.vertices(cell)[[1,2,3]], Ferrite.vertices(cell)[[3,4,1]],
                                                          Ferrite.vertices(cell)[[1,5,6]], Ferrite.vertices(cell)[[6,2,1]],
                                                          Ferrite.vertices(cell)[[2,6,7]], Ferrite.vertices(cell)[[7,3,2]],
                                                          Ferrite.vertices(cell)[[3,7,8]], Ferrite.vertices(cell)[[8,4,3]],
                                                          Ferrite.vertices(cell)[[1,4,8]], Ferrite.vertices(cell)[[8,5,1]],
                                                          Ferrite.vertices(cell)[[5,8,7]], Ferrite.vertices(cell)[[7,6,5]]]

"""
TODO this looks faulty...think harder.
"""
# Helper to count triangles e.g. for preallocations.
ntriangles(cell::Ferrite.AbstractCell{2,3,3}) = 1 # Tris in 2D
ntriangles(cell::Ferrite.AbstractCell{3,3,1}) = 1 # Tris in 3D
ntriangles(cell::Ferrite.AbstractCell{dim,N,4}) where {dim,N} = 4 # Quads in 2D and 3D
ntriangles(cell::Ferrite.AbstractCell{3,N,1}) where N = 4 # Tets as a special case of a Quad, obviously :)
ntriangles(cell::Ferrite.AbstractCell{3,N,6}) where N = 6*4 # Hex

"""
Get the vertices represented as a list of coordinates of a cell.

@TODO refactor into Ferrite core.
"""
function vertices(grid::Ferrite.AbstractGrid, cell::Ferrite.AbstractCell{dim,N,M}) where {dim,N,M}
    Ferrite.getnodes(grid)[[Ferrite.vertices(cell)...]]
end

"""
Decompose a triangle into a coordinates and a triangle index list to disconnect it properly. Guarantees to preserve orderings and orientations.
"""
function decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Union{Ferrite.AbstractCell{2,N,3}, Ferrite.AbstractCell{3,3,1}}) where {N}
    for (i,v) in enumerate(vertices(grid, cell))
        coord_matrix[coord_offset, :] = Ferrite.getcoordinates(v)
        ref_coord_matrix[coord_offset, 1:2] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i]
        triangle_matrix[triangle_offset, i] = coord_offset
        coord_offset+=1
    end
    triangle_offset+=1
    (coord_offset, triangle_offset)
end

"""
Decompose a quadrilateral into a coordinates and a triangle index list to disconnect it properly. Guarantees to preserve orderings and orientations.

Assumes a CCW quadrilateral numbering:
4---3
|   |
1---2
and creates the decomposition
4-------3
| \\ C / |
|  \\ /  |
|D  5  B|
|  / \\  |
| / A \\ |
1-------2
where A=(1,2,5),B=(2,3,5),C=(3,4,5),D=(4,1,5) are the generated triangles in this order.
"""
function decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Union{Ferrite.AbstractCell{2,N,4}, Ferrite.AbstractCell{3,4,1}}) where {N}
    # A bit more complicated. The default diagonal decomposition into 2 triangles is missing a solution mode.
    # To resolve this we make a more expensive decomposition in 4 triangles which correctly displays the solution mode (in linear problems).
    coord_offset_initial = coord_offset
    vts = vertices(grid, cell)
    space_dim = size(coord_matrix,2)

    # Compute coordinate of vertex 5
    center = [ntuple(x->0.0, space_dim)...]
    for v in vts
        center += Ferrite.getcoordinates(v)
    end
    center /= 4.0

    # Generate triangles in order
    for i = 1:length(vts)
        v1 = vts[i]
        # next index on the outer chain CCW
        i2 = i+1
        if i2 > length(vts)
            i2 = 1 # dunno how to modulo this correctly :/
        end
        v2 = vts[i2]

        # current vertex
        coord_matrix[coord_offset, :] = Ferrite.getcoordinates(v1)
        ref_coord_matrix[coord_offset, 1:2] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i]
        coord_offset+=1

        # next vertex in chain
        coord_matrix[coord_offset, :] = Ferrite.getcoordinates(v2)
        ref_coord_matrix[coord_offset, 1:2] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i2]
        coord_offset+=1

        # center vertex (5)
        coord_matrix[coord_offset, :] = center
        ref_coord_matrix[coord_offset, 1:2] = [ntuple(x->0.0, 2)...]
        coord_offset+=1

        # collect indices
        triangle_matrix[triangle_offset, :] = (0:2) .+ (coord_offset_initial+3*(i-1))
        triangle_offset+=1
    end

    (coord_offset, triangle_offset)
end

"""
Decompose volumetric objects via their faces.
"""
function decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Ferrite.AbstractCell{3,N,M}) where {N,M}
    # Just 6 quadrilaterals :)
    for face_index ∈ 1:M
        face_coord_offset = coord_offset
        (coord_offset, triangle_offset) = decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, linear_face_cell(cell, face_index))
        for ci ∈ face_coord_offset:(coord_offset-1)
            new_coord = transfer_quadrature_face_to_cell(ref_coord_matrix[ci, 1:2], cell, face_index)
            ref_coord_matrix[ci, :] = new_coord
        end
    end
    (coord_offset, triangle_offset)
end


refshape(cell::Ferrite.AbstractCell) = typeof(Ferrite.default_interpolation(typeof(cell))).parameters[2]

x₁(x) = x[1]
x₂(x) = x[2]
x₃(x) = x[3]
l2(x) = LinearAlgebra.norm(x,2)
l1(x) = LinearAlgebra.norm(x,1)

midpoint(cell::Ferrite.AbstractCell{2,N,3}, points) where N = Point2f((1/3) * (points[cell.nodes[1],:] + points[cell.nodes[2],:] + points[cell.nodes[3],:]))
midpoint(cell::Ferrite.AbstractCell{2,N,4}, points) where N = Point2f(0.5 * (points[cell.nodes[1],:] + points[cell.nodes[3],:]))
midpoint(cell::Ferrite.AbstractCell{3,N,4}, points) where N = Point3f((1/4) * (points[cell.nodes[1],:] + points[cell.nodes[2],:] + points[cell.nodes[3],:] + points[cell.nodes[4],:]))
midpoint(cell::Ferrite.AbstractCell{3,N,6}, points) where N = Point3f(0.5 * (points[cell.nodes[1],:] + points[cell.nodes[7],:]))

function postprocess(node_values)
    dim = length(node_values)
    if dim == 1
        return node_values
    else
        return sqrt(sum(node_values.^2))
    end
end

"""
Transfer the solution of a plotter to the new mesh in 2D.

@precondition: Assumes that the number of basis function for each dof is equal.

@TODO Refactor. This is peak inefficiency.
"""
function transfer_solution(plotter::MakiePlotter{2}, u::Vector; field_idx::Int=1, process::Function=FerriteViz.postprocess)
    n_vertices_per_tri = 3 # we have 3 vertices per triangle...

    # select objects from plotter
    dh = plotter.dh
    ref_coords = plotter.reference_coords
    grid = dh.grid

    # field related variables
    field_name = Ferrite.getfieldnames(dh)[field_idx]
    field_dim = Ferrite.getfielddim(dh, field_idx)
    ip_field = dh.field_interpolations[field_idx]

    # actual data
    local_dof_range = Ferrite.dof_range(dh, field_name)

    cell_geo_ref = Ferrite.getcells(grid, 1)
    ip_geo = Ferrite.default_interpolation(typeof(cell_geo_ref))
    pv = (field_dim == 1) ? Ferrite.PointScalarValues(ip_field, ip_geo) : Ferrite.PointVectorValues(ip_field, ip_geo)

    data = fill(0.0, num_vertices(plotter), field_dim)
    current_vertex_index = 1
    for (cell_index, cell_geo) in enumerate(Ferrite.getcells(grid))
        _celldofs_field = reshape(Ferrite.celldofs(dh,cell_index)[local_dof_range], (field_dim, Ferrite.getnbasefunctions(ip_field)))

        # Loop over all local triangle vertices
        for i in 1:(ntriangles(cell_geo)*n_vertices_per_tri)
            ξ = Tensors.Vec{2}(ref_coords[current_vertex_index, :])
            for d in 1:field_dim
                for node_idx ∈ 1:Ferrite.getnbasefunctions(ip_field)
                    data[current_vertex_index, d] += Ferrite.value(ip_field, node_idx, ξ) ⋅ u[_celldofs_field[d, node_idx]]
                end
            end
            current_vertex_index += 1
        end
    end
    return mapslices(process, data, dims=[2])
end

"""
Transfer the solution of a plotter to the new mesh in 3D.

@TODO Refactor. This is peak inefficiency.
"""
function transfer_solution(plotter::MakiePlotter{3}, u::Vector; field_idx::Int=1, process::Function=FerriteViz.postprocess)
    n_vertices_per_tri = 3 # we have 3 vertices per triangle...

    # select objects from plotter
    dh = plotter.dh
    ref_coords = plotter.reference_coords
    grid = dh.grid

    # field related variables
    field_name = Ferrite.getfieldnames(dh)[field_idx]
    field_dim = Ferrite.getfielddim(dh, field_idx)

    # TODO this should be moved inside the loop below to gt the correct interpolator for the current cell.
    ip_field = dh.field_interpolations[field_idx]

    # actual data
    local_dof_range = Ferrite.dof_range(dh, field_name)

    cell_geo_ref = Ferrite.getcells(grid, 1)
    ip_geo = Ferrite.default_interpolation(typeof(cell_geo_ref))
    pv = (field_dim == 1) ? Ferrite.PointScalarValues(ip_field, ip_geo) : Ferrite.PointVectorValues(ip_field, ip_geo)

    current_vertex_index = 1
    data = fill(0.0, num_vertices(plotter), field_dim)
    for cell in Ferrite.CellIterator(plotter.dh)
        cell_idx = Ferrite.cellid(cell)
        cell_geo = Ferrite.getcells(grid, cell_idx)
        # This should make the loop work for mixed grids
        if typeof(cell_geo_ref) != typeof(cell_geo)
            ip_geo = Ferrite.default_interpolation(typeof(cell_geo))
            pv = (field_dim == 1) ? Ferrite.PointScalarValues(ip_field, ip_geo) : Ferrite.PointVectorValues(ip_field, ip_geo)
            cell_geo_ref = cell_geo
        end
        _local_celldofs = Ferrite.celldofs(cell)[local_dof_range]
        _celldofs_field = reshape(_local_celldofs, (field_dim, Ferrite.getnbasefunctions(ip_field)))
        # TODO replace this with a triangle-to-cell map.
        for (local_face_idx,_) in enumerate(Ferrite.faces(cell_geo))
            # Construct face values to evaluate
            face_geo = linear_face_cell(cell_geo, local_face_idx)
            # TODO Optimize for mixed geometries
            nfvertices = ntriangles(face_geo)*n_vertices_per_tri

            # Loop over vertices
            for i in 1:nfvertices
                Ferrite.reinit!(pv, Ferrite.getcoordinates(grid,cell_idx), Tensors.Vec(ref_coords[current_vertex_index,:]...))
                # val = Ferrite.function_value(cv, i, u[_local_celldofs])
                val = Ferrite.function_value(pv, 1, u[_local_celldofs])
                for d in 1:field_dim
                    data[current_vertex_index, d] += val[d] #current_vertex_index
                end
                current_vertex_index += 1
            end
        end
    end

    return mapslices(process, data, dims=[2])
end

function transfer_scalar_celldata(plotter::MakiePlotter{3}, u::Vector; process::Function=FerriteViz.postprocess)
    n_vertices = 3 # we have 3 vertices per triangle...

    # select objects from plotter
    dh = plotter.dh
    grid = dh.grid

    current_vertex_index = 1
    data = fill(0.0, num_vertices(plotter), 1)
    for (cell_index, cell_geo) in enumerate(Ferrite.getcells(grid))
        for (local_face_idx,_) in enumerate(Ferrite.faces(cell_geo))
            face_geo = linear_face_cell(cell_geo, local_face_idx)
            # Loop over vertices
            for i in 1:(ntriangles(face_geo)*n_vertices)
                data[current_vertex_index, 1] = u[cell_index]
                current_vertex_index += 1
            end
        end
    end

    return mapslices(process, data, dims=[2])
end

function transfer_scalar_celldata(plotter::MakiePlotter{2}, u::Vector;  process::Function=FerriteViz.postprocess)
    n_vertices = 3 # we have 3 vertices per triangle...

    # select objects from plotter
    dh = plotter.dh
    grid = dh.grid

    current_vertex_index = 1
    data = fill(0.0, num_vertices(plotter), 1)
    for (cell_index, cell_geo) in enumerate(Ferrite.getcells(grid))
        for i in 1:(ntriangles(cell_geo)*n_vertices)
            data[current_vertex_index, 1] = u[cell_index]
            current_vertex_index += 1
        end
    end

    return mapslices(process, data, dims=[2])
end

function transfer_scalar_celldata(grid::Ferrite.AbstractGrid{3}, num_vertices::Number, u::Vector; process::Function=FerriteViz.postprocess)
    n_vertices = 3 # we have 3 vertices per triangle...
    current_vertex_index = 1
    data = fill(0.0, num_vertices, 1)
    for (cell_index, cell_geo) in enumerate(Ferrite.getcells(grid))
        for (local_face_idx,_) in enumerate(Ferrite.faces(cell_geo))
            face_geo = linear_face_cell(cell_geo, local_face_idx)
            # Loop over vertices
            for i in 1:(ntriangles(face_geo)*n_vertices)
                data[current_vertex_index, 1] = u[cell_index]
                current_vertex_index += 1
            end
        end
    end
    return mapslices(process, data, dims=[2])
end

function transfer_scalar_celldata(grid::Ferrite.AbstractGrid{2}, num_vertices::Number, u::Vector;  process::Function=FerriteViz.postprocess)
    n_vertices = 3 # we have 3 vertices per triangle...
    current_vertex_index = 1
    data = fill(0.0, num_vertices, 1)
    for (cell_index, cell_geo) in enumerate(Ferrite.getcells(grid))
        for i in 1:(ntriangles(cell_geo)*n_vertices)
            data[current_vertex_index, 1] = u[cell_index]
            current_vertex_index += 1
        end
    end
    return mapslices(process, data, dims=[2])
end

function dof_to_node(dh::Ferrite.AbstractDofHandler, u::Array{T,1}; field::Int=1, process::Function=postprocess) where T
    fieldnames = Ferrite.getfieldnames(dh)
    field_dim = Ferrite.getfielddim(dh, field)
    data = fill(NaN, Ferrite.getnnodes(dh.grid), field_dim)
    offset = Ferrite.field_offset(dh, fieldnames[field])

    for cell in Ferrite.CellIterator(dh)
        _celldofs = Ferrite.celldofs(cell)
        counter = 1
        for node in cell.nodes
            for d in 1:field_dim
                data[node, d] = u[_celldofs[counter + offset]]
                counter += 1
            end
        end
    end
    return mapslices(process, data, dims=[2])
end

get_gradient_interpolation(::Ferrite.Lagrange{dim,shape,order}) where {dim,shape,order} = Ferrite.DiscontinuousLagrange{dim,shape,order-1}()
get_gradient_interpolation_type(::Type{Ferrite.Lagrange{dim,shape,order}}) where {dim,shape,order} = Ferrite.DiscontinuousLagrange{dim,shape,order-1}
# TODO remove if Knuth's PR on this gets merged (Ferrite PR 552)
getgrid(dh::Ferrite.DofHandler) = dh.grid

function ε(x::Vector{T}) where T
    ngrad = length(x)
    dim = isqrt(ngrad)
    ∇u = Tensor{2,dim,T,ngrad}(x)
    return symmetric(∇u)
end


"""
This is a helper to access the correct value in Tensors.jl entities, because the gradient index is the outermost one.
"""
@inline _tensorsjl_gradient_accessor(v::Tensors.Vec{dim}, field_dim_idx::Int, spatial_dim_idx::Int) where {dim} = v[spatial_dim_idx]
@inline _tensorsjl_gradient_accessor(m::Tensors.Tensor{2,dim}, field_dim_idx::Int, spatial_dim_idx::Int) where {dim} = m[field_dim_idx, spatial_dim_idx]

"""
    interpolate_gradient_field(dh::DofHandler, u::AbstractVector, field_name::Symbol)

Compute the piecewise discontinuous gradient field for `field_name`. Returns the flux dof handler and the corresponding flux dof values.
"""
function interpolate_gradient_field(dh::Ferrite.DofHandler{spatial_dim}, u::AbstractVector, field_name::Symbol) where {spatial_dim}
    # Get some helpers
    field_idx = Ferrite.find_field(dh, field_name)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)

    # Create dof handler for gradient field
    dh_gradient = Ferrite.DofHandler(getgrid(dh))
    ip_gradient = get_gradient_interpolation(ip)
    field_dim = Ferrite.getfielddim(dh,field_name)
    push!(dh_gradient, :gradient, field_dim*spatial_dim, ip_gradient) # field dim × spatial dim components
    Ferrite.close!(dh_gradient)

    num_base_funs = Ferrite.getnbasefunctions(ip_gradient)

    # FIXME this does not work for mixed grids
    ip_geom = Ferrite.default_interpolation(typeof(Ferrite.getcells(getgrid(dh), 1)))
    ref_coords_gradient = Ferrite.reference_coordinates(ip_gradient)
    qr_gradient = Ferrite.QuadratureRule{spatial_dim, refshape(Ferrite.getcells(getgrid(dh), 1)), Float64}(ones(length(ref_coords_gradient)), ref_coords_gradient)
    cv = (field_dim == 1) ? Ferrite.CellScalarValues(qr_gradient, ip, ip_geom) : Ferrite.CellVectorValues(qr_gradient, ip, ip_geom)

    # Buffer for the dofs
    cell_dofs = zeros(Int, Ferrite.ndofs_per_cell(dh))
    cell_dofs_gradient = zeros(Int, Ferrite.ndofs_per_cell(dh_gradient))

    # Allocate storage for the fluxes to store
    u_gradient = zeros(Ferrite.ndofs(dh_gradient))
    # In general uᵉ_gradient is an order 3 tensor [field_dim, spatial_dim, num_base_funs]
    uᵉ_gradient = zeros(length(cell_dofs_gradient))
    uᵉ_gradient_view = reshape(uᵉ_gradient, (spatial_dim, field_dim, num_base_funs))
    uᵉ = zeros(field_dim*Ferrite.getnbasefunctions(ip))

    for (cell_num, cell) in enumerate(Ferrite.CellIterator(dh))
        # Get element dofs on parent field
        Ferrite.celldofs!(cell_dofs, dh, cell_num)
        uᵉ .= u[cell_dofs[Ferrite.dof_range(dh, field_name)]]

        # And initialize cellvalues for the cell to evaluate the gradient at the basis functions 
        # of the gradient field
        Ferrite.reinit!(cv, cell)

        # Now we simply loop over all basis functions of the gradient field and evaluate the gradient
        for i ∈ 1:num_base_funs
            uᵉgradi = Ferrite.function_gradient(cv, i, uᵉ)
            for ds in 1:spatial_dim
                for df in 1:field_dim
                    uᵉ_gradient_view[ds, df, i] = _tensorsjl_gradient_accessor(uᵉgradi, df, ds)
                end
            end
        end

        # We finally write back the result to the global dof vector of the gradient field
        Ferrite.celldofs!(cell_dofs_gradient, dh_gradient, cell_num)
        u_gradient[cell_dofs_gradient] .+= uᵉ_gradient
    end
    return dh_gradient, u_gradient
end

"""
Mapping from 2D triangle to 3D face of a tetrahedon.
"""
function transfer_quadrature_face_to_cell(point::AbstractVector, cell::Ferrite.AbstractCell{3,N,4}, face::Int) where {N}
    x,y = point
    face == 1 && return [ 1-x-y,  y,  0]
    face == 2 && return [ y,  0,  1-x-y]
    face == 3 && return [ x,  y,  1-x-y]
    face == 4 && return [ 0,  1-x-y,  y]
end

"""
Mapping from 2D quadrilateral to 3D face of a hexahedron.
"""
function transfer_quadrature_face_to_cell(point::AbstractVector, cell::Ferrite.AbstractCell{3,N,6}, face::Int) where {N}
    x,y = point
    face == 1 && return [ y,  x, -1]
    face == 2 && return [ x, -1,  y]
    face == 3 && return [ 1,  x,  y]
    face == 4 && return [-x,  1,  y]
    face == 5 && return [-1,  y,  x]
    face == 6 && return [ x,  y,  1]
end
