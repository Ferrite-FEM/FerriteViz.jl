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
face_cell(cell::Ferrite.Tetrahedron, i::Int) =  Ferrite.Cell{3,3,1}(Ferrite.faces(cell)[i])
face_cell(cell::Ferrite.Hexahedron, i::Int) = Ferrite.Quadrilateral3D(Ferrite.faces(cell)[i])

struct MakiePlotter{dim,DH<:Ferrite.AbstractDofHandler} <: AbstractPlotter
    dh::DH
    u::Vector # original solution on the original mesh (i.e. dh.mesh)

    cells_connectivity::Vector

    physical_coords::Matrix{AbstractFloat} # coordinates in physical space of a vertex
    triangles::Matrix{Int} # each row carries a triple with the indices into the coords matrix defining the triangle
    reference_coords::Matrix{AbstractFloat} # coordinates on the reference cell for the corresponding vertex
end

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
    physical_coords = Matrix{Float64}(undef, num_triangles*3, dim)
    reference_coords = Matrix{Float64}(undef, num_triangles*3, 2) # for now no 1D

    # Decompose does the heavy lifting for us
    coord_offset = 1
    triangle_offset = 1
    for cell in cells
        (coord_offset, triangle_offset) = decompose!(coord_offset, physical_coords, reference_coords, triangle_offset, triangles, dh.grid, cell)
    end
    return MakiePlotter{dim,typeof(dh)}(dh,u,[],physical_coords,triangles,reference_coords);
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
ntriangles(cell::Ferrite.AbstractCell{2,3,3}) where {N} = 1 # Tris in 2D
ntriangles(cell::Ferrite.AbstractCell{3,3,1}) where {N} = 1 # Tris in 3D
ntriangles(cell::Ferrite.AbstractCell{dim,N,4}) where {dim,N} = 4 # Quads in 2D and 3D
ntriangles(cell::Ferrite.AbstractCell{3,N,1}) where {dim,N} = 4 # Tets as a special case of a Quad, obviously :)
ntriangles(cell::Ferrite.AbstractCell{3,N,6}) where {dim,N} = 6*4 # Hex

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
        ref_coord_matrix[coord_offset, :] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i]
        triangle_matrix[triangle_offset, i] = coord_offset
        coord_offset+=1
    end
    triangle_offset+=1
    (coord_offset, triangle_offset)
end

"""
Decompose a tetrahedron into a coordinates and a triangle index list to disconnect it properly. Guarantees to preserve orderings and orientations.
"""
function decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Ferrite.AbstractCell{3,N,4}) where {N}
    # Is just 4 triangles :)
    for ti ∈ Ferrite.faces(cell)
        (coord_offset, triangle_offset) = decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, Ferrite.Cell{3,3,1}(ti))
    end
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
        ref_coord_matrix[coord_offset, :] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i]
        coord_offset+=1

        # next vertex in chain
        coord_matrix[coord_offset, :] = Ferrite.getcoordinates(v2)
        ref_coord_matrix[coord_offset, :] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i2]
        coord_offset+=1

        # center vertex (5)
        coord_matrix[coord_offset, :] = center
        ref_coord_matrix[coord_offset, :] = [ntuple(x->0.0, 2)...]
        coord_offset+=1

        # collect indices
        triangle_matrix[triangle_offset, :] = (0:2) .+ (coord_offset_initial+3*(i-1))
        triangle_offset+=1
    end

    (coord_offset, triangle_offset)
end

"""
Decompose a hexahedron into a coordinates and a triangle index list to disconnect it properly. Guarantees to preserve orderings and orientations.
"""
function decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Ferrite.AbstractCell{3,N,6}) where {N}
    # Just 6 quadrilaterals :)
    for ti ∈ Ferrite.faces(cell)
        (coord_offset, triangle_offset) = decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, Ferrite.Cell{3,4,1}(ti))
    end
    (coord_offset, triangle_offset)
end


refshape(cell::Ferrite.AbstractCell) = typeof(Ferrite.default_interpolation(typeof(cell))).parameters[2]

function postprocess(node_values)
    dim = length(node_values)
    if dim == 1
        return node_values
    else
        return sqrt(sum(node_values.^2))
    end
end


using Tensors
"""
Transfer the solution of a plotter to the new mesh in 2D.

@precondition: Assumes that the number of basis function for each dof is equal.

@TODO Refactor. This is peak inefficiency.
"""
function transfer_solution(plotter::MakiePlotter{2}; field_idx::Int=1, process::Function=FerriteVis.postprocess)
    n_vertices = 3 # we have 3 vertices per triangle...

    # select objects from plotter
    dh = plotter.dh
    u = plotter.u
    ref_coords = plotter.reference_coords
    grid = dh.grid

    # field related variables
    field_name = Ferrite.getfieldnames(dh)[field_idx]
    field_dim = Ferrite.getfielddim(dh, field_idx)
    ip = dh.field_interpolations[field_idx]

    # actual data
    data = fill(0.0, num_vertices(plotter), field_dim)
    local_dof_range = Ferrite.dof_range(dh, field_name)

    current_vertex_index = 1
    for cell in Ferrite.CellIterator(dh)
        cell_idx = cell.current_cellid.x

        cell_geo = dh.grid.cells[cell_idx]
        _celldofs_field = reshape(Ferrite.celldofs(cell)[local_dof_range], (field_dim, Ferrite.getnbasefunctions(ip)))

        # Loop over vertices
        for i in 1:(ntriangles(cell_geo)*n_vertices)
            ξ = Vec(ref_coords[current_vertex_index, :]...)
            for d in 1:field_dim
                for node_idx ∈ 1:Ferrite.getnbasefunctions(ip)
                    data[current_vertex_index, d] += Ferrite.value(ip, node_idx, ξ) ⋅ u[_celldofs_field[d, node_idx]]
                end
            end
            current_vertex_index += 1
        end
    end
    return mapslices(process, data, dims=[2])
end

"""
Transfer the solution of a plotter to the new mesh in 3D. We just evaluate the faces.

@TODO Refactor. This is peak inefficiency.
"""
function transfer_solution(plotter::MakiePlotter{3}; field_idx::Int=1, process::Function=FerriteVis.postprocess)
    # select objects from plotter
    dh = plotter.dh
    u = plotter.u
    ref_coords = plotter.reference_coords
    grid = dh.grid

    # field related variables
    field_name = Ferrite.getfieldnames(dh)[field_idx]
    field_dim = Ferrite.getfielddim(dh, field_idx)
    ip_cell = dh.field_interpolations[field_idx]
    faces_ = Ferrite.faces(ip_cell) # faces of the cell with local dofs

    # face related variables
    ip_face = Ferrite.getlowerdim(ip_cell)

    # actual data
    data = fill(0.0, num_vertices(plotter), field_dim)
    offset = Ferrite.field_offset(dh, field_name)

    current_vertex_index = 1
    for (cell_index, cell) in enumerate(Ferrite.CellIterator(dh))
        cell_geo = Ferrite.getcells(grid)[cell.current_cellid.x]
        celldofs_ = Ferrite.celldofs(cell)
        for (local_face_idx,_) in enumerate(Ferrite.faces(cell_geo))
            # extract face vertex dofs
            face_vertex_incides = faces_[local_face_idx][1:2]
            facedofs_ = celldofs_[[face_vertex_incides...].+ offset]

            face_geo = face_cell(cell_geo, local_face_idx)

            # Loop over vertices
            for i in 1:(ntriangles(face_geo)*3)
                ξ = Vec(ref_coords[current_vertex_index, :]...)
                for d in 1:field_dim
                    ndofs = length(facedofs_) ÷ field_dim

                    field_dof_indices = [field_dim*(dofi-1)+d for dofi in 1:ndofs]

                    actual_dofs = facedofs_[field_dof_indices]
                    for i in 1:ndofs
                        data[current_vertex_index, d] += Ferrite.value(ip_face, field_dim*(i-1)+d, ξ) ⋅ u[actual_dofs[i]]
                    end
                end
                current_vertex_index += 1
            end
        end
    end
    return mapslices(process, data, dims=[2])
end