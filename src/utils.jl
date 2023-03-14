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

struct MakiePlotter{dim,DH<:Ferrite.AbstractDofHandler,T1,TOP<:Union{Nothing,Ferrite.AbstractTopology},T2,M,TRI} <: AbstractPlotter
    dh::DH
    u::Makie.Observable{Vector{T1}} # original solution on the original mesh (i.e. dh.mesh)
    topology::TOP
    visible::Vector{Bool} #TODO change from per cell to per triangle
    gridnodes::Vector{GeometryBasics.Point{dim,T2}} # coordinates of grid nodes in matrix form
    physical_coords::Vector{GeometryBasics.Point{dim,T2}} # coordinates in physical space of a vertex
    all_triangles::Vector{TRI}
    vis_triangles::ShaderAbstractions.Buffer{TRI,Vector{TRI}}
    triangle_cell_map::Vector{Int}
    reference_coords::Matrix{T1} # coordinates on the associated reference cell for the corresponding triangle vertex
    mesh::M
end

"""
    MakiePlotter(dh::Ferrite.AbstractDofHandler, u::Vector)
    MakiePlotter(dh::Ferrite.AbstractDofHandler, u::Vector, topology::TOP) where {TOP<:Ferrite.AbstractTopology}

Builds a static triangulation of the underlying `grid` in `dh.grid` for rendering via Makie.
The triangulation acts as a "L2" triangulation, i.e. each triangle node is doubled.
For large 3D grids, prefer to use the second constructor if you have already a `topology`.
Otherwise, it will be rebuilt which is time consuming.
"""
function MakiePlotter(dh::Ferrite.AbstractDofHandler, u::Vector, topology::TOP) where {TOP<:Union{Nothing,Ferrite.AbstractTopology}}
    cells = Ferrite.getcells(dh.grid)
    dim = Ferrite.getdim(dh.grid)
    visible = zeros(Bool,length(cells))
    if dim > 2
        boundaryfaces = findall(isempty,topology.face_neighbor)
        boundaryelements = Ferrite.getindex.(boundaryfaces,1)
    else
        boundaryelements = collect(1:Ferrite.getncells(dh.grid))
    end
    visible[boundaryelements] .= true

    # We do not take any assumptions on the mesh, so first we have to loopkup
    num_triangles = 0
    for cell in cells
        num_triangles += ntriangles(cell)
    end

    # Preallocate the matrices carrying the triangulation information
    triangles = Matrix{Int}(undef, num_triangles, 3)
    triangle_cell_map = Vector{Int}(undef, num_triangles)
    physical_coords = Vector{GeometryBasics.Point{dim,Float32}}(undef, num_triangles*3)
    gridnodes = [GeometryBasics.Point{dim,Float32}(node.data) for node in Ferrite.getcoordinates.(Ferrite.getnodes(dh.grid))]
    reference_coords = Matrix{Float64}(undef, num_triangles*3,dim)

    # Decompose does the heavy lifting for us
    coord_offset = 1
    triangle_offset = 1
    for (cell_id,cell) ∈ enumerate(cells)
        triangle_offset_begin = triangle_offset
        (coord_offset, triangle_offset) = decompose!(coord_offset, physical_coords, reference_coords, triangle_offset, triangles, dh.grid, cell)
        triangle_cell_map[triangle_offset_begin:(triangle_offset-1)] .= cell_id
    end
    all_triangles = Makie.to_triangles(triangles)
    vis_triangles = copy(all_triangles)
    n_visible = sum(visible[triangle_cell_map])
    n_notvisible = length(all_triangles) - n_visible
    vis_triangles[ .! visible[triangle_cell_map]] .= (GeometryBasics.GLTriangleFace(1,1,1) for i in 1:n_notvisible)
    vis_triangles = ShaderAbstractions.Buffer(Makie.Observable(vis_triangles))
    mesh = GeometryBasics.Mesh(physical_coords,vis_triangles)
    return MakiePlotter{dim,typeof(dh),eltype(u),typeof(topology),Float32,typeof(mesh),eltype(vis_triangles)}(dh,Observable(u),topology,visible,gridnodes,physical_coords,all_triangles,vis_triangles,triangle_cell_map,reference_coords,mesh)
end
MakiePlotter(dh,u) = MakiePlotter(dh,u,Ferrite.getdim(dh.grid) > 2 ? Ferrite.ExclusiveTopology(dh.grid.cells) : nothing)

# triangle_to_cell -> visible -> triangle access
visible(plotter::MakiePlotter{3}) = @views plotter.triangles[plotter.visible[plotter.triangle_cell_map],:]
visible(plotter::MakiePlotter{2}) = plotter.triangles

"""
Clip plane described by the normal and its distance to the coordinate origin.
"""
struct ClipPlane{T}
    normal::Tensors.Vec{3,T}
    distance::T
end

"""
Binary decision function to clip a cell with a plane for the crinkle clip.
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
    crinkle_clip!(plotter::MakiePlotter{3}, decision_fun)
Crinkle clip updates the visibility of the triangles, based on an
implicit description of the clipping surface. Here `decision_fun` takes the grid and
a cell index as input and returns whether the cell is visible or not.
Note that chained calls to `crinkle_clip!` won't work.
"""
function crinkle_clip!(plotter::MakiePlotter{3,DH,T}, decision_fun::DF) where {DH,T,DF}
    dh = plotter.dh
    u = plotter.u
    grid = dh.grid

    # We iterate over all triangles and check if the corresponding cell is visible.
    for (cell_id, cell) ∈ enumerate(Ferrite.getcells(plotter.dh.grid))
        dfun_visible = decision_fun(grid, cell_id)
        if dfun_visible
            cell_neighbors = Ferrite.getneighborhood(plotter.topology, grid, Ferrite.CellIndex(cell_id))
            plotter.visible[cell_id] = !all(decision_fun.((grid,),cell_neighbors)) || plotter.visible[cell_id]
        else
            plotter.visible[cell_id] = false
        end
    end
end

"""
    crinkle_clip(plotter::MakiePlotter{3}, decision_fun) -> MakiePlotter
Crinkle clip generates a new plotter with updated visibility of the triangles.
Non-mutating version of `crinkle_clip!`.
Note that chained calls to `crinkle_clip` won't work.
"""
function crinkle_clip(plotter::MakiePlotter{3,DH,T}, decision_fun) where {DH,T}
    plotter_clipped = MakiePlotter{3,DH,T,typeof(plotter.topology)}(plotter.dh,plotter.u,plotter.topology,copy(plotter.visible),plotter.gridnodes,plotter.physical_coords,plotter.triangles,plotter.triangle_cell_map,plotter.reference_coords);
    crinkle_clip!(plotter_clipped,decision_fun)
    return plotter_clipped
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
        coord_matrix[coord_offset] = GeometryBasics.Point(Ferrite.getcoordinates(v)...)
        ref_coord_matrix[coord_offset,1:2] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i]
        triangle_matrix[triangle_offset,i] = coord_offset
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
function decompose!(coord_offset, coord_matrix::Vector{Point{space_dim,T}}, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Union{Ferrite.AbstractCell{2,N,4}, Ferrite.AbstractCell{3,4,1}}) where {N,space_dim,T}
    # A bit more complicated. The default diagonal decomposition into 2 triangles is missing a solution mode.
    # To resolve this we make a more expensive decomposition in 4 triangles which correctly displays the solution mode (in linear problems).
    coord_offset_initial = coord_offset
    vts = vertices(grid, cell)

    # Compute coordinate of vertex 5
    center = zeros(space_dim)
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
        coord_matrix[coord_offset] = GeometryBasics.Point(Ferrite.getcoordinates(v1)...)
        ref_coord_matrix[coord_offset, 1:2] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i]
        coord_offset+=1

        # next vertex in chain
        coord_matrix[coord_offset] = GeometryBasics.Point(Ferrite.getcoordinates(v2)...)
        ref_coord_matrix[coord_offset, 1:2] = Ferrite.reference_coordinates(Ferrite.default_interpolation(typeof(cell)))[i2]
        coord_offset+=1

        # center vertex (5)
        coord_matrix[coord_offset] = GeometryBasics.Point(center...)
        ref_coord_matrix[coord_offset, 1:2] .= ntuple(x->0.0, 2)
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

"""
    postprocess(node_values::Vector{T}) -> T
Takes the nodal dof vector and maps it either to the scalar or to the
euclidean norm (in the vectorial case)
"""
function postprocess(node_values)
    dim = length(node_values)
    if dim == 1
        return node_values[1] #scalar values vectors with length 1
    else
        return sqrt(sum(node_values.^2))
    end
end

"""
    transfer_solution(plotter::MakiePlotter{dim,DH,T}, u::Vector; field_idx::Int=1, process::Function=FerriteViz.postprocess) where {dim,DH<:Ferrite.AbstractDofHandler,T}
Transfer the solution of a plotter to the tessellated mesh in `dim`.

@TODO Refactor. This is peak inefficiency.
"""
function transfer_solution(plotter::MakiePlotter{dim,DH,T}, u::Vector; field_idx::Int=1, process::FUN=FerriteViz.postprocess) where {dim,DH<:Ferrite.AbstractDofHandler,T,FUN}
    # select objects from plotter
    dh = plotter.dh
    grid = dh.grid

    # field related variables
    field_name = Ferrite.getfieldnames(dh)[field_idx]
    field_dim = Ferrite.getfielddim(dh, field_idx)

    # interpolation extraction
    ip_field = dh.field_interpolations[field_idx]
    cell_geo_ref = Ferrite.getcells(grid, 1)
    ip_geo = Ferrite.default_interpolation(typeof(cell_geo_ref))
    val_buffer = zeros(T,field_dim)
    val = process(val_buffer)

    return _transfer_solution(ip_geo,ip_field,val_buffer,val,field_name,field_dim,plotter,u,process) #function barrier for ip_field and thus pointvalues
end

function _transfer_solution(ip_geo,ip_field,val_buffer,val,field_name,field_dim,plotter::MakiePlotter{dim,DH,T}, u::Vector, process::FUN) where {dim,DH<:Ferrite.AbstractDofHandler,T,FUN}
    n_vertices_per_tri = 3 # we have 3 vertices per triangle...
    dh = plotter.dh
    ref_coords = plotter.reference_coords
    grid = dh.grid
    # actual data
    local_dof_range = Ferrite.dof_range(dh, field_name)

    pv = Ferrite.PointScalarValues(ip_field, ip_geo)

    current_vertex_index = 1

    Ferrite.reinit!(pv, Ferrite.getcoordinates(grid,1), Tensors.Vec{dim}(ref_coords[current_vertex_index,:]))
    n_basefuncs = Ferrite.getnbasefunctions(pv)
    _processreturn = length(process(val_buffer))

    data = fill(0.0, num_vertices(plotter),_processreturn)
    _local_coords = Ferrite.getcoordinates(grid,1)
    _local_celldofs = Ferrite.celldofs(dh,1)
    _celldofs_field = reshape(_local_celldofs[local_dof_range], (field_dim, n_basefuncs))
    _local_ref_coords = Tensors.Vec{dim}(ref_coords[1,:])

    for (isvisible,(cell_idx,cell_geo)) in zip(plotter.visible,enumerate(Ferrite.getcells(dh.grid)))
        # This should make the loop work for mixed grids
        if !isvisible
            current_vertex_index += ntriangles(cell_geo)*n_vertices_per_tri
            continue
        end
        #if typeof(cell_geo_ref) != typeof(cell_geo)
        #    ip_geo = Ferrite.default_interpolation(typeof(cell_geo))
        #    pv = getpointvalues(Val(field_dim),ip_field,ip_geo)
        #    cell_geo_ref = cell_geo
        #end
        Ferrite.getcoordinates!(_local_coords,grid,cell_idx)
        Ferrite.celldofs!(_local_celldofs,dh,cell_idx)
        _celldofs_field = reshape(@view(_local_celldofs[local_dof_range]), (field_dim, n_basefuncs))
        ncvertices = ntriangles(cell_geo)*n_vertices_per_tri
        # TODO replace this with a triangle-to-cell map.
        for i in 1:ncvertices
            _local_ref_coords = Tensors.Vec{dim}(@view(ref_coords[current_vertex_index,:]))
            Ferrite.reinit!(pv, _local_coords, _local_ref_coords)
            for d in 1:field_dim
                val_buffer[d] = Ferrite.function_value(pv, 1, @view(u[_celldofs_field[d,:]]))
            end
            val = process(val_buffer)
            for d in 1:_processreturn
                data[current_vertex_index,d] = val[d]
            end
            current_vertex_index += 1
        end
    end
    return data
end

function transfer_scalar_celldata(plotter::MakiePlotter{dim,DH,T}, u::Vector; process::FUN=FerriteViz.postprocess) where {dim,DH<:Ferrite.AbstractDofHandler,T,FUN<:Function}
    n_vertices_per_tri = 3 # we have 3 vertices per triangle...

    # select objects from plotter
    dh = plotter.dh
    grid = dh.grid

    current_vertex_index = 1
    data = fill(0.0, num_vertices(plotter))
    for (isvisible,(cell_idx,cell_geo)) in zip(plotter.visible,enumerate(Ferrite.getcells(dh.grid)))
        if !isvisible
            current_vertex_index += ntriangles(cell_geo)*n_vertices_per_tri
            continue
        end
        ncvertices = ntriangles(cell_geo)*n_vertices_per_tri
        for i in 1:ncvertices
            data[current_vertex_index] = process(u[cell_idx])
            current_vertex_index += 1
        end
    end

    return data::Vector{T}
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

# maps the dof vector in nodal order, only needed for wireframe nodal deformation (since we display the original nodes)
function dof_to_node(dh::Ferrite.AbstractDofHandler, u::Vector{T}; field::Int=1) where T
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
    return data::Matrix{T}
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
