# These functions generate the corresponding first order cells of an interpolation.
# Triangle
for_nodes(::Union{Ferrite.Lagrange{2,Ferrite.RefTetrahedron,1},Ferrite.DiscontinuousLagrange{2,Ferrite.RefTetrahedron,1},Ferrite.Triangle}) = (
    (3,1,2),
)
# Quadratic Triangle
for_nodes(::Union{Ferrite.Lagrange{2,Ferrite.RefTetrahedron,2},Ferrite.DiscontinuousLagrange{2,Ferrite.RefTetrahedron,2},Ferrite.QuadraticTriangle}) = (
    (6,1,4),
    (5,6,4),
    (3,6,5),
    (5,4,2),
)
# Cubic Triangle
for_nodes(::Union{Ferrite.Lagrange{2,Ferrite.RefTetrahedron,3},Ferrite.DiscontinuousLagrange{2,Ferrite.RefTetrahedron,3},Ferrite.Cell{2,10,3}}) = (
    (3,8,7),
    (7,8,10),
    (8,9,10),
    (10,9,4),
    (9,1,4),
    (7,10,6),
    (6,10,5),
    (6,5,2),
    (10,4,5),
)
# Biquadratic Triangle
for_nodes(::Union{Ferrite.Lagrange{2,Ferrite.RefTetrahedron,4},Ferrite.DiscontinuousLagrange{2,Ferrite.RefTetrahedron,4},Ferrite.Cell{2,15,3}}) = (
    (3,10,9),
    (13,9,10),
    (10,11,13),
    (14,13,11),
    (11,12,14),
    (4,14,12),
    (12,1,4),
    (9,13,8),
    (15,8,13),
    (13,14,15),
    (5,15,14),
    (14,4,5),
    (8,15,7),
    (6,7,15),
    (15,5,6),
    (7,6,2),
)
# Quintic Triangle
for_nodes(::Union{Ferrite.Lagrange{2,Ferrite.RefTetrahedron,5},Ferrite.DiscontinuousLagrange{2,Ferrite.RefTetrahedron,5},Ferrite.Cell{2,20,3}}) = (
    (3,12,11),
    (16,11,12),
    (12,13,16),
    (17,16,13),
    (13,14,17),
    (18,17,14),
    (14,15,18),
    (4,18,15),
    (15,1,4),
    (11,16,10),
    (19,10,16),
    (16,17,19),
    (20,19,17),
    (17,18,20),
    (5,20,18),
    (18,4,5),
    (10,19,9),
    (21,9,19),
    (19,20,21),
    (6,21,20),
    (20,5,6),
    (9,21,8),
    (7,8,21),
    (21,6,7),
    (8,7,2),
)
# Tetrahedron
for_nodes(::Union{Ferrite.Lagrange{3,Ferrite.RefTetrahedron,1},Ferrite.DiscontinuousLagrange{3,Ferrite.RefTetrahedron,1},Ferrite.Tetrahedron}) = (
    (3,1,2,4),
)
# Quadratic Tetrahedron
for_nodes(::Union{Ferrite.Lagrange{3,Ferrite.RefTetrahedron,2},Ferrite.DiscontinuousLagrange{3,Ferrite.RefTetrahedron,2},Ferrite.QuadraticTetrahedron}) = (
    (5,2,6,9),
    (7,6,3,10),
    (8,9,10,4),
    (8,5,6,9),
    (8,6,7,10),
    (5,8,1,6),
    (7,6,1,8),
    (9,10,8,6),
)
# Quadrilateral
for_nodes(::Union{Ferrite.Lagrange{2,Ferrite.RefCube,1},Ferrite.DiscontinuousLagrange{2,Ferrite.RefCube,1},Ferrite.Quadrilateral}) = (
    (1,2,3,4),
)
# Quadratic Quadrilateral
for_nodes(::Union{Ferrite.Lagrange{2,Ferrite.RefCube,2},Ferrite.DiscontinuousLagrange{2,Ferrite.RefCube,2},Ferrite.QuadraticQuadrilateral}) = (
    (1,5,9,8),
    (5,2,6,9),
    (9,6,3,7),
    (8,9,7,4),
)

"""
Get the interpolation of the first order refinement. 
"""
for_interpolation(ip::Ferrite.Lagrange{dim,shape,order}) where {dim,shape,order} = Ferrite.Lagrange{dim,shape,1}()

for_base_geometry_type(ip::Ferrite.Lagrange{2,Ferrite.RefCube,order}) where {order} = Ferrite.Quadrilateral
for_base_geometry_type(ip::Ferrite.Lagrange{3,Ferrite.RefCube,order}) where {order} = Ferrite.Hexahedron
for_base_geometry_type(ip::Ferrite.Lagrange{2,Ferrite.RefTetrahedron,order}) where {order} = Ferrite.Triangle
for_base_geometry_type(ip::Ferrite.Lagrange{3,Ferrite.RefTetrahedron,order}) where {order} = Ferrite.Tetrahedron

# TODO move into ferrite core
function Ferrite.field_offset(dh::Ferrite.DofHandler, field_name::Int)
    offset = 0
    for i in 1:field_name-1
        offset += Ferrite.getnbasefunctions(dh.field_interpolations[i])::Int * dh.field_dims[i]
    end
    return offset
end

# TODO move into ferrite core
function Ferrite.dof_range(dh::Ferrite.DofHandler, field_idx::Int)
    offset = Ferrite.field_offset(dh, field_idx)
    n_field_dofs = Ferrite.getnbasefunctions(dh.field_interpolations[field_idx])::Int * dh.field_dims[field_idx]
    return (offset+1):(offset+n_field_dofs)
end

# TODO move to ferrite core
getfieldname(dh, field_idx) = dh.field_names[field_idx]

"""
Create a first order discretization w.r.t. a field and transfer
the solution.
"""
function for_discretization(dh, u)
    field_idx=1

    # Some helpers
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    field_dim = Ferrite.getfielddim(dh, field_idx)
    spatial_dim = Ferrite.getdim(dh.grid)

    # TODO Dofs for fields are not continuous. Think harder.
    @assert Ferrite.nfields(dh) == 1 "Multiple fields not supported yet"
    # # Get dof range, the hard way
    # dof_min = dh.ndofs.x
    # dof_max = 0
    ncells = 0
    for cell ∈ Ferrite.CellIterator(dh)
        # celldofs = Ferrite.celldofs(cell)
        # dof_max = max(dof_max, maximum(celldofs))
        # dof_min = min(dof_min, minimum(celldofs))
        ncells += length(for_nodes(ip))
    end

    # Preallocate
    nodes = Vector{typeof(dh.grid.nodes[1])}(undef, Ferrite.ndofs(dh)) #(dof_max-dof_min+1)÷field_dim)
    cells = Vector{Ferrite.getcelltype(dh.grid)}(undef, ncells)

    ref_coords = Ferrite.reference_coordinates(ip)
    # Starting here we assume a single type of cell being present
    # TODO improve this.
    ip_geo = Ferrite.default_interpolation(typeof(dh.grid.cells[1]))
    nodes_per_cell = length(ref_coords)
    qr = Ferrite.QuadratureRule{spatial_dim, Ferrite.getrefshape(ip)}(zeros(nodes_per_cell), ref_coords)
    cv = Ferrite.CellScalarValues(qr, ip, ip_geo)
    cellidx = 1
    for cell ∈ Ferrite.CellIterator(dh)
        Ferrite.reinit!(cv, cell)
        coords = Ferrite.getcoordinates(cell)
        dofs_f = Ferrite.celldofs(cell)[Ferrite.dof_range(dh, field_idx)]

        # Extract coordinates
        for q ∈ 1:nodes_per_cell
            nodes[dofs_f[q]] = Ferrite.Node(Ferrite.spatial_coordinate(cv, q, coords))
        end
        # And the 
        for subcellnodes ∈ for_nodes(ip)
            # Splatting sorcery to extract the global node indices.
            cells[cellidx] = for_base_geometry_type(ip)(((dofs_f[[subcellnodes...]]...,)))
            cellidx += 1
        end
    end

    # Generate a new dof handler.
    grid_new = Ferrite.Grid(cells, nodes)
    dh_new = Ferrite.DofHandler(grid_new)
    Ferrite.push!(dh_new, getfieldname(dh, field_idx), Ferrite.getfielddim(dh, field_idx), for_interpolation(ip))
    Ferrite.close!(dh_new);

    # Transfer solution the dumb way.
    # TODO this can be optimized.
    u_new = zeros(Ferrite.ndofs(dh_new))
    for cell_idx ∈ 1:length(dh_new.grid.cells)
        dh_dof_range = dh_new.cell_dofs_offset[cell_idx]:(dh_new.cell_dofs_offset[cell_idx+1]-1)
        dofs = dh_new.cell_dofs[dh_dof_range][Ferrite.dof_range(dh_new, field_idx)]
        u_new[dofs] .= u[[dh_new.grid.cells[cell_idx].nodes...]]
    end

    return dh_new, u_new
end
