# These functions generate the corresponding first order cells of an interpolation.
# Triangle
for_nodes(::Union{Lagrange{RefTriangle,1},DiscontinuousLagrange{RefTriangle,1},Triangle}) = (
    (3,1,2),
)
# Quadratic Triangle
for_nodes(::Union{Lagrange{RefTriangle,2},DiscontinuousLagrange{RefTriangle,2},QuadraticTriangle}) = (
    (6,1,4),
    (5,6,4),
    (3,6,5),
    (5,4,2),
)
# Cubic Triangle
for_nodes(::Union{Lagrange{RefTriangle,3},DiscontinuousLagrange{RefTriangle,3}}) = (
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
for_nodes(::Union{Lagrange{RefTriangle,4},DiscontinuousLagrange{RefTriangle,4}}) = (
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
for_nodes(::Union{Lagrange{RefTriangle,5},DiscontinuousLagrange{RefTriangle,5}}) = (
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
for_nodes(::Union{Lagrange{RefTetrahedron,1},DiscontinuousLagrange{RefTetrahedron,1},Tetrahedron}) = (
    (1,2,3,4),
)
# Quadratic Tetrahedron
for_nodes(::Union{Lagrange{RefTetrahedron,2},DiscontinuousLagrange{RefTetrahedron,2},QuadraticTetrahedron}) = (
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
for_nodes(::Union{Lagrange{RefQuadrilateral,1},DiscontinuousLagrange{RefQuadrilateral,1},Quadrilateral}) = (
    (1,2,3,4),
)
# Quadratic Quadrilateral
for_nodes(::Union{Lagrange{RefQuadrilateral,2},DiscontinuousLagrange{RefQuadrilateral,2},QuadraticQuadrilateral}) = (
    (1,5,9,8),
    (5,2,6,9),
    (9,6,3,7),
    (8,9,7,4),
)
# Hexahedron
for_nodes(::Union{Lagrange{RefHexahedron,1},DiscontinuousLagrange{RefHexahedron,1},Hexahedron}) = (
    (1,2,3,4,5,6,7,8),
)
# Quadratic Hexahedron
for_nodes(::Union{Lagrange{RefHexahedron,2},DiscontinuousLagrange{RefHexahedron,2},QuadraticHexahedron}) = (
    (1,9,21,12,17,22,27,25),
    (17,22,27,25,5,13,26,16),
    (9,2,10,21,22,18,23,27),
    (22,18,23,27,13,6,14,26),
    (12,21,11,4,25,27,24,20),
    (25,27,24,20,16,26,15,8),
    (21,10,3,11,27,23,19,24),
    (27,23,19,24,26,14,7,15)
)

"""
Get the interpolation of the first order refinement. 
"""
for_interpolation(ip::Lagrange{shape,order}) where {shape,order} = Lagrange{shape,1}()

for_base_geometry_type(ip::Lagrange{RefQuadrilateral,order}) where {order} = Quadrilateral
for_base_geometry_type(ip::Lagrange{RefHexahedron,order}) where {order} = Hexahedron
for_base_geometry_type(ip::Lagrange{RefTriangle,order}) where {order} = Triangle
for_base_geometry_type(ip::Lagrange{RefTetrahedron,order}) where {order} = Tetrahedron

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
    @assert length(dh.subdofhandlers) == 1 "Subdomains not supported yet"
    # TODO Dofs for fields are not continuous. Think harder.
    @assert Ferrite.nfields(dh) == 1 "Multiple fields not supported yet"
    sdh = dh.subdofhandlers[1]
    field_idx=1
    grid = dh.grid

    # Some helpers
    ip = Ferrite.getfieldinterpolation(sdh, field_idx)
    field_dim = Ferrite.n_components(sdh, field_idx)
    spatial_dim = Ferrite.getspatialdim(grid)

    # # Get dof range, the hard way
    # dof_min = dh.ndofs.x
    # dof_max = 0
    ncells = 0
    for cell ∈ Ferrite.CellIterator(sdh)
        # celldofs = Ferrite.celldofs(cell)
        # dof_max = max(dof_max, maximum(celldofs))
        # dof_min = min(dof_min, minimum(celldofs))
        ncells += length(for_nodes(ip))
    end

    # Preallocate
    nodes = Vector{typeof(grid.nodes[1])}(undef, Ferrite.ndofs(dh)) #(dof_max-dof_min+1)÷field_dim)
    cells = Vector{Ferrite.getcelltype(grid)}(undef, ncells)

    ref_coords = Ferrite.reference_coordinates(ip)
    # Starting here we assume a single type of cell being present
    # TODO improve this.
    ip_geo = Ferrite.geometric_interpolation(typeof(grid.cells[1]))
    nodes_per_cell = length(ref_coords)
    qr = Ferrite.QuadratureRule{Ferrite.getrefshape(ip)}(zeros(nodes_per_cell), ref_coords)
    cv = Ferrite.CellValues(qr, ip, ip_geo)
    cellidx = 1
    for cell ∈ Ferrite.CellIterator(sdh)
        Ferrite.reinit!(cv, cell)
        coords = Ferrite.getcoordinates(cell)
        dofs_f = Ferrite.celldofs(cell)[Ferrite.dof_range(sdh, field_idx)]

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
    grid_new = Grid(cells, nodes)
    dh_new = DofHandler(grid_new)
    add!(dh_new, getfieldname(sdh, field_idx), Ferrite.n_components(sdh, field_idx), for_interpolation(ip))
    close!(dh_new);

    # Transfer solution the dumb way.
    # TODO this can be optimized.
    u_new = zeros(Ferrite.ndofs(dh_new))
    dh_dof_range = Ferrite.celldofs(dh_new, 1)
    for cell_idx ∈ 1:length(dh_new.grid.cells)
        Ferrite.celldofs!(dh_dof_range, dh_new, cell_idx)
        dofs = dh_dof_range[Ferrite.dof_range(dh_new.subdofhandlers[1], 1)]
        u_new[dofs] .= u[[dh_new.grid.cells[cell_idx].nodes...]]
    end

    return dh_new, u_new
end
