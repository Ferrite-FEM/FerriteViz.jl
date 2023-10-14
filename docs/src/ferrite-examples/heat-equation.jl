using Ferrite, SparseArrays

function assemble_heat_element!(Ke::Matrix, fe::Vector, cellvalues::CellValues, coords::Vector, rhs::Function)
    n_basefuncs = getnbasefunctions(cellvalues)

    fill!(Ke, 0)
    fill!(fe, 0)

    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        x = spatial_coordinate(cellvalues, q_point, coords)

        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)

            fe[i] += rhs(x) * δu * dΩ

            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke, fe
end

function assemble_steady_heat_global(cellvalues::CellValues, K::SparseMatrixCSC, dh::DofHandler, rhs::Function)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    f = zeros(ndofs(dh))

    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        assemble_heat_element!(Ke, fe, cellvalues, getcoordinates(cell), rhs)
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

function manufactured_heat_problem(element_type, ip, num_elements_per_dim)
    dim = Ferrite.getdim(ip)
    grid = generate_grid(element_type, ntuple(x->num_elements_per_dim, dim));
    ip_geo = Ferrite.default_interpolation(typeof(grid.cells[1]))
    qr = QuadratureRule{Ferrite.getrefshape(ip)}(2*Ferrite.getorder(ip))
    cellvalues = CellValues(qr, ip, ip_geo);

    ∂Ω = union(
        getfaceset(grid, "left"),
        getfaceset(grid, "right"),
        getfaceset(grid, "top"),
        getfaceset(grid, "bottom"),
    );

    if dim == 3
        ∂Ω = union(
            ∂Ω,
            getfaceset(grid, "front"),
            getfaceset(grid, "back")
        )
    end

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh);

    K = create_sparsity_pattern(dh)

    ch = ConstraintHandler(dh);
    dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
    add!(ch, dbc);

    close!(ch)
    update!(ch, 0.0);

    K, f = assemble_steady_heat_global(cellvalues, K, dh, x->(π/2)^2 * dim * prod(cos, x*π/2));
    apply!(K, f, ch)
    u = K \ f;

    return dh, u
end
