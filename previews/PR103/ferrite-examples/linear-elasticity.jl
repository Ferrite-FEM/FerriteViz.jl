using Ferrite, SparseArrays

function assemble_linear_elastic_element!(Ke::Matrix, fe::Vector, cellvalues::CellValues, coords::Vector, rhs::Function)
    n_basefuncs = getnbasefunctions(cellvalues)

    E = 200e3 # Young's modulus [MPa]
    ν = 0.3 # Poisson's ratio [-]

    G = E/2(1+ν)
    K = E/3(1-2ν)

    fill!(Ke, 0)
    fill!(fe, 0)

    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            Nᵢ = shape_value(cellvalues, q_point, i)
            x = spatial_coordinate(cellvalues, q_point, coords)
            fe[i] += rhs(x) ⋅ Nᵢ * dΩ
            ∇ˢʸᵐNᵢ = shape_symmetric_gradient(cellvalues, q_point, i)
            for j in 1:n_basefuncs
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                σ = 2G * dev(∇ˢʸᵐNⱼ) + K * tr(∇ˢʸᵐNⱼ) * one(∇ˢʸᵐNⱼ)
                Ke[i, j] += σ ⊡ ∇ˢʸᵐNᵢ * dΩ
            end
        end
    end

    return Ke, fe
end

function assemble_steady_linear_elastic_global(cellvalues::CellValues, K::SparseMatrixCSC, dh::DofHandler, rhs::Function)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)

    f = zeros(ndofs(dh))

    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        assemble_linear_elastic_element!(Ke, fe, cellvalues, getcoordinates(cell), rhs)
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

function manufactured_linear_elastic_problem(element_type, ip, num_elements_per_dim, component)
    dim = Ferrite.getrefdim(ip)
    grid = generate_grid(element_type, ntuple(x->num_elements_per_dim, dim));
    ip_geo = Ferrite.geometric_interpolation(typeof(grid.cells[1]))
    qr = QuadratureRule{Ferrite.getrefshape(ip)}(2*Ferrite.getorder(ip))
    cellvalues = CellValues(qr, ip, ip_geo);

    ∂Ω = union(
        getfacetset(grid, "left"),
        getfacetset(grid, "right"),
        getfacetset(grid, "top"),
        getfacetset(grid, "bottom"),
    );

    if dim == 3
        ∂Ω = union(
            ∂Ω,
            getfacetset(grid, "front"),
            getfacetset(grid, "back")
        )
    end

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh);

    K = allocate_matrix(dh)

    ch = ConstraintHandler(dh);
    dbc = Dirichlet(:u, ∂Ω, (x, t) -> [0.0 for i ∈ 1:dim])
    add!(ch, dbc);

    close!(ch)
    update!(ch, 0.0);

    function analytical_rhs(x_eval)
        E = 200e3 # Young's modulus [MPa]
        ν = 0.3 # Poisson's ratio [-]

        G = E/2(1+ν)
        K = E/3(1-2ν)

        u(x)  = Vec{2}(i->i==component ? (π/2)^2 * 2 * prod(cos, x*π/2) : 0.0)
        ∇u(x) = Tensors.gradient(u, x)
        ε(x)  = (∇u(x)  + transpose(∇u(x)))/2
        σ(x)  = 2G * dev(ε(x)) + K * tr(ε(x)) * one(ε(x))

        return Tensors.divergence(σ, x_eval)
    end

    K, f = assemble_steady_linear_elastic_global(cellvalues, K, dh, analytical_rhs);
    apply!(K, f, ch)
    u = K \ f;

    return dh, u
end

Tensors.tr(t::Tensor{3,dim}) where dim = Vec{dim}(i->sum(t[:,i,:]))
