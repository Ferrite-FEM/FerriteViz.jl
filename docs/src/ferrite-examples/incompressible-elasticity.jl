using Ferrite
using BlockArrays, SparseArrays, LinearAlgebra

function create_cook_grid(nx, ny)
    corners = [Tensors.Vec{2}((0.0,   0.0)),
               Tensors.Vec{2}((48.0, 44.0)),
               Tensors.Vec{2}((48.0, 60.0)),
               Tensors.Vec{2}((0.0,  44.0))]
    grid = generate_grid(Quadrilateral, (nx, ny), corners);
    # facesets for boundary conditions
    addfacetset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0);
    addfacetset!(grid, "traction", x -> norm(x[1]) ≈ 48.0);
    return grid
end;

function create_values(interpolation_u, interpolation_p)
    # quadrature rules
    qr      = QuadratureRule{RefQuadrilateral}(3)
    face_qr = FacetQuadratureRule{RefQuadrilateral}(3)

    # geometric interpolation
    interpolation_geom = Lagrange{RefQuadrilateral,1}()

    # cell and facevalues for u
    cellvalues_u = CellValues(qr, interpolation_u, interpolation_geom)
    facevalues_u = FacetValues(face_qr, interpolation_u, interpolation_geom)

    # cellvalues for p
    cellvalues_p = CellValues(qr, interpolation_p, interpolation_geom)

    return cellvalues_u, cellvalues_p, facevalues_u
end;

function create_dofhandler(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu) # displacement
    add!(dh, :p, ipp) # pressure
    close!(dh)
    return dh
end;

function create_bc(dh)
    dbc = ConstraintHandler(dh)
    add!(dbc, Dirichlet(:u, getfacetset(dh.grid, "clamped"), (x,t) -> zero(Tensors.Vec{2}), [1,2]))
    close!(dbc)
    t = 0.0
    update!(dbc, t)
    return dbc
end;

struct LinearElasticity{T}
    G::T
    K::T
end

function doassemble(cellvalues_u::CellValues, cellvalues_p::CellValues,
                    facevalues_u::FacetValues, K::SparseMatrixCSC, grid::Grid{sdim},
                    dh::DofHandler, mp::LinearElasticity) where {sdim}

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    fe = BlockedArray(zeros(nu + np), [nu, np]) # local force vector
    ke = BlockedArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix

    # traction vector
    t = Tensors.Vec{2}((0.0, 1/16))
    # cache ɛdev outside the element routine to avoid some unnecessary allocations
    ɛdev = [zero(SymmetricTensor{2, sdim}) for i in 1:getnbasefunctions(cellvalues_u)]

    for cell in CellIterator(dh)
        fill!(ke, 0)
        fill!(fe, 0)
        assemble_up!(ke, fe, cell, cellvalues_u, cellvalues_p, facevalues_u, grid, mp, ɛdev, t)
        assemble!(assembler, celldofs(cell), ke, fe)
    end

    return K, f
end;

function assemble_up!(Ke, fe, cell, cellvalues_u, cellvalues_p, facevalues_u, grid, mp, ɛdev, t)

    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    u▄, p▄ = 1, 2
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_p, cell)

    # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    @inbounds for q_point in 1:getnquadpoints(cellvalues_u)
        for i in 1:n_basefuncs_u
            ɛdev[i] = dev(symmetric(shape_gradient(cellvalues_u, q_point, i)))
        end
        dΩ = getdetJdV(cellvalues_u, q_point)
        for i in 1:n_basefuncs_u
            divδu = shape_divergence(cellvalues_u, q_point, i)
            δu = shape_value(cellvalues_u, q_point, i)
            for j in 1:i
                Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * mp.G * ɛdev[i] ⊡ ɛdev[j] * dΩ
            end
        end

        for i in 1:n_basefuncs_p
            δp = shape_value(cellvalues_p, q_point, i)
            for j in 1:n_basefuncs_u
                divδu = shape_divergence(cellvalues_u, q_point, j)
                Ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
            end
            for j in 1:i
                p = shape_value(cellvalues_p, q_point, j)
                Ke[BlockIndex((p▄, p▄), (i, j))] += - 1/mp.K * δp * p * dΩ
            end

        end
    end

    symmetrize_lower!(Ke)

    # We integrate the Neumann boundary using the facevalues.
    # We loop over all the faces in the cell, then check if the face
    # is in our `"traction"` faceset.
    @inbounds for face in 1:nfacets(cell)
        if (cellid(cell), face) ∈ getfacetset(grid, "traction")
            reinit!(facevalues_u, cell, face)
            for q_point in 1:getnquadpoints(facevalues_u)
                dΓ = getdetJdV(facevalues_u, q_point)
                for i in 1:n_basefuncs_u
                    δu = shape_value(facevalues_u, q_point, i)
                    fe[i] += (δu ⋅ t) * dΓ
                end
            end
        end
    end
end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

function solve(interpolation_u, interpolation_p, mp)
    # grid, dofhandler, boundary condition
    n = 50
    grid = create_cook_grid(n, n)
    dh = create_dofhandler(grid, interpolation_u, interpolation_p)
    dbc = create_bc(dh)

    # cellvalues
    cellvalues_u, cellvalues_p, facevalues_u = create_values(interpolation_u, interpolation_p)

    # assembly and solve
    K = allocate_matrix(dh);
    K, f = doassemble(cellvalues_u, cellvalues_p, facevalues_u, K, grid, dh, mp);
    apply!(K, f, dbc)
    u = K \ f;

    # export
    # filename = "cook_" * (isa(interpolation_u, Lagrange{RefQuadrilateral,1}) ? "linear" : "quadratic") *
    #                      "_linear"
    # vtk_grid(filename, dh) do vtkfile
    #     vtk_point_data(vtkfile, dh, u)
    # end
    return u,dh
end

linear    = Lagrange{RefQuadrilateral,1}()
quadratic = Lagrange{RefQuadrilateral,2}()

ν = 0.4999999
Emod = 1.
Gmod = Emod / 2(1 + ν)
Kmod = Emod * ν / ((1+ν) * (1-2ν))
mp = LinearElasticity(Gmod, Kmod)
u_linear,dh_linear = solve(linear^2, linear, mp);
u_quadratic,dh_quadratic = solve(quadratic^2, linear, mp);
