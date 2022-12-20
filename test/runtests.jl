using FerriteViz
using Test

include("../docs/src/ferrite-examples/heat-equation.jl")

# TODO move this into Ferrite core

@testset "gradient fields (scalar)" begin
    # Check scalar problems
    for (size, geo, ip) ∈ [(21, Triangle, Lagrange{2,RefTetrahedron,1}()),
                           (7,Triangle, Lagrange{2,RefTetrahedron,2}()),
                           (5,Triangle, Lagrange{2,RefTetrahedron,3}()),
                           #(21,Tetrahedron, Lagrange{3,RefTetrahedron,1}()), #converges rather slowly in the gradient
                           (7,Tetrahedron, Lagrange{3,RefTetrahedron,2}()),
                           (21,Quadrilateral, Lagrange{2,RefCube,1}()),
                           (7,Quadrilateral, Lagrange{2,RefCube,2}()),
                           #(21,Hexahedron, Lagrange{3,RefCube,1}()), # slows down the pipeline quite a bit, so left out
                           (7,Hexahedron, Lagrange{3,RefCube,2}())]
        @testset "($size, $geo, $ip)" begin
            # Compute solution
            dh, u = manufactured_heat_problem(geo, ip, size)
            ip_geo = Ferrite.default_interpolation(typeof(Ferrite.getcells(FerriteViz.getgrid(dh), 1)))
            dim = Ferrite.getdim(ip)
            qr = QuadratureRule{dim, Ferrite.getrefshape(ip)}(1)

            # Check solution
            cellvalues = CellScalarValues(qr, Ferrite.getfieldinterpolation(dh, 1), ip_geo);
            for cell in CellIterator(dh)
                reinit!(cellvalues, cell)
                n_basefuncs = getnbasefunctions(cellvalues)
                coords = getcoordinates(cell)
                uₑ = u[celldofs(cell)]
                for q_point in 1:getnquadpoints(cellvalues)
                    x = spatial_coordinate(cellvalues, q_point, coords)
                    for i in 1:n_basefuncs
                        uₐₙₐ    = prod(cos, x*π/2)
                        uₐₚₚᵣₒₓ = function_value(cellvalues, q_point, uₑ)
                        @test isapprox(uₐₙₐ, uₐₚₚᵣₒₓ; atol=1e-2)
                    end
                end
            end

            # Compute gradient/flux field
            (dh_grad, u_grad) = FerriteViz.interpolate_gradient_field(dh, u, :u)

            # Check gradient of solution
            cellvalues_grad = CellVectorValues(qr, Ferrite.getfieldinterpolation(dh_grad, 1), ip_geo);
            for cell in CellIterator(dh_grad)
                reinit!(cellvalues_grad, cell)
                n_basefuncs = getnbasefunctions(cellvalues_grad)
                coords = getcoordinates(cell)
                uₑ = u_grad[celldofs(cell)]
                for q_point in 1:getnquadpoints(cellvalues_grad)
                    x = spatial_coordinate(cellvalues_grad, q_point, coords)
                    for i ∈ 1:n_basefuncs
                        uₐₚₚᵣₒₓ = function_value(cellvalues_grad, q_point, uₑ)
                        for d ∈ 1:dim
                            uₐₙₐ = π/2
                            for j ∈ 1:(d-1)
                                uₐₙₐ *= cos(x[j]*π/2)
                            end
                            uₐₙₐ *= -sin(x[d]*π/2)
                            for j ∈ (d+1):dim
                                uₐₙₐ *= cos(x[j]*π/2)
                            end
                            @test isapprox(uₐₙₐ, uₐₚₚᵣₒₓ[d]; atol=1e-1)
                        end
                    end
                end
            end
        end
    end
end

# TODO add gradient field test for vector valued problems via manufactured solution
