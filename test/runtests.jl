using FerriteViz, Ferrite
using Test

_test_tolerance(ip::Interpolation{<:Any,1}) = 5e-1
_test_tolerance(ip::Interpolation) = 1e-8

@testset "gradient fields" begin
    # Check scalar problems
    for (num_elements_per_dim, geo, ip) ∈ [
                           (4, Triangle, Lagrange{RefTriangle,1}()),
                           (2,Triangle, Lagrange{RefTriangle,2}()),
                           (2,Triangle, Lagrange{RefTriangle,3}()),
                           (5,Tetrahedron, Lagrange{RefTetrahedron,1}()),
                           (3,Tetrahedron, Lagrange{RefTetrahedron,2}()),
                           (4,Quadrilateral, Lagrange{RefQuadrilateral,1}()),
                           (2,Quadrilateral, Lagrange{RefQuadrilateral,2}()),
                           (4,Hexahedron, Lagrange{RefHexahedron,1}()),
                           (2,Hexahedron, Lagrange{RefHexahedron,2}())
        ]
        @testset failfast=true "scalar($num_elements_per_dim, $geo, $ip)" begin
            # Get solution
            dim = Ferrite.getdim(ip)
            grid = generate_grid(geo, ntuple(x->num_elements_per_dim, dim));

            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh);

            u = Vector{Float64}(undef, ndofs(dh))
            f_ana(x) = sum(0.5 * x.^2)
            Ferrite.apply_analytical!(u, dh, :u, f_ana)

            # Compute gradient/flux field
            (dh_grad, u_grad) = FerriteViz.interpolate_gradient_field(dh, u, :u)

            # Check gradient of solution
            qr = QuadratureRule{Ferrite.getrefshape(ip)}(2)
            ip_geo = Ferrite.default_interpolation(geo)
            ip_grad = Ferrite.getfieldinterpolation(dh_grad, Ferrite.find_field(dh_grad, :gradient))
            cellvalues_grad = Ferrite.PointValuesInternal(qr.points[1], ip_grad)
            cellvalues_geo = Ferrite.PointValuesInternal(qr.points[1], ip_geo)
            for cell in CellIterator(dh_grad)
                coords = getcoordinates(cell)
                uₑ = u_grad[celldofs(cell)]
                for q_point in 1:getnquadpoints(cellvalues_grad)
                    x = function_value(cellvalues_geo, q_point, coords)
                    uₐₚₚᵣₒₓ = function_value(cellvalues_grad, q_point, uₑ)
                    uₐₙₐ = Tensors.gradient(f_ana, x)
                    @test isapprox(uₐₙₐ, uₐₚₚᵣₒₓ;atol=_test_tolerance(ip))
                end
            end
        end

        @testset "vector($num_elements_per_dim, $geo, $ip)" begin
            # Get solution
            dim = Ferrite.getdim(ip)
            grid = generate_grid(geo, ntuple(x->num_elements_per_dim, dim));

            dh = DofHandler(grid)
            add!(dh, :u, ip^dim)
            close!(dh);

            f_ana(x) = Vec{dim}(i->(sum(x.^(i-1) /i)))
            u = Vector{Float64}(undef, ndofs(dh))
            Ferrite.apply_analytical!(u, dh, :u, f_ana)

            # Compute gradient/flux field
            (dh_grad, u_grad) = FerriteViz.interpolate_gradient_field(dh, u, :u)

            # Check gradient of solution
            qr = QuadratureRule{Ferrite.getrefshape(ip)}(2)
            ip_geo = Ferrite.default_interpolation(geo)
            ip_grad = Ferrite.getfieldinterpolation(dh_grad, Ferrite.find_field(dh_grad, :gradient))
            cellvalues_grad = Ferrite.PointValuesInternal(qr.points[1], ip_grad)
            cellvalues_geo = Ferrite.PointValuesInternal(qr.points[1], ip_geo)
            for cell in CellIterator(dh_grad)
                coords = getcoordinates(cell)
                uₑ = u_grad[celldofs(cell)]
                for q_point in 1:getnquadpoints(cellvalues_grad)
                    x = function_value(cellvalues_geo, q_point, coords)
                    uₐₚₚᵣₒₓ = function_value(cellvalues_grad, q_point, uₑ)
                    uₐₙₐ = Tensors.gradient(f_ana, x)
                    @test isapprox(uₐₙₐ, uₐₚₚᵣₒₓ;atol=_test_tolerance(ip))
                end
            end
        end
    end
end
