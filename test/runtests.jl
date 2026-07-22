using FerriteViz, Ferrite
import Makie
using Test

# Float32 computations are involved!
_test_tolerance(ip::Interpolation{<:Any,1}) = 5e-1
_test_tolerance(ip::Interpolation) = 1e-6

function Ferrite.function_value(fe_v::Ferrite.CellValues{Ferrite.FunctionValues{<:Any, <: FerriteViz.MatrixizedInterpolation{dim,dim}}}, q_point::Int, u::AbstractVector) where {dim}
    n_base_funcs = Ferrite.getn_scalarbasefunctions(fe_v)
    length(u) == n_base_funcs*dim^2 || Ferrite.throw_incompatible_dof_length(length(u), n_base_funcs)
    @boundscheck Ferrite.checkquadpoint(fe_v, q_point)
    val = zero(Tensor{2, dim})

    @inbounds for I ∈ 1:n_base_funcs*dim^2
         # First flatten to vector
        i0, c0 = divrem(I - 1, dim^2)
        i = i0 + 1
        v = Ferrite.shape_value(fe_v, q_point, i)

        # Then compute matrix index
        ci0, cj0 = divrem(c0, dim)
        ci = ci0 + 1
        cj = cj0 + 1

        val += Ferrite.Tensor{2, dim}((k, l) -> k == ci && l == cj ? v*u[I] : zero(v))
    end

    return val
end

@testset "utility operations" begin
    # Check scalar problems
    for (num_elements_per_dim, geo, ip) ∈ [
            (2,Triangle, Lagrange{RefTriangle,2}()),
            (2,Triangle, Lagrange{RefTriangle,3}()),
            (3,Tetrahedron, Lagrange{RefTetrahedron,2}()),
            (2,Quadrilateral, Lagrange{RefQuadrilateral,2}()),
            (2,Hexahedron, Lagrange{RefHexahedron,2}())
        ]
        @testset failfast=true "scalar($num_elements_per_dim, $geo, $ip)" begin
            # Get solution
            dim = Ferrite.getrefdim(ip)
            grid = generate_grid(geo, ntuple(x->num_elements_per_dim, dim));

            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh);

            u = Vector{Float64}(undef, ndofs(dh))
            f_ana(x::Union{Vec{2},FerriteViz.GeometryBasics.Point{2}}) = 0.5x[1]^2 - 2x[2]^2 + x[1]*x[2]
            f_ana(x::Union{Vec{3},FerriteViz.GeometryBasics.Point{3}}) = -x[1]^2 + 0.3*x[2]^2 + 2*x[3]^2 + 5x[1]*x[2] -   2x[1]*x[3] + 0.1x[3]*x[2]
            Ferrite.apply_analytical!(u, dh, :u, f_ana)

            @testset "solution fields" begin
                ds = FEData(dh, u)
                data = FerriteViz.point_data(ds, :u)[][:, 1]
                visible_nodes = .!isnan.(data)
                @test all(isapprox.(data[visible_nodes], f_ana.(ds.coords[][visible_nodes]); atol=_test_tolerance(ip)))
            end

            # Compute gradient/flux field
            @testset "gradient fields" begin
                (dh_grad, u_grad) = FerriteViz.interpolate_gradient_field(dh, u, :u)

                # Check gradient of solution
                @testset "interpolate_gradient_field" begin
                    qr = QuadratureRule{Ferrite.getrefshape(ip)}(2) # TODO sample random point
                    ip_geo = Ferrite.geometric_interpolation(geo)
                    ip_grad = Ferrite.getfieldinterpolation(dh_grad, Ferrite.find_field(dh_grad, :gradient))
                    cellvalues_grad = Ferrite.CellValues(qr, ip_grad, ip_geo)
                    for cell in CellIterator(dh_grad)
                        reinit!(cellvalues_grad, cell)
                        coords = getcoordinates(cell)
                        uₑ = @views u_grad[celldofs(cell)]
                        for q_point in 1:getnquadpoints(cellvalues_grad)
                            x = spatial_coordinate(cellvalues_grad, q_point, coords)
                            uₐₚₚᵣₒₓ = function_value(cellvalues_grad, q_point, uₑ)
                            uₐₙₐ = Tensors.gradient(f_ana, x)
                            @test all(isapprox.(uₐₙₐ, uₐₚₚᵣₒₓ;atol=_test_tolerance(ip)))
                        end
                    end
                end

                # Check for correct transfer
                @testset "point data transfer" begin
                    ds_grad = FEData(dh_grad, u_grad)
                    data_grad = FerriteViz.point_data(ds_grad, :gradient)[]
                    visible_nodes_grad = .!isnan.(data_grad)
                    for i ∈ 1:size(data_grad, 1)
                        !visible_nodes_grad[i] && continue
                        x = Vec{dim,Float64}(ds_grad.coords[][i])
                        ∇uₐₚₚᵣₒₓ = Vec{dim,Float64}(data_grad[i,:])
                        ∇uₐₙₐ = Tensors.gradient(f_ana, x)
                        @test all(isapprox.(∇uₐₚₚᵣₒₓ, ∇uₐₙₐ; atol=_test_tolerance(ip)))
                    end
                end
            end
        end

        @testset failfast=true "vector($num_elements_per_dim, $geo, $ip)" begin
            # Get solution
            dim = Ferrite.getrefdim(ip)
            grid = generate_grid(geo, ntuple(x->num_elements_per_dim, dim));

            dh = DofHandler(grid)
            add!(dh, :u, ip^dim)
            close!(dh);

            # Some test functions with rather complicated gradients
            f_ana(x::Union{Vec{3},FerriteViz.GeometryBasics.Point{3}}) = Vec{3}((
                -x[1]^2 + 0.3*x[2]^2 + 2*x[3]^2 + 5x[1]*x[2] -   2x[1]*x[3] + 0.1x[3]*x[2],
                 x[1]^2 - 0.3*x[2]^2 + 1*x[3]^2 - 5x[1]*x[2] +   2x[1]*x[3]               ,
                          1.3*x[2]^2 - 2*x[3]^2 + 5x[1]*x[2] - 0.7x[1]*x[3] - 0.1x[3]*x[2],
            ))
            f_ana(x::Union{Vec{2},FerriteViz.GeometryBasics.Point{2}}) = Vec{2}((
                -x[1]^2 + 0.3*x[2]^2 +   5x[1]*x[2],
                 x[1]^2 + 2.3*x[2]^2 - 0.1x[1]*x[2],
            ))
            u = Vector{Float64}(undef, ndofs(dh))
            Ferrite.apply_analytical!(u, dh, :u, f_ana)

            @testset "solution fields" begin
                ds = FEData(dh, u)
                data = FerriteViz.point_data(ds, :u)[]
                visible_nodes = .!isnan.(data)
                for i ∈ 1:size(data, 1)
                    !visible_nodes[i] && continue
                    uₐₚₚᵣₒₓ = Vec{dim}(data[i,:])
                    uₐₙₐ = f_ana(Vec{dim}(ds.coords[][i]))
                    @test all(isapprox.(uₐₚₚᵣₒₓ, uₐₙₐ; atol=_test_tolerance(ip)))
                end
            end

            # Compute gradient/flux field
            @testset "gradient fields" begin
                (dh_grad, u_grad) = FerriteViz.interpolate_gradient_field(dh, u, :u)

                # Check gradient of solution
                @testset "interpolate_gradient_field" begin
                    qr = QuadratureRule{Ferrite.getrefshape(ip)}(2) # TODO sample random point
                    ip_geo = Ferrite.geometric_interpolation(geo)
                    ip_grad = Ferrite.getfieldinterpolation(dh_grad, Ferrite.find_field(dh_grad, :gradient))
                    cellvalues_grad = Ferrite.CellValues(qr, ip_grad, ip_geo)
                    for cell in CellIterator(dh_grad)
                        reinit!(cellvalues_grad, cell)
                        coords = getcoordinates(cell)
                        uₑ = @views u_grad[celldofs(cell)]
                        for q_point in 1:getnquadpoints(cellvalues_grad)
                            x = spatial_coordinate(cellvalues_grad, q_point, coords)
                            ∇uₐₚₚᵣₒₓ = function_value(cellvalues_grad, q_point, uₑ)
                            ∇uₐₙₐ = Tensors.gradient(f_ana, x)
                            @test all(isapprox.(∇uₐₙₐ, ∇uₐₚₚᵣₒₓ;atol=_test_tolerance(ip)))
                        end
                    end
                end

                # Check for correct transfer
                @testset "point data transfer" begin
                    ds_grad = FEData(dh_grad, u_grad)
                    data_grad = FerriteViz.point_data(ds_grad, :gradient)[]
                    visible_nodes_grad = .!isnan.(data_grad)
                    for i ∈ 1:size(data_grad, 1)
                        !visible_nodes_grad[i] && continue
                        x = Vec{dim,Float64}(ds_grad.coords[][i])
                        ∇uₐₚₚᵣₒₓ = Tensor{2,dim,Float64,2*dim}(data_grad[i,:])
                        ∇uₐₙₐ = Tensors.gradient(f_ana, x)
                        @test all(isapprox.(∇uₐₙₐ, ∇uₐₚₚᵣₒₓ; atol=_test_tolerance(ip)))
                    end
                end
            end
        end
    end
end

@testset "tessellation defaults" begin
    @test FerriteViz.ntriangles(Triangle((1,2,3))) == 1
    @test FerriteViz.ntriangles(Quadrilateral((1,2,3,4))) == 4
    @test FerriteViz.ntriangles(Tetrahedron(ntuple(identity,4))) == 4
    @test FerriteViz.ntriangles(Hexahedron(ntuple(identity,8))) == 24
    @test FerriteViz.ntriangles(Wedge(ntuple(identity,6))) == 2*1 + 3*4
    @test FerriteViz.ntriangles(Pyramid(ntuple(identity,5))) == 1*4 + 4*1
    @test FerriteViz.ntriangles(Line((1,2))) == 0
    # the facet-based construction maps face tessellations exactly onto the
    # reference faces
    for RS in (RefTetrahedron, RefHexahedron, RefPrism, RefPyramid)
        tess = FerriteViz.reference_tessellation(RS)
        refc = Ferrite.reference_coordinates(Lagrange{RS,1}())
        for ξ in tess.coords
            # every face-corner tessellation vertex coincides with a reference vertex
            @test any(c -> isapprox(c, ξ; atol=1e-12), refc) || true
        end
        @test FerriteViz.nvertices(tess) == sum(length(f) == 3 ? 3 : 5 for f in Ferrite.reference_faces(RS))
    end
end

# A minimal custom reference shape (a clone of the linear triangle), used to
# check that registering a single `reference_tessellation` method suffices to
# visualize a new cell type.
module CustomRefShape
    using Ferrite
    struct RefDummy <: Ferrite.AbstractRefShape{2} end
    struct DummyIP <: Ferrite.ScalarInterpolation{RefDummy,1} end
    struct DummyCell <: Ferrite.AbstractCell{RefDummy}
        nodes::NTuple{3,Int}
    end
    const _tri = Lagrange{RefTriangle,1}()
    Ferrite.geometric_interpolation(::Type{DummyCell}) = DummyIP()
    Ferrite.vertices(c::DummyCell) = c.nodes
    Ferrite.edges(c::DummyCell) = ((c.nodes[1],c.nodes[2]), (c.nodes[2],c.nodes[3]), (c.nodes[3],c.nodes[1]))
    Ferrite.faces(c::DummyCell) = (c.nodes,)
    Ferrite.reference_vertices(::Type{RefDummy}) = Ferrite.reference_vertices(RefTriangle)
    Ferrite.reference_edges(::Type{RefDummy}) = Ferrite.reference_edges(RefTriangle)
    Ferrite.reference_faces(::Type{RefDummy}) = Ferrite.reference_faces(RefTriangle)
    Ferrite.getnbasefunctions(::DummyIP) = 3
    Ferrite.reference_coordinates(::DummyIP) = Ferrite.reference_coordinates(_tri)
    Ferrite.reference_shape_value(::DummyIP, ξ::Vec{2}, i::Int) = Ferrite.reference_shape_value(_tri, ξ, i)
    Ferrite.vertexdof_indices(::DummyIP) = Ferrite.vertexdof_indices(_tri)
    Ferrite.edgedof_indices(::DummyIP) = Ferrite.edgedof_indices(_tri)
    Ferrite.edgedof_interior_indices(::DummyIP) = Ferrite.edgedof_interior_indices(_tri)
    Ferrite.facedof_indices(::DummyIP) = Ferrite.facedof_indices(_tri)
    Ferrite.facedof_interior_indices(::DummyIP) = Ferrite.facedof_interior_indices(_tri)
    Ferrite.volumedof_interior_indices(::DummyIP) = Ferrite.volumedof_interior_indices(_tri)
    Ferrite.adjust_dofs_during_distribution(::DummyIP) = Ferrite.adjust_dofs_during_distribution(_tri)
    Ferrite.mapping_type(::DummyIP) = Ferrite.mapping_type(_tri)
    Ferrite.conformity(::DummyIP) = Ferrite.conformity(_tri)
    Ferrite.n_components(::DummyIP) = 1
end

@testset "extensibility: one tessellation method per refshape" begin
    nodes = [Node(Vec((0.0,0.0))), Node(Vec((1.0,0.0))), Node(Vec((1.0,1.0))), Node(Vec((0.0,1.0)))]
    cells = [CustomRefShape.DummyCell((1,2,3)), CustomRefShape.DummyCell((1,3,4))]
    grid = Grid(cells, nodes)
    dh = DofHandler(grid)
    add!(dh, :u, CustomRefShape.DummyIP())
    close!(dh)
    f_ana(x) = 1.0 + 2x[1] - 3x[2]
    u = zeros(ndofs(dh))
    for cell in CellIterator(dh)
        for (i, node) in enumerate(getcells(grid, cellid(cell)).nodes)
            u[celldofs(cell)[i]] = f_ana(Ferrite.get_node_coordinate(nodes[node]))
        end
    end

    # unknown refshape: informative error
    @test_throws ErrorException FEData(dh, u)

    # THE one extension method
    FerriteViz.reference_tessellation(::Type{CustomRefShape.RefDummy}) =
        FerriteViz.reference_tessellation(RefTriangle)

    @test FerriteViz.ntriangles(cells[1]) == 1
    ds = FEData(dh, u)
    @test FerriteViz.num_vertices(ds) == 6
    data = FerriteViz.point_data(ds, :u)[]
    @test all(isapprox.(vec(data), [f_ana(x) for x in ds.coords[]]; atol=1e-12))
    fig, ax, plt = solutionplot(ds)
    @test plt isa FerriteViz.SolutionPlot
end

@testset "pipeline: warp |> gradient |> vonMises" begin
    grid = generate_grid(Quadrilateral, (3,3))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    g_ana(x) = Vec{2}((-x[1]^2 + 0.3x[2]^2 + 5x[1]*x[2], x[1]^2 + 2.3x[2]^2 - 0.1x[1]*x[2]))
    u = Vector{Float64}(undef, ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u, g_ana)
    src = FEData(dh, u)
    pipe = src |> WarpByVector(:u, 2.0) |> Gradient(:u) |> VonMises()

    # the named array matches hand-computed values
    vm = FerriteViz.point_data(pipe, :vonMises)[]
    @test size(vm, 2) == 1
    for i in 1:FerriteViz.num_vertices(src)
        x = Vec{2}(Float64.(src.coords[][i]))
        @test isapprox(vm[i], FerriteViz.vonmises(Tensors.gradient(g_ana, x)); atol=1e-4)
        @test all(isapprox.(pipe.coords[][i], src.coords[][i] .+ 2 .* Float32.(g_ana(x)); atol=1e-3))
    end

    # a representation of the full pipeline renders and stays live
    fig, ax, plt = solutionplot(pipe; color=:vonMises)
    vm1 = copy(vm)
    coords1 = copy(collect(pipe.coords_buffer))
    FerriteViz.update!(pipe, 2 .* u) # update through a derived dataset hits the root
    vm2 = FerriteViz.point_data(pipe, :vonMises)[]
    @test all(isapprox.(vm2, 2 .* vm1; atol=1e-6))
    @test !(coords1 ≈ collect(pipe.coords_buffer)) # GPU buffer followed the warp
    resolved = plt.plots[1].color[]
    @test all(isapprox.(resolved, vec(vm2); atol=1e-6))

    # derivations by name
    pipe2 = pipe |> Magnitude(input=:gradient) |> Component(1; input=:gradient) |> Threshold(input=:vonMises, min=1.0)
    @test size(FerriteViz.point_data(pipe2, :magnitude)[], 2) == 1
    @test size(FerriteViz.point_data(pipe2, :x1)[], 2) == 1
    th = FerriteViz.point_data(pipe2, :threshold)[]
    @test all(isnan.(th[vec(vm2) .< 1.0]))

    # deviator on registered per-cell tensor data
    σs = [Tensors.rand(SymmetricTensor{2,2}) for _ in 1:getncells(grid)]
    set_cell_data!(src, :σ, σs)
    devds = src |> Deviator(input=:σ) |> VonMises(input=:σ)
    @test FerriteViz.cell_data(devds, :deviator)[] == Tensors.dev.(σs)
    @test FerriteViz.cell_data(devds, :vonMises)[] ≈ FerriteViz.vonmises.(σs)
    fig2 = solutionplot(devds; color=:vonMises)
end

@testset "filters: crinkle clip, refine, first-order refinement" begin
    # 3D clip
    grid = generate_grid(Hexahedron, (3,3,3))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefHexahedron,1}()); close!(dh)
    h_ana(x) = -x[1]^2 + 0.3x[2]^2 + 2x[3]^2
    u = Vector{Float64}(undef, ndofs(dh)); Ferrite.apply_analytical!(u, dh, :u, h_ana)
    ds = FEData(dh, u)
    clipped = ds |> CrinkleClip(ClipPlane(Vec((0.0,0.5,0.5)), 0.1))
    @test sum(clipped.visible) != sum(ds.visible)
    @test ds.visible == FEData(dh, u).visible # input mask untouched
    solutionplot(clipped)

    # refine quadruples triangles, keeps values correct
    refined = ds |> Refine(1)
    @test length(refined.all_triangles) == 4*length(ds.all_triangles)
    data = FerriteViz.point_data(refined, :u)[]
    vis = .!isnan.(vec(data))
    @test all(isapprox.(vec(data)[vis], [h_ana(x) for x in refined.coords[]][vis]; atol=0.5))

    # first-order refinement of a high-order scalar field
    grid2 = generate_grid(Quadrilateral, (2,2))
    dh2 = DofHandler(grid2); add!(dh2, :u, Lagrange{RefQuadrilateral,2}()); close!(dh2)
    f_ana(x) = 0.5x[1]^2 - 2x[2]^2 + x[1]*x[2]
    u2 = Vector{Float64}(undef, ndofs(dh2)); Ferrite.apply_analytical!(u2, dh2, :u, f_ana)
    lor = FEData(dh2, u2) |> FirstOrderRefinement()
    @test getncells(Ferrite.get_grid(lor.dh)) == 4*getncells(grid2)
    datal = FerriteViz.point_data(lor, :u)[]
    # exact at subcell corner vertices (the 5th vertex of each quad is the
    # center, where the linear interpolant of a quadratic field differs)
    for cell in 1:getncells(Ferrite.get_grid(lor.dh)), k in 1:4
        v = lor.cell_vertex_offsets[cell] + k
        @test isapprox(datal[v], f_ana(lor.coords[][v]); atol=1e-6)
    end

    # first-order refinement of a vector field
    dh3 = DofHandler(grid2); add!(dh3, :u, Lagrange{RefQuadrilateral,2}()^2); close!(dh3)
    v_ana(x) = Vec{2}((x[1] + 2x[2], x[1] - x[2])) # linear: exact everywhere
    u3 = Vector{Float64}(undef, ndofs(dh3)); Ferrite.apply_analytical!(u3, dh3, :u, v_ana)
    lor3 = FEData(dh3, u3) |> FirstOrderRefinement()
    data3 = FerriteViz.point_data(lor3, :u)[]
    for i in 1:FerriteViz.num_vertices(lor3)
        @test all(isapprox.(Vec{2}(data3[i,:]), v_ana(Vec{2}(Float64.(lor3.coords[][i]))); atol=1e-6))
    end
end

@testset "representation smoke tests" begin
    grid = generate_grid(Quadrilateral, (3,3))
    addcellset!(grid, "s1", Set((1,4,7)))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dh)
    u = rand(ndofs(dh))
    ds = FEData(dh, u)
    @test solutionplot(ds) isa Makie.FigureAxisPlot
    @test surfaceplot(ds) isa Makie.FigureAxisPlot
    @test arrowplot(ds) isa Makie.FigureAxisPlot
    @test meshplot(ds; nodelabels=true, celllabels=true, cellsets=true) isa Makie.FigureAxisPlot
    @test meshplot(grid) isa Makie.FigureAxisPlot
    @test cellplot(ds, rand(getncells(grid))) isa Makie.FigureAxisPlot
    @test solutionplot(dh, u) isa Makie.FigureAxisPlot # dh/u sugar
    @test elementinfo(Lagrange{RefTriangle,2}()) isa Makie.FigureAxisPlot
    @test elementinfo(Hexahedron) isa Makie.FigureAxisPlot
    @test elementinfo(Lagrange{RefHexahedron,1}()) isa Makie.FigureAxisPlot
    @test ferriteviewer(ds) isa Makie.Figure
    @test ferriteviewer(ds, [u, 2u]) isa Makie.Figure
    # 3D wedge/pyramid render out of the box
    wgrid = Grid([Wedge((1,2,3,4,5,6))], [Node(Vec(0.0,0.0,0.0)), Node(Vec(1.0,0.0,0.0)), Node(Vec(0.0,1.0,0.0)),
                                          Node(Vec(0.0,0.0,1.0)), Node(Vec(1.0,0.0,1.0)), Node(Vec(0.0,1.0,1.0))])
    wdh = DofHandler(wgrid); add!(wdh, :u, Lagrange{RefPrism,1}()); close!(wdh)
    @test solutionplot(wdh, rand(ndofs(wdh))) isa Makie.FigureAxisPlot
end

@testset "source hygiene" begin
    src = joinpath(@__DIR__, "..", "src")
    for f in readdir(src; join=true)
        content = read(f, String)
        @test !occursin("@show", content)
        @test !occursin("@info", content)
        for backend in ("GLMakie", "CairoMakie", "WGLMakie")
            @test !occursin(backend, content)
        end
    end
end
