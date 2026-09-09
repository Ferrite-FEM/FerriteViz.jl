using FerriteViz, Ferrite
import Makie
import GeometryBasics
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

include("isubd.jl")

@testset "tessellation defaults" begin
    @test FerriteViz.ntriangles(Triangle((1,2,3))) == 1
    @test FerriteViz.ntriangles(Quadrilateral((1,2,3,4))) == 4
    @test FerriteViz.ntriangles(Tetrahedron(ntuple(identity,4))) == 4
    @test FerriteViz.ntriangles(Hexahedron(ntuple(identity,8))) == 24
    @test FerriteViz.ntriangles(Wedge(ntuple(identity,6))) == 2*1 + 3*4
    @test FerriteViz.ntriangles(Pyramid(ntuple(identity,5))) == 1*4 + 4*1
    @test FerriteViz.ntriangles(Line((1,2))) == 0
    # the facet-based construction maps each face's corners exactly onto the
    # reference-face vertices, preserving orientation
    for RS in (RefTetrahedron, RefHexahedron, RefPrism, RefPyramid)
        tess = FerriteViz.reference_tessellation(RS)
        refc = Ferrite.reference_coordinates(Lagrange{RS,1}())
        for (fi, face) in enumerate(Ferrite.reference_faces(RS))
            face2d = Ferrite.reference_coordinates(length(face) == 3 ? Lagrange{RefTriangle,1}() : Lagrange{RefQuadrilateral,1}())
            for (k, v) in enumerate(face)
                @test Ferrite.facet_to_element_transformation(face2d[k], RS, fi) ≈ refc[v] atol=1e-12
            end
        end
        # face vertices plus 2 dedicated endpoint vertices per wireframe edge
        @test FerriteViz.nvertices(tess) == sum(length(f) == 3 ? 3 : 5 for f in Ferrite.reference_faces(RS)) +
                                            2 * length(Ferrite.reference_edges(RS))
        @test FerriteViz.nedges(tess) == length(Ferrite.reference_edges(RS))
        # every edge connects the reference coordinates of its FE edge's vertices
        for (e, (v1, v2)) in enumerate(Ferrite.reference_edges(RS))
            a, b = tess.edges[e]
            @test tess.coords[a] ≈ refc[v1] && tess.coords[b] ≈ refc[v2]
        end
    end
end

@testset "reference-space subdivision" begin
    for RS in (RefTriangle, RefQuadrilateral, RefTetrahedron, RefHexahedron)
        base = FerriteViz.reference_tessellation(RS)
        for n in 1:2
            tess = FerriteViz.subdivide(base, n)
            @test FerriteViz.ntriangles(tess) == 4^n * FerriteViz.ntriangles(base)
            @test FerriteViz.nedges(tess) == 2^n * FerriteViz.nedges(base)
            # subdivided edge segments join up into the original edge polylines:
            # each original edge's segments cover it with 2^n equal pieces
            for (e, (a, b)) in enumerate(base.edges)
                segs = tess.edges[(2^n*(e-1)+1):(2^n*e)]
                @test tess.coords[segs[1][1]] ≈ base.coords[a]
                @test tess.coords[segs[end][2]] ≈ base.coords[b]
                for k in 1:(length(segs)-1) # consecutive segments share their joint
                    @test segs[k][2] == segs[k+1][1]
                end
            end
        end
        # midpoints are deduplicated: neighbouring triangles share them
        t1 = FerriteViz.subdivide(base, 1)
        naive = FerriteViz.nvertices(base) + 3 * FerriteViz.ntriangles(base) + FerriteViz.nedges(base)
        @test FerriteViz.nvertices(t1) < naive
    end
    # edge-only rounds refine the wireframe without touching the surface
    base = FerriteViz.reference_tessellation(RefQuadrilateral)
    fine = FerriteViz._subdivided_tessellation(base, 0, 2)
    @test FerriteViz.ntriangles(fine) == FerriteViz.ntriangles(base)
    @test FerriteViz.nedges(fine) == 4 * FerriteViz.nedges(base)
end

@testset "wireframe edges follow the pipeline" begin
    # curved wireframe: a quadratic displacement on a linear grid bends edges
    grid = generate_grid(Quadrilateral, (1,1))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2); close!(dh)
    g_ana(x) = Vec{2}((0.0, x[1]^2))
    u = zeros(ndofs(dh)); Ferrite.apply_analytical!(u, dh, :u, g_ana)
    # the constructor builds the flat base; Refine()'s automatic mode picks
    # the subdivision per cell type (3 edge rounds for order 2)
    ds = FEData(dh, u; adaptivity=false) |> Refine()
    @test length(ds.all_edges) == 4 * 2^3
    warped = ds |> WarpByVector(:u)
    pts = warped.coords[][FerriteViz._visible_edge_indices(warped)]
    # every warped edge vertex satisfies y = y₀ + x² exactly (up to Float32);
    # a straight-chord wireframe would interpolate linearly between corners
    for (p0, p) in zip(ds.coords[][FerriteViz._visible_edge_indices(ds)], pts)
        @test isapprox(p[2], p0[2] + p0[1]^2; atol=1e-5)
    end
    @test meshplot(warped) isa Makie.FigureAxisPlot

    # the constructor's static tessellation is the flat base
    ds0 = FEData(dh, u)
    @test length(ds0.all_triangles) == 4 && length(ds0.all_edges) == 4
    # ... and the Refine filter subdivides it with explicit control
    @test length((ds0 |> Refine(1)).all_triangles) == 16
    dsl = FEData(DofHandler(grid), Float64[]) |> Refine(1) # linear grid, forced
    @test length(dsl.all_triangles) == 16
    @test length((ds0 |> Refine(surface=0, edges=2)).all_edges) == 4 * 4
    # the automatic filter mode matches the explicit pipe above
    dsauto = ds0 |> Refine()
    @test length(dsauto.all_triangles) == length(ds.all_triangles)
    @test dsauto.reference_coords == ds.reference_coords
    # data across the rebuild: cell data survives, point data is dropped,
    # dof fields transfer onto the new vertices
    set_cell_data!(ds0, :c, [1.0])
    set_point_data!(ds0, :pd, ones(FerriteViz.num_vertices(ds0)))
    sub = ds0 |> Refine(1)
    @test FerriteViz.cell_data(sub, :c)[] == [1.0]
    @test !haskey(sub.point_data, :pd)
    @test size(FerriteViz.point_data(sub, :u)[], 1) == FerriteViz.num_vertices(sub)

    # 3D: the wireframe is restricted to the visible cells, so a crinkle clip
    # hides the clipped cells' edges
    grid3 = generate_grid(Hexahedron, (3,3,3))
    dh3 = DofHandler(grid3); add!(dh3, :u, Lagrange{RefHexahedron,1}()); close!(dh3)
    ds3 = FEData(dh3, rand(ndofs(dh3)))
    @test length(ds3.all_edges) == getncells(grid3) * 12
    clipped = ds3 |> CrinkleClip(ClipPlane(Vec((0.0,0.5,0.5)), 0.1))
    @test length(FerriteViz._visible_edge_indices(clipped)) != length(FerriteViz._visible_edge_indices(ds3))
    @test 2 * length(ds3.all_edges) > length(FerriteViz._visible_edge_indices(ds3)) # interior cells hidden
    @test meshplot(clipped) isa Makie.FigureAxisPlot
    # node markers/labels follow the visibility too
    @test length(FerriteViz._visible_node_ids(clipped)) < length(FerriteViz._visible_node_ids(ds3))
    @test length(FerriteViz._visible_node_ids(FEData(DofHandler(grid), Float64[]))) == getnnodes(grid)
    # Refine preserves the input's visibility (e.g. downstream of a clip)
    subclipped = clipped |> Refine(1)
    @test subclipped.visible == clipped.visible
    @test length(subclipped.all_triangles) == 4 * length(clipped.all_triangles)
    @test meshplot(subclipped) isa Makie.FigureAxisPlot

    # Refine splits the wireframe with the triangles
    refined = ds3 |> Refine(1)
    @test length(refined.all_edges) == 2 * length(ds3.all_edges)
    @test all(1:getncells(grid3)) do c
        length(FerriteViz.edges_on_cell(refined, c)) == 2 * length(FerriteViz.edges_on_cell(ds3, c))
    end
    @test meshplot(refined) isa Makie.FigureAxisPlot

    # AddQuadraturePointData rebuilds the geometry but keeps the FE edges
    qr = QuadratureRule{RefQuadrilateral}(2)
    gridq = generate_grid(Quadrilateral, (3,3))
    dhq = DofHandler(gridq); add!(dhq, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dhq)
    dsq = FEData(dhq, rand(ndofs(dhq)))
    vals = [[Float64(10c + q) for q in 1:getnquadpoints(qr)] for c in 1:getncells(gridq)]
    pipe = dsq |> AddQuadraturePointData(qr, vals; output=:iv)
    @test length(pipe.all_edges) == getncells(gridq) * 4
    @test meshplot(pipe) isa Makie.FigureAxisPlot
    @test meshplot(pipe |> WarpByVector(:u, 0.1)) isa Makie.FigureAxisPlot
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
    # this testset checks the *static* resolution end to end
    src = FEData(dh, u; adaptivity=false)
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
    pipe2 = pipe |> Magnitude(input=:gradient) |> ExtractComponent(1; input=:gradient) |> Threshold(input=:vonMises, min=1.0)
    @test size(FerriteViz.point_data(pipe2, :magnitude)[], 2) == 1
    @test size(FerriteViz.point_data(pipe2, :x1)[], 2) == 1
    th = FerriteViz.point_data(pipe2, :threshold)[]
    @test all(isnan.(th[vec(vm2) .< 1.0]))

    # Derive over several inputs: one argument per name, each wrapped on its own
    nv = FerriteViz.num_vertices(src)
    set_point_data!(src, :a, reshape(collect(1.0:nv), nv, 1))
    set_point_data!(src, :b, reshape(collect(1.0:nv) .* 10, nv, 1))
    multi = src |> Derive((x, y) -> x + y; input=[:a, :b], output=:sum)
    @test vec(FerriteViz.point_data(multi, :sum)[]) ≈ collect(1.0:nv) .* 11
    # a single Symbol keeps working and is equivalent to a one-element vector
    @test vec(FerriteViz.point_data(src |> Derive(x -> 2x; input=:a, output=:d), :d)[]) ≈
          vec(FerriteViz.point_data(src |> Derive(x -> 2x; input=[:a], output=:d), :d)[])
    # inputs of different component counts are wrapped independently
    mixed = src |> Gradient(:u; copy_fields=[:u]) |>
            Derive((∇u, uu) -> norm(∇u) * norm(uu); input=[:gradient, :u], output=:mix)
    G = FerriteViz.point_data(mixed, :gradient)[]; U = FerriteViz.point_data(mixed, :u)[]
    @test vec(FerriteViz.point_data(mixed, :mix)[]) ≈
          [norm(FerriteViz._wrap_row(view(G, i, :), 2)) *
           norm(FerriteViz._wrap_row(view(U, i, :), 2)) for i in axes(G, 1)]
    # multiple cell-data inputs
    ncell = getncells(grid)
    set_cell_data!(src, :ca, collect(1.0:ncell))
    set_cell_data!(src, :cb, collect(1.0:ncell) .* 5)
    @test FerriteViz.cell_data(src |> Derive((x, y) -> x * y; input=[:ca, :cb], output=:cp), :cp)[] ≈
          collect(1.0:ncell) .^ 2 .* 5
    @test_throws ErrorException src |> Derive((x, y) -> x + y; input=[:a, :nope], output=:z)
    @test_throws ErrorException src |> Derive((x, y) -> x + y; input=[:a, :ca], output=:z) # point/cell mix

    # deviator on registered per-cell tensor data
    σs = [Tensors.rand(SymmetricTensor{2,2}) for _ in 1:getncells(grid)]
    set_cell_data!(src, :σ, σs)
    devds = src |> Deviator(input=:σ) |> VonMises(input=:σ)
    @test FerriteViz.cell_data(devds, :deviator)[] == Tensors.dev.(σs)
    @test FerriteViz.cell_data(devds, :vonMises)[] ≈ FerriteViz.vonmises.(σs)
    fig2 = solutionplot(devds; color=:vonMises)
    fig3 = cellplot(devds; color=:vonMises) # named cell-data variant

    # warp by a registered (non-dof) point-data array: tessellation vertices
    # move, grid nodes stay
    disp = ones(FerriteViz.num_vertices(src), 2)
    set_point_data!(src, :d, disp)
    warped = src |> WarpByVector(:d)
    @test warped.coords[] ≈ [c .+ Float32.((1, 1)) for c in src.coords[]]
    @test warped.gridnodes[] == src.gridnodes[]

    # named data cannot shadow dof fields
    @test_throws ErrorException set_point_data!(src, :u, disp)
    @test_throws ErrorException set_cell_data!(src, :u, zeros(getncells(grid)))
end

@testset "quadrature point data (Voronoi)" begin
    # area of a triangle given in reference space (2D or 3D)
    function _triarea(p, q, r)
        u = q - p; v = r - p
        length(p) == 2 && return abs(u[1] * v[2] - u[2] * v[1]) / 2
        c = (u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1])
        return sqrt(sum(abs2, c)) / 2
    end
    function region_areas(t, nqp)
        a = zeros(nqp)
        for tri in t.triangles
            qp = t.vertex_qp[tri[1]]
            @test t.vertex_qp[tri[2]] == qp && t.vertex_qp[tri[3]] == qp # regions never mix
            a[qp] += _triarea(t.coords[tri[1]], t.coords[tri[2]], t.coords[tri[3]])
        end
        return a
    end

    # the reference partition is exact: it tiles the shape without gaps/overlaps
    qr = QuadratureRule{RefQuadrilateral}(2)
    tq = FerriteViz.qp_voronoi_tessellation(RefQuadrilateral, qr)
    @test length(tq.triangles) == 8                       # 4 quadrants, 2 triangles each
    @test all(a -> isapprox(a, 1.0; atol=1e-10), region_areas(tq, getnquadpoints(qr)))
    # each region stays on its quadrature point's side of the mid-lines
    for (i, ξ) in enumerate(Ferrite.getpoints(qr)), v in eachindex(tq.coords)
        tq.vertex_qp[v] == i || continue
        @test sign(tq.coords[v][1]) == sign(ξ[1]) || abs(tq.coords[v][1]) < 1e-12
        @test sign(tq.coords[v][2]) == sign(ξ[2]) || abs(tq.coords[v][2]) < 1e-12
    end
    qrh = QuadratureRule{RefHexahedron}(2)
    @test all(a -> isapprox(a, 3.0; atol=1e-10),
              region_areas(FerriteViz.qp_voronoi_tessellation(RefHexahedron, qrh), getnquadpoints(qrh)))
    qrt = QuadratureRule{RefTetrahedron}(2)
    @test sum(region_areas(FerriteViz.qp_voronoi_tessellation(RefTetrahedron, qrt), getnquadpoints(qrt))) ≈
          3 * 0.5 + sqrt(3) / 2 atol = 1e-10
    qr1 = QuadratureRule{RefTriangle}(1)                  # single point -> whole cell
    @test sum(region_areas(FerriteViz.qp_voronoi_tessellation(RefTriangle, qr1), 1)) ≈ 0.5 atol = 1e-12
    # prism/pyramid: the regions tile exactly the boundary surface drawn by the
    # reference tessellation
    for RS in (RefPrism, RefPyramid)
        qrv = QuadratureRule{RS}(2)
        surf = FerriteViz.reference_tessellation(RS)
        expected = sum(_triarea(surf.coords[tri[1]], surf.coords[tri[2]], surf.coords[tri[3]]) for tri in surf.triangles)
        @test sum(region_areas(FerriteViz.qp_voronoi_tessellation(RS, qrv), getnquadpoints(qrv))) ≈ expected atol = 1e-10
    end

    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dh)
    u = rand(ndofs(dh)); ds = FEData(dh, u)
    nqp = getnquadpoints(qr); ncells = getncells(grid)
    vals = [[Float64(10c + q) for q in 1:nqp] for c in 1:ncells]

    pipe = ds |> AddQuadraturePointData(qr, vals; output=:iv)
    A = FerriteViz.point_data(pipe, :iv)[]
    @test size(A) == (FerriteViz.num_vertices(pipe), 1)
    # piecewise constant: a cell shows exactly its quadrature point values
    for c in 1:ncells
        @test sort(unique(vec(A[collect(FerriteViz.vertices_on_cell(pipe, c)), 1]))) ≈ sort(vals[c])
    end
    # Matrix input is equivalent to the ragged Vector{Vector} form
    M = [Float64(10c + q) for c in 1:ncells, q in 1:nqp]
    @test FerriteViz.point_data(ds |> AddQuadraturePointData(qr, M; output=:iv), :iv)[] ≈ A

    # symmetric tensors expand to full components so VonMises composes
    symv = [[SymmetricTensor{2,2}((1.0c, 0.5q, 2.0)) for q in 1:nqp] for c in 1:ncells]
    pσ = ds |> AddQuadraturePointData(qr, symv; output=:σ) |> VonMises(input=:σ, output=:σvM)
    @test size(FerriteViz.point_data(pσ, :σ)[], 2) == 4
    @test all(isfinite, FerriteViz._scalar_data(pσ, :σvM)[])
    # extract pulls the value out of a state-like struct
    states = [[(a=v, b=0.0) for v in cell] for cell in vals]
    @test FerriteViz.point_data(ds |> AddQuadraturePointData(qr, states; output=:iv, extract=s -> s.a), :iv)[] ≈ A

    # geometry is rebuilt but dof fields stay usable, so warping still works after
    @test size(FerriteViz.point_data(pipe |> WarpByVector(:u, 2.0), :u)[], 2) == 2
    @test solutionplot(pipe; color=:iv) isa Makie.FigureAxisPlot

    # values are reactive
    obs = Makie.Observable(vals)
    pr = ds |> AddQuadraturePointData(qr, obs; output=:iv)
    before = copy(FerriteViz.point_data(pr, :iv)[])
    obs[] = [[3v for v in cell] for cell in vals]
    @test FerriteViz.point_data(pr, :iv)[] ≈ 3 .* before

    # two rules with the same quadrature rule share a vertex layout, so their
    # arrays survive and can be combined in one Derive
    both = ds |> AddQuadraturePointData(qr, vals; output=:a) |>
                 AddQuadraturePointData(qr, [[2v for v in c] for c in vals]; output=:b) |>
                 Derive((x, y) -> x + y; input=[:a, :b], output=:s)
    @test sort(collect(keys(both.point_data))) == [:a, :b, :s]
    @test vec(FerriteViz.point_data(both, :s)[]) ≈ 3 .* vec(FerriteViz.point_data(both, :a)[])
    # ... but a genuine geometry change upstream still invalidates them
    stale = ds |> Refine(1)
    set_point_data!(stale, :old, ones(FerriteViz.num_vertices(stale)))
    @test !haskey((stale |> AddQuadraturePointData(qr, vals; output=:a)).point_data, :old)

    # mixed cell types via a per-reference-shape rule mapping
    mnodes = [Node((0.0, 0.0)), Node((1.0, 0.0)), Node((1.0, 1.0)), Node((0.0, 1.0)), Node((2.0, 0.0)), Node((2.0, 1.0))]
    mcells = Ferrite.AbstractCell[Quadrilateral((1, 2, 3, 4)), Triangle((2, 5, 3)), Triangle((5, 6, 3))]
    mds = FEData(DofHandler(Grid(mcells, mnodes)), Float64[])
    qrs = Dict(RefQuadrilateral => qr, RefTriangle => QuadratureRule{RefTriangle}(2))
    mvals = [[Float64(10c + q) for q in 1:getnquadpoints(qrs[Ferrite.getrefshape(mcells[c])])] for c in 1:3]
    @test cellplot(mds |> AddQuadraturePointData(qrs, mvals; output=:iv); color=:iv) isa Makie.FigureAxisPlot

    @test_throws ErrorException ds |> AddQuadraturePointData(qr, vals[1:2]; output=:iv)          # wrong ncells
    @test_throws ErrorException mds |> AddQuadraturePointData(Dict(RefQuadrilateral => qr), mvals; output=:iv) # missing rule
end

@testset "subdomain restrictions error clearly" begin
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid)
    sdh = SubDofHandler(dh, Set(1:2)) # partial coverage
    add!(sdh, :u, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    u = zeros(ndofs(dh))
    ds = FEData(dh, u)
    @test_throws ErrorException ds |> Gradient(:u)
    @test_throws ErrorException FerriteViz.interpolate_gradient_field(dh, u, :u)
end

@testset "filters: crinkle clip, refine" begin
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
end

@testset "derivation filters: Norm1, point-data Deviator, tuple output, Threshold bounds" begin
    grid = generate_grid(Quadrilateral, (3,3))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dh)
    src = FEData(dh, rand(ndofs(dh)))
    nv = FerriteViz.num_vertices(src)
    U = FerriteViz.point_data(src, :u)[]

    # Norm1: rowwise 1-norm
    N1 = FerriteViz.point_data(src |> Norm1(input=:u), :norm1)[]
    @test size(N1, 2) == 1
    @test vec(N1) ≈ [sum(abs, view(U, i, :)) for i in 1:nv]

    # Deviator on point data: tensor-valued output written back componentwise
    devp = src |> Gradient(:u) |> Deviator(input=:gradient)
    G = FerriteViz.point_data(devp, :gradient)[]
    D = FerriteViz.point_data(devp, :deviator)[]
    @test size(D, 2) == 4
    for i in 1:size(D, 1)
        any(isnan, view(G, i, :)) && continue
        @test all(isapprox.(Tuple(view(D, i, :)),
                            FerriteViz._components(Tensors.dev(FerriteViz._wrap_row(view(G, i, :), 2))); atol=1e-12))
    end

    # Derive may return a Tuple -> one column per entry
    set_point_data!(src, :a, collect(1.0:nv))
    T = FerriteViz.point_data(src |> Derive(x -> (x, 2x); input=:a, output=:t), :t)[]
    @test size(T, 2) == 2
    @test T[:, 2] ≈ 2 .* T[:, 1]

    # Threshold on cell data honours both bounds
    ncells = getncells(grid)
    set_cell_data!(src, :cv, collect(1.0:ncells))
    th = FerriteViz.cell_data(src |> Threshold(input=:cv, min=2.0, max=5.0), :threshold)[]
    @test isnan(th[1]) && all(isnan, th[6:end])
    @test th[2:5] == collect(2.0:5.0)

    # the unary derivations work on cell data too
    @test FerriteViz.cell_data(src |> ExtractComponent(1; input=:cv, output=:c1), :c1)[] == collect(1.0:ncells)

    # A row may hold a tensor whose dimension differs from the grid's: shells and
    # plane-strain problems carry 3D stresses on a 2D grid. These used to reach
    # `Tensors.Vec{n}`, which does not exist beyond n = 3, and died with a
    # MethodError from inside Tensors.
    @test FerriteViz._wrap_row(collect(1.0:9), 2) isa Tensors.Tensor{2,3}
    @test FerriteViz._wrap_row(collect(1.0:6), 2) isa Tensors.SymmetricTensor{2,3}
    @test FerriteViz._wrap_row(collect(1.0:6), 3) isa Tensors.SymmetricTensor{2,3}
    @test FerriteViz._wrap_row(collect(1.0:4), 3) isa Tensors.Tensor{2,2}
    # dimension-matched cases still win
    @test FerriteViz._wrap_row(collect(1.0:4), 2) isa Tensors.Tensor{2,2}
    @test FerriteViz._wrap_row(collect(1.0:9), 3) isa Tensors.Tensor{2,3}
    @test FerriteViz._wrap_row(collect(1.0:2), 2) isa Tensors.Vec{2}
    @test FerriteViz._wrap_row(collect(1.0:3), 2) isa Tensors.Vec{3}
    @test FerriteViz._wrap_row([7.0], 3) === 7.0
    # and an uninterpretable width errors clearly instead of throwing from Tensors
    @test_throws ErrorException FerriteViz._wrap_row(collect(1.0:5), 2)
    @test_throws ErrorException FerriteViz._wrap_row(collect(1.0:12), 3)

    # end to end: a 3D stress field stored on a 2D grid reduces with VonMises
    set_point_data!(src, :sigma3d, repeat(reshape(collect(1.0:9), 1, 9), nv, 1))
    vm = FerriteViz.point_data(src |> VonMises(input=:sigma3d), :vonMises)[]
    @test size(vm) == (nv, 1)
    @test all(vm .≈ FerriteViz.vonmises(Tensors.Tensor{2,3}(NTuple{9,Float64}(1.0:9))))

    # Threshold masks values but must not touch the geometry (it is not a clip)
    thp = src |> Threshold(input=:a, min=2.0, max=4.0)
    @test FerriteViz.num_vertices(thp) == nv
    @test FerriteViz.ShaderAbstractions.data(thp.vis_triangles) ==
          FerriteViz.ShaderAbstractions.data(src.vis_triangles)
    @test count(isnan, FerriteViz.point_data(thp, :threshold)[]) == nv - 3
end

@testset "dataset validation and reactivity" begin
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dh)
    ds = FEData(dh, rand(ndofs(dh)))
    nv = FerriteViz.num_vertices(ds)
    ncells = getncells(grid)

    # update! validates against the root solution length
    @test_throws ErrorException FerriteViz.update!(ds, zeros(ndofs(dh) + 1))

    # registered data is validated against the tessellation/grid size
    @test_throws ErrorException set_point_data!(ds, :p, ones(nv + 1))
    @test_throws ErrorException set_cell_data!(ds, :c, ones(ncells + 1))

    # Observable-backed registration stays live
    obs = Makie.Observable(ones(nv))
    set_point_data!(ds, :p, obs)
    @test vec(FerriteViz.point_data(ds, :p)[]) == ones(nv)
    obs[] = fill(2.0, nv)
    @test vec(FerriteViz.point_data(ds, :p)[]) == fill(2.0, nv)
    cobs = Makie.Observable(ones(ncells))
    set_cell_data!(ds, :c, cobs)
    cobs[] = fill(3.0, ncells)
    @test FerriteViz.cell_data(ds, :c)[] == fill(3.0, ncells)

    # unknown names and unreduced data error informatively
    @test_throws ErrorException FerriteViz.point_data(ds, :nope)
    @test_throws ErrorException FerriteViz.cell_data(ds, :nope)
    @test_throws ErrorException FerriteViz._scalar_data(ds, :u) # 2 components, not reduced
    set_cell_data!(ds, :σ, [Tensors.rand(SymmetricTensor{2,2}) for _ in 1:ncells])
    @test_throws ErrorException FerriteViz._scalar_data(ds, :σ) # non-scalar cell data

    # a dof handler without fields has no default field to resolve
    @test_throws ErrorException FerriteViz.point_data(FEData(DofHandler(grid), Float64[]), :default)

    # :default resolves to the first field, and a filter without an explicit
    # input picks the same one
    dhd = DofHandler(grid)
    add!(dhd, :u, Lagrange{RefQuadrilateral,1}()^2)
    add!(dhd, :p, Lagrange{RefQuadrilateral,1}())
    close!(dhd)
    dsd = FEData(dhd, rand(ndofs(dhd)))
    @test FerriteViz._resolve_name(dsd, :default) === :u
    @test FerriteViz._resolve_name(dsd, :p) === :p
    @test size(FerriteViz.point_data(dsd, :default)[], 2) == 2
    let U = FerriteViz.point_data(dsd, :u)[]
        @test vec(FerriteViz.point_data(dsd |> Magnitude(), :magnitude)[]) ≈
              [sqrt(U[i, 1]^2 + U[i, 2]^2) for i in 1:size(U, 1)]
    end

    # `:default` is reserved: a dof field of that name could never be addressed,
    # since naming it resolves to the first field instead. Reject it at
    # construction rather than silently reading the wrong array.
    dhr = DofHandler(grid)
    add!(dhr, :u, Lagrange{RefQuadrilateral,1}()^2)
    add!(dhr, :default, Lagrange{RefQuadrilateral,1}())
    close!(dhr)
    @test_throws ErrorException FEData(dhr, rand(ndofs(dhr)))
    @test_throws ErrorException FEData(dhr, Makie.Observable(rand(ndofs(dhr))))
    # rejected in first position too, so the rule is unconditional
    dhr1 = DofHandler(grid)
    add!(dhr1, :default, Lagrange{RefQuadrilateral,1}())
    close!(dhr1)
    @test_throws ErrorException FEData(dhr1, rand(ndofs(dhr1)))
end

@testset "warp reactivity" begin
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dh)
    ds = FEData(dh, rand(ndofs(dh)))
    nv = FerriteViz.num_vertices(ds)
    set_point_data!(ds, :d1, [ones(nv) zeros(nv)])
    set_point_data!(ds, :d2, [zeros(nv) ones(nv)])

    # an Observable scale streams into the coordinates
    scale = Makie.Observable(1.0)
    w = ds |> WarpByVector(:d1, scale)
    @test w.coords[] ≈ [c .+ Float32.((1, 0)) for c in ds.coords[]]
    scale[] = 3.0
    @test w.coords[] ≈ [c .+ Float32.((3, 0)) for c in ds.coords[]]

    # an Observable field name rewires the displacement source
    fname = Makie.Observable(:d1)
    wf = ds |> WarpByVector(fname, 1.0)
    @test wf.coords[] ≈ [c .+ Float32.((1, 0)) for c in ds.coords[]]
    fname[] = :d2
    @test wf.coords[] ≈ [c .+ Float32.((0, 1)) for c in ds.coords[]]

    # the deformation field's component count is validated
    dhs = DofHandler(grid); add!(dhs, :t, Lagrange{RefQuadrilateral,1}()); close!(dhs)
    @test_throws ErrorException FEData(dhs, rand(ndofs(dhs))) |> WarpByVector(:t)
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
    # composable viewer (SpecApi): defaults reproduce a solutionplot panel
    @test ferriteviewer(ds) isa Makie.Figure
    @test ferriteviewer(ds, [u, 2u]) isa Makie.Figure
    # viewer on a scalar-only dataset (no deformable field)
    dhs = DofHandler(grid); add!(dhs, :t, Lagrange{RefQuadrilateral,1}()); close!(dhs)
    @test ferriteviewer(FEData(dhs, rand(ndofs(dhs)))) isa Makie.Figure
    # spec building blocks are composable and typed
    @test FerriteViz.solutionplotspec(ds; color=:default) isa Makie.PlotSpec
    @test FerriteViz.panelspec(FerriteViz.solutionplotspec(ds; color=:default); dim=2) isa Makie.GridLayoutSpec
    @test FerriteViz.default_layout(ds, (field=:u, process="magnitude", colormap=:inferno, labels=false)) isa Makie.GridLayoutSpec
    # regression: realizing a panel whose Colorbar links to a spec WITHOUT an
    # explicit colormap must not trip Makie's lookup_default on our recipes,
    # and the bar must take its range from the named data instead of (0, 1)
    let sol = FerriteViz.solutionplotspec(ds; color=:default)
        @test Makie.plot(FerriteViz.panelspec(sol; colorbar=sol, dim=2)) isa Union{Makie.FigureAxisPlot,Makie.Figure}
        vals = FerriteViz._scalar_data(ds, :default; reduce_default=true)[]
        @test FerriteViz._spec_colorrange(sol) == (minimum(vals), maximum(vals))
        # an explicit colorrange wins (PlotSpec stores it converted), and a
        # plain colour yields no range at all
        explicit = FerriteViz.solutionplotspec(ds; color=:default, colorrange=(0, 2))
        @test FerriteViz._colorbar_kw(explicit, (;))[:colorrange] == explicit.kwargs[:colorrange]
        @test FerriteViz._spec_colorrange(FerriteViz.solutionplotspec(ds; color=:red)) === nothing
    end
    # custom layout + custom pluggable control, with a reactive structural update
    flagobs = Ref{Any}(nothing)
    probe = FerriteViz.Control() do fig, d
        tog = Makie.Toggle(fig)
        flagobs[] = tog.active
        FerriteViz.ControlResult(Any[tog]; structural=[:flag => tog.active])
    end
    mylayout(d, s) = FerriteViz.panelspec(
        FerriteViz.solutionplotspec(d; color=:default, colormap=(s.flag ? :inferno : :viridis)); dim=2)
    vfig = ferriteviewer(ds; layout=mylayout, controls=FerriteViz.Control[probe])
    @test vfig isa Makie.Figure
    flagobs[][] = true                 # flip control -> re-diff spec through the whole pipeline
    @test vfig isa Makie.Figure        # (no error thrown by the spec update)
    # 3D wedge/pyramid render out of the box
    wgrid = Grid([Wedge((1,2,3,4,5,6))], [Node(Vec(0.0,0.0,0.0)), Node(Vec(1.0,0.0,0.0)), Node(Vec(0.0,1.0,0.0)),
                                          Node(Vec(0.0,0.0,1.0)), Node(Vec(1.0,0.0,1.0)), Node(Vec(0.0,1.0,1.0))])
    wdh = DofHandler(wgrid); add!(wdh, :u, Lagrange{RefPrism,1}()); close!(wdh)
    @test solutionplot(wdh, rand(ndofs(wdh))) isa Makie.FigureAxisPlot
    # 1D grids: meshplot pads coordinates to 2D
    lgrid = generate_grid(Line, (3,))
    @test meshplot(lgrid; nodelabels=true, celllabels=true) isa Makie.FigureAxisPlot
end

@testset "viewer building blocks" begin
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dh)
    u = rand(ndofs(dh))
    ds = FEData(dh, u)

    # every spec helper returns a PlotSpec
    for spec in (meshplotspec(ds), cellplotspec(ds), surfaceplotspec(ds), arrowplotspec(ds))
        @test spec isa Makie.PlotSpec
    end

    @test_throws ErrorException FerriteViz.ControlResult(Any[]; placement=:left)

    # the default pipeline warps by the observable deformation scale ...
    scale = Makie.Observable(0.0)
    p = default_pipeline(ds, Dict{Symbol,Makie.Observable}(:deform_scale => scale))
    @test p.coords[] ≈ ds.coords[]
    scale[] = 1.0
    @test !(p.coords[] ≈ ds.coords[])
    # ... and passes a dataset without deformable field through untouched
    dhs = DofHandler(grid); add!(dhs, :t, Lagrange{RefQuadrilateral,1}()); close!(dhs)
    dss = FEData(dhs, rand(ndofs(dhs)))
    @test default_pipeline(dss, Dict{Symbol,Makie.Observable}(:deform_scale => Makie.Observable(1.0))) === dss
    @test length(default_controls(ds)) == length(default_controls(dss)) + 1 # DeformationToggle

    # process reductions register derived point-data arrays under stable names
    U = FerriteViz.point_data(ds, :u)[]
    n1 = FerriteViz._viewer_color_array!(ds, :u, "x₁")
    @test n1 === Symbol("u_x₁")
    @test vec(FerriteViz.point_data(ds, n1)[]) ≈ U[:, 1]
    nm = FerriteViz._viewer_color_array!(ds, :u, "magnitude")
    @test vec(FerriteViz.point_data(ds, nm)[]) ≈ [norm(view(U, i, :)) for i in 1:FerriteViz.num_vertices(ds)]
    @test FerriteViz._viewer_color_array!(dss, :t, "magnitude") === :t # scalar short-circuit
    @test FerriteViz.default_layout(ds, (field=:u, process="x₂", colormap=:viridis, labels=true, wireframe=false)) isa Makie.GridLayoutSpec

    # the TimeSlider streams solutions through update!
    fig = Makie.Figure()
    r = TimeSlider([u, 2u]).make(fig, ds)
    @test r.placement === :below
    Makie.set_close_to!(r.content[1].sliders[1], 2)
    @test ds.u[] ≈ 2 .* u

    # a 3D dataset builds its panel on an LScene
    grid3 = generate_grid(Hexahedron, (2,2,2))
    dh3 = DofHandler(grid3); add!(dh3, :u, Lagrange{RefHexahedron,1}()^3); close!(dh3)
    @test ferriteviewer(FEData(dh3, rand(ndofs(dh3)))) isa Makie.Figure
end

@testset "representation color resolution" begin
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dh)
    ds = FEData(dh, rand(ndofs(dh)))
    nv = FerriteViz.num_vertices(ds)
    set_point_data!(ds, :s1, collect(1.0:nv))
    set_point_data!(ds, :s2, 2 .* collect(1.0:nv))

    fig, ax, plt = solutionplot(ds; color=:s1)
    mesh = plt.plots[1]
    @test mesh.color[] ≈ collect(1.0:nv)
    plt.color = :s2                          # switch to another named array
    @test mesh.color[] ≈ 2 .* collect(1.0:nv)
    plt.color = :s1                          # and back to the first one
    @test mesh.color[] ≈ collect(1.0:nv)
    ds.point_data[:s1][] = fill(7.0, nv, 1)  # updating the array flows into the plot
    @test mesh.color[] ≈ fill(7.0, nv)
    # a plain (non-data-name) color passes through at construction; switching a
    # live plot from a data array to a plain color is not supported, since the
    # mesh's color input is typed by its initial value
    figc, axc, pltc = solutionplot(ds; color=:red)
    @test !(pltc.plots[1].color[] isa AbstractVector)

    # arrowplot accepts named scalar arrays and plain colors ...
    @test arrowplot(ds; color=:s1) isa Makie.FigureAxisPlot
    @test arrowplot(ds; color=:orange) isa Makie.FigureAxisPlot
    # ... but the arrow field itself must be vector-valued
    dhs = DofHandler(grid); add!(dhs, :t, Lagrange{RefQuadrilateral,1}()); close!(dhs)
    @test_throws ErrorException arrowplot(FEData(dhs, rand(ndofs(dhs))))
    # surfaceplot needs a data array, not a plain color
    @test_throws ErrorException surfaceplot(ds; color=:red)

    # CairoMakie shim units: Buffers unwrap to their vectors, plain data passes through
    @test FerriteViz._buffer_data(ds.coords_buffer) isa Vector
    @test FerriteViz._buffer_data([1, 2]) == [1, 2]
end

# ---- exact cutting: shared helpers -----------------------------------------
_vv3(p) = Vec{3,Float64}(NTuple{3,Float64}(p))
_vv2(p) = Vec{2,Float64}(NTuple{2,Float64}(p))
_tetv(a, b, c, d) = abs(FerriteViz._signed_tet_volume(a, b, c, d))
# signed on purpose: the simplices invariant is positive orientation, so any
# negative tet would make the analytic volume checks fail. Only solid cells
# count — CrinkleClip removes cells by mask without rebuilding the geometry.
function _kept_volume(ds::FerriteViz.FEData{3})
    pc = _vv3.(ds.coords[])
    return sum(FerriteViz._signed_tet_volume(pc[s[1]], pc[s[2]], pc[s[3]], pc[s[4]])
               for (i, s) in enumerate(ds.simplices) if ds.solid[ds.simplex_cell_map[i]]; init=0.0)
end
_all_positive(ds::FerriteViz.FEData{3}) = begin
    pc = _vv3.(ds.coords[])
    all(FerriteViz._signed_tet_volume(pc[s[1]], pc[s[2]], pc[s[3]], pc[s[4]]) > 0 for s in ds.simplices)
end
function _cap_area(ds, n, d; tol=1e-6)
    pc = _vv3.(ds.coords[])
    a = 0.0
    for t in ds.all_triangles
        i, j, k = convert(Int, t[1]), convert(Int, t[2]), convert(Int, t[3])
        if all(abs(pc[v] ⋅ n - d) <= tol for v in (i, j, k))
            a += FerriteViz._tri_area(pc[i], pc[j], pc[k])
        end
    end
    return a
end
function _surface_area(ds::FerriteViz.FEData{3})
    pc = _vv3.(ds.coords[])
    return sum(FerriteViz._tri_area(pc[convert(Int, t[1])], pc[convert(Int, t[2])], pc[convert(Int, t[3])])
               for t in ds.all_triangles; init=0.0)
end
function _line_length(ds::FerriteViz.FEData{2})
    pc = _vv2.(ds.coords[])
    return sum(norm(pc[e[2]] - pc[e[1]]) for e in ds.all_edges; init=0.0)
end
function _hexds(nel; ip=Lagrange{RefHexahedron,1}(), f=nothing)
    grid = generate_grid(Hexahedron, (nel, nel, nel))
    dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
    u = zeros(ndofs(dh))
    f !== nothing && Ferrite.apply_analytical!(u, dh, :u, f)
    return FEData(dh, u), dh, u
end

@testset "volume simplices: infrastructure" begin
    # fan tets tile every reference shape exactly
    for (CT, ip, nel, vol) in [(Tetrahedron, Lagrange{RefTetrahedron,1}(), 2, 8.0),
                               (Hexahedron, Lagrange{RefHexahedron,1}(), 2, 8.0)]
        grid = generate_grid(CT, ntuple(_ -> nel, 3))
        dh = DofHandler(grid); add!(dh, :u, ip); close!(dh)
        ds = FEData(dh, zeros(ndofs(dh)))
        @test _kept_volume(ds) ≈ vol atol = 1e-9
        @test _all_positive(ds)
        @test all(ds.solid) && ds.cells_intact
        @test length(ds.simplices) == length(ds.all_triangles)
    end
    wgrid = Grid([Wedge((1, 2, 3, 4, 5, 6))],
                 [Node(Vec(0.0, 0.0, 0.0)), Node(Vec(1.0, 0.0, 0.0)), Node(Vec(0.0, 1.0, 0.0)),
                  Node(Vec(0.0, 0.0, 1.0)), Node(Vec(1.0, 0.0, 1.0)), Node(Vec(0.0, 1.0, 1.0))])
    wdh = DofHandler(wgrid); add!(wdh, :u, Lagrange{RefPrism,1}()); close!(wdh)
    @test _kept_volume(FEData(wdh, zeros(ndofs(wdh)))) ≈ 0.5 atol = 1e-9
    pgrid = Grid([Pyramid((1, 2, 3, 4, 5))],
                 [Node(Vec(0.0, 0.0, 0.0)), Node(Vec(1.0, 0.0, 0.0)), Node(Vec(0.0, 1.0, 0.0)),
                  Node(Vec(1.0, 1.0, 0.0)), Node(Vec(0.5, 0.5, 1.0))])
    pdh = DofHandler(pgrid); add!(pdh, :u, Lagrange{RefPyramid,1}()); close!(pdh)
    @test _kept_volume(FEData(pdh, zeros(ndofs(pdh)))) ≈ 1 / 3 atol = 1e-9

    # 2D: simplices alias the triangles, no extra vertices
    grid2 = generate_grid(Quadrilateral, (3, 3))
    dh2 = DofHandler(grid2); add!(dh2, :u, Lagrange{RefQuadrilateral,1}()); close!(dh2)
    ds2 = FEData(dh2, zeros(ndofs(dh2)))
    @test length(ds2.simplices) == length(ds2.all_triangles)
    @test FerriteViz.num_vertices(ds2) == 5 * getncells(grid2)

    # transfer is solid-gated: hidden interior cells carry real values
    f(x) = x[1] + 2x[2] - x[3]
    ds3, _, _ = _hexds(3; f)
    @test !all(ds3.visible)
    data = FerriteViz.point_data(ds3, :u)[]
    @test !any(isnan, data)
    @test all(isapprox.(vec(data), [f(_vv3(x)) for x in ds3.coords[]]; atol=1e-4))

    # CrinkleClip maintains solid; registered arrays survive, dof caches drop
    ds4, dh4, _ = _hexds(3)
    set_point_data!(ds4, :reg, ones(FerriteViz.num_vertices(ds4)))
    FerriteViz.point_data(ds4, :u)
    c4 = ds4 |> CrinkleClip(ClipPlane(Vec((1.0, 0.0, 0.0)), 0.0))
    @test sum(c4.solid) < getncells(Ferrite.get_grid(dh4)) && all(ds4.solid)
    @test haskey(c4.point_data, :reg) && !haskey(c4.point_data, :u)
    c5 = c4 |> CrinkleClip(ClipPlane(Vec((0.0, 1.0, 0.0)), 0.0))
    @test all(c4.solid[c5.solid])

    # refine rebuilds simplices and preserves the volume
    ds6, _, _ = _hexds(2)
    r6 = ds6 |> Refine(1)
    @test length(r6.simplices) == 4 * length(ds6.simplices)
    @test _kept_volume(r6) ≈ 8.0 atol = 1e-6
    @test _all_positive(r6)

    # embedded 2D cells in a 3D grid: no volume simplices, Clip cuts the
    # triangles without fabricating caps
    snodes = [Node(Vec(0.0, 0.0, 0.0)), Node(Vec(1.0, 0.0, 0.3)), Node(Vec(1.0, 1.0, 0.5)),
              Node(Vec(0.0, 1.0, 0.2)), Node(Vec(2.0, 0.0, 0.6)), Node(Vec(2.0, 1.0, 0.8))]
    sgrid = Grid([Quadrilateral((1, 2, 3, 4)), Quadrilateral((2, 5, 6, 3))], snodes)
    sdh = DofHandler(sgrid); add!(sdh, :u, Lagrange{RefQuadrilateral,1}()); close!(sdh)
    sds = FEData(sdh, rand(ndofs(sdh)))
    @test isempty(sds.simplices)
    scl = sds |> Clip(ClipPlane(Vec(1.0, 0.0, 0.0), 0.7))
    @test isempty(scl.simplices) && !isempty(scl.all_triangles)
    @test all(p[1] <= 0.7 + 1e-5 for p in scl.coords[])
end

@testset "cutting core: conservation and side rules" begin
    sv(p, t) = FerriteViz._signed_tet_volume(p.pos[t[1]], p.pos[t[2]], p.pos[t[3]], p.pos[t[4]])
    ta(p, t) = FerriteViz._tri_area(p.pos[t[1]], p.pos[t[2]], p.pos[t[3]])
    for trial in 1:50
        pts = [GeometryBasics.Point{3,Float32}(rand(3)...) for _ in 1:4]
        pf = [_vv3(x) for x in pts]
        vol = _tetv(pf...)
        vol < 1e-4 && continue
        n = _vv3(rand(3) .- 0.5); n /= norm(n)
        d = 0.2 + 0.6 * rand()
        parts = Float64[]
        caps = Float64[]
        for s in (1.0, -1.0)
            pool = FerriteViz.CutVertexPool(pts)
            g = FerriteViz.snap!([s * (pf[i] ⋅ n - d) for i in 1:4], 1e-12)
            tets = NTuple{4,Int}[]; tc = Int[]; cp = NTuple{3,Int}[]; cc = Int[]
            FerriteViz.clip_tet!(pool, tets, tc, cp, cc, (1, 2, 3, 4), 1, g, 1e-14, 1e-14)
            push!(parts, sum(abs(sv(pool, t)) for t in tets; init=0.0))
            push!(caps, sum(ta(pool, t) for t in cp; init=0.0))
            for t in cp, v in t   # cap vertices on the plane
                @test abs(pool.pos[v] ⋅ n - d) <= 1e-9
            end
        end
        @test sum(parts) ≈ vol atol = 1e-10
        @test caps[1] ≈ caps[2] atol = 1e-10   # both sides share the cross-section
    end
    # exhaustive sign patterns: the emitted cap triangles are exactly the kept
    # tets' boundary faces on the cut plane (consistent quad diagonals), and
    # kept tets are positively oriented
    let pts = [GeometryBasics.Point{3,Float32}(0, 0, 0), GeometryBasics.Point{3,Float32}(1, 0, 0),
               GeometryBasics.Point{3,Float32}(0, 1, 0), GeometryBasics.Point{3,Float32}(0, 0, 1)]
        for pattern in Iterators.product((-1.0, 0.0, 1.0), (-1.0, 0.0, 1.0), (-1.0, 0.0, 1.0), (-1.0, 0.0, 1.0))
            g = FerriteViz.snap!(collect(pattern) .* 0.37, 1e-12)
            pool = FerriteViz.CutVertexPool(pts)
            tets = NTuple{4,Int}[]; tc = Int[]; cp = NTuple{3,Int}[]; cc = Int[]
            FerriteViz.clip_tet!(pool, tets, tc, cp, cc, (1, 2, 3, 4), 1, g, 1e-14, 1e-14)
            @test all(sv(pool, t) > 0 for t in tets)
            gout = vec(FerriteViz.combine_rows(pool.combos, reshape(g, :, 1)))
            onplane(tri) = all(abs(gout[v]) <= 1e-12 for v in tri)
            # compare by positions: t=0/1 cuts create coincident duplicates of
            # kept vertices (distinct indices, the documented conformity
            # contract), so index equality would be too strict
            facekey(f) = Tuple(sort!([round.(Float64.(Tuple(pool.pos[v])); digits=9) for v in f]))
            tetfaces = Set()
            for t in tets, f in ((t[1], t[2], t[3]), (t[1], t[2], t[4]), (t[1], t[3], t[4]), (t[2], t[3], t[4]))
                (onplane(f) && ta(pool, f) > 1e-14) && push!(tetfaces, facekey(f))
            end
            capset = Set(facekey(f) for f in cp)
            if any(<(0), g) && any(>(0), g)   # genuinely cut
                @test capset == tetfaces
            elseif !any(>(0), g)              # nothing outside: kept whole, no caps
                @test isempty(capset)
            else                              # nothing strictly inside: no kept volume
                # (the Clip filter removes such cells before reaching the core)
                @test isempty(tets)
            end
        end
    end
    # a face exactly at the level is emitted exactly once, and not at all
    # inside a plateau
    pts = [GeometryBasics.Point{3,Float32}(0, 0, 0), GeometryBasics.Point{3,Float32}(1, 0, 0),
           GeometryBasics.Point{3,Float32}(0, 1, 0), GeometryBasics.Point{3,Float32}(0, 0, 1),
           GeometryBasics.Point{3,Float32}(0, 0, -1)]
    pool = FerriteViz.CutVertexPool(pts)
    g = FerriteViz.snap!([0.0, 0.0, 0.0, 1.0, -1.0], 1e-12)
    tris = NTuple{3,Int}[]; cells = Int[]
    FerriteViz.march_tet!(pool, tris, cells, (1, 2, 3, 4), 1, g, 1e-14, 1)
    FerriteViz.march_tet!(pool, tris, cells, (1, 2, 3, 5), 1, g, 1e-14, 1)
    @test length(tris) == 1
    pool2 = FerriteViz.CutVertexPool(pts)
    g2 = FerriteViz.snap!([0.0, 0.0, 0.0, -1.0, -1.0], 1e-12)
    tris2 = NTuple{3,Int}[]; cells2 = Int[]
    FerriteViz.march_tet!(pool2, tris2, cells2, (1, 2, 3, 4), 1, g2, 1e-14, 1)
    FerriteViz.march_tet!(pool2, tris2, cells2, (1, 2, 3, 5), 1, g2, 1e-14, 1)
    @test isempty(tris2)
end

@testset "Clip: geometry and corner cases" begin
    f(x) = x[1] + 2x[2] - x[3]
    ds, dh, _ = _hexds(3; f)
    grid = Ferrite.get_grid(dh)
    n = Vec(1.0, 0.0, 0.0)
    clipped = ds |> Clip(ClipPlane(n, 0.1))
    @test _kept_volume(clipped) ≈ 1.1 * 4 atol = 1e-6
    @test _cap_area(clipped, n, 0.1) ≈ 4.0 atol = 1e-6
    @test !clipped.cells_intact
    pc = _vv3.(clipped.coords[])
    @test all(p ⋅ n <= 0.1 + 1e-5 for p in pc)
    # dof data is exact at the cut positions
    data = FerriteViz.point_data(clipped, :u)[]
    @test all(isapprox(data[i], f(pc[i]); atol=1e-4) for i in eachindex(pc))
    # cap triangles stay connected to the cell they cut through
    for (t, tri) in enumerate(clipped.all_triangles)
        i, j, k = convert(Int, tri[1]), convert(Int, tri[2]), convert(Int, tri[3])
        all(abs(pc[v][1] - 0.1) < 1e-6 for v in (i, j, k)) || continue
        xs = [x[1] for x in Ferrite.getcoordinates(grid, clipped.triangle_cell_map[t])]
        @test minimum(xs) - 1e-6 <= sum(pc[v][1] for v in (i, j, k)) / 3 <= maximum(xs) + 1e-6
    end
    # complementary clips partition the volume
    @test _kept_volume(clipped) + _kept_volume(ds |> Clip(ClipPlane(-n, -0.1))) ≈ 8.0 atol = 1e-6
    # the previously hidden interior cell is revealed by the cut
    interior = findfirst(i -> !ds.visible[i], 1:getncells(grid))
    @test clipped.visible[interior]
    # the meshplot wireframe is clipped along with the surface
    @test !isempty(clipped.all_edges)
    @test all(max(pc[e[1]][1], pc[e[2]][1]) <= 0.1 + 1e-5 for e in clipped.all_edges)
    # oblique plane: half volume, hexagonal cross-section
    no = Vec(1.0, 1.0, 1.0) / sqrt(3.0)
    ob = ds |> Clip(ClipPlane(no, 0.0))
    @test _kept_volume(ob) ≈ 4.0 atol = 1e-6
    @test _cap_area(ob, no, 0.0) ≈ 3 * sqrt(3.0) atol = 1e-5

    # node-aligned plane: whole-cell removal, cells stay intact, faces exposed
    aligned = ds |> Clip(ClipPlane(n, 1 / 3))
    @test aligned.cells_intact && sum(aligned.solid) == 18
    @test _kept_volume(aligned) ≈ (4 / 3) * 4 atol = 1e-6
    @test _cap_area(aligned, n, 1 / 3) ≈ 4.0 atol = 1e-6
    # plane outside the bbox: identity / empty
    @test _kept_volume(ds |> Clip(ClipPlane(n, 5.0))) ≈ 8.0 atol = 1e-6
    emptied = ds |> Clip(ClipPlane(n, -5.0))
    @test isempty(emptied.all_triangles) && !any(emptied.solid)
    # errors
    @test_throws ErrorException FEData(DofHandler(generate_grid(Quadrilateral, (2, 2))), Float64[]) |> Clip(ClipPlane(n, 0.0))
    @test_throws ErrorException ds |> Clip(ClipPlane(Vec(0.0, 0.0, 0.0), 0.0))
    @test_throws ErrorException ds |> Clip(ClipPlane(n, Inf))
    # translated Float32 geometry
    tgrid0 = generate_grid(Hexahedron, (2, 2, 2))
    tg = Grid(tgrid0.cells, [Node(nd.x + Vec(1000.0, 2000.0, -500.0)) for nd in tgrid0.nodes])
    tdh = DofHandler(tg); add!(tdh, :u, Lagrange{RefHexahedron,1}()); close!(tdh)
    tclip = FEData(tdh, zeros(ndofs(tdh))) |> Clip(ClipPlane(n, 1000.1))
    @test _kept_volume(tclip) ≈ 1.1 * 4 rtol = 1e-3
    # scale invariance: uniformly tiny and huge grids clip identically
    for s in (1e-6, 1e6)
        sg = Grid(tgrid0.cells, [Node(nd.x * s) for nd in tgrid0.nodes])
        sdh = DofHandler(sg); add!(sdh, :u, Lagrange{RefHexahedron,1}()); close!(sdh)
        sclip = FEData(sdh, zeros(ndofs(sdh))) |> Clip(ClipPlane(n, 0.1 * s))
        @test _kept_volume(sclip) ≈ 1.1 * 4 * s^3 rtol = 1e-3
    end
    # graded mesh: a cell far smaller than the domain is cut, not discarded
    # (degeneracy pruning is cell-local; only classification is global)
    gnodes = [Node(Vec(x, y, z)) for x in (0.0, 1.0) for y in (0.0, 1.0) for z in (0.0, 1.0)]
    tiny(v) = Vec(2.0, 0.0, 0.0) + 1e-4 * v
    append!(gnodes, [Node(tiny(Vec(x, y, z))) for x in (0.0, 1.0) for y in (0.0, 1.0) for z in (0.0, 1.0)])
    order = (1, 5, 7, 3, 2, 6, 8, 4)  # (x,y,z) loop order -> hex node order
    gcells = [Hexahedron(ntuple(i -> order[i], 8)), Hexahedron(ntuple(i -> order[i] + 8, 8))]
    ggrid = Grid(gcells, gnodes)
    gdh = DofHandler(ggrid); add!(gdh, :u, Lagrange{RefHexahedron,1}()); close!(gdh)
    gclip = FEData(gdh, zeros(ndofs(gdh))) |> Clip(ClipPlane(n, 2.0 + 5e-5))
    @test gclip.solid[2]
    pcg = _vv3.(gclip.coords[])
    tinyvol = sum(FerriteViz._signed_tet_volume(pcg[s[1]], pcg[s[2]], pcg[s[3]], pcg[s[4]])
                  for s in gclip.simplices[FerriteViz.simplices_on_cell(gclip, 2)]; init=0.0)
    @test tinyvol ≈ 0.5e-4 * 1e-8 rtol = 1e-2
end

@testset "Clip: composition" begin
    fv(x) = Vec{3}((0.1x[1], 0.05x[2], -0.1x[3]))
    grid = generate_grid(Hexahedron, (3, 3, 3))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefHexahedron,1}()^3); close!(dh)
    u = zeros(ndofs(dh)); Ferrite.apply_analytical!(u, dh, :u, fv)
    ds = FEData(dh, u)
    n = Vec(1.0, 0.0, 0.0)

    # clip of clip: slab, and oblique clips commute
    slab = ds |> Clip(ClipPlane(n, 0.5)) |> Clip(ClipPlane(-n, 0.25))
    @test _kept_volume(slab) ≈ 0.75 * 4 atol = 1e-6
    p1 = ClipPlane(Vec(1.0, 1.0, 0.0) / sqrt(2.0), 0.2)
    p2 = ClipPlane(Vec(0.3, -1.0, 0.9) / norm(Vec(0.3, -1.0, 0.9)), -0.1)
    p3 = ClipPlane(Vec(0.0, 0.0, 1.0), 0.2)
    @test _kept_volume(ds |> Clip(p1) |> Clip(p2) |> Clip(p3)) ≈
          _kept_volume(ds |> Clip(p3) |> Clip(p2) |> Clip(p1)) atol = 1e-6
    # three successive clips with an analytic volume, the third cutting
    # through the earlier caps: {x<=0.5, y<=0.25, x+y<=0} of the (-1,1)^3 box
    tri3 = ds |> Clip(ClipPlane(n, 0.5)) |> Clip(ClipPlane(Vec(0.0, 1.0, 0.0), 0.25)) |>
                 Clip(ClipPlane(Vec(1.0, 1.0, 0.0) / sqrt(2.0), 0.0))
    @test _kept_volume(tri3) ≈ 2 * (1.5 * 1.25 - 0.75^2 / 2) atol = 1e-6
    @test _all_positive(tri3)

    # crinkle∘clip and clip∘crinkle: removed cells stay removed
    ck = ds |> CrinkleClip(ClipPlane(Vec(0.0, 1.0, 0.0), 0.0))
    ckc = ck |> Clip(ClipPlane(n, 0.1))
    @test !any(ckc.solid .& .!ck.solid)
    @test _kept_volume(ckc) ≈ 1.1 * (2 / 3) * 2 atol = 1e-5
    cck = ds |> Clip(ClipPlane(n, 0.1)) |> CrinkleClip(ClipPlane(Vec(0.0, 1.0, 0.0), 0.0))
    @test sum(cck.solid) < sum((ds |> Clip(ClipPlane(n, 0.1))).solid)
    @test _kept_volume(cck) ≈ 1.1 * (2 / 3) * 2 atol = 1e-5

    # warp→clip cuts the warped geometry (interior cells included), clip→warp
    # evaluates the displacement at the cut vertices
    w = ds |> WarpByVector(:u, 1.0)
    wc = w |> Clip(ClipPlane(n, 0.0))
    @test all(p[1] <= 1e-4 for p in wc.coords[])
    @test _kept_volume(wc) ≈ 0.5 * (2 * 1.1) * (2 * 1.05) * (2 * 0.9) atol = 1e-4
    cw = ds |> Clip(ClipPlane(n, 0.0)) |> WarpByVector(:u, 1.0)
    @test all(p[1] <= 1e-4 for p in cw.coords[])

    # downstream pipeline
    pipe = ds |> Clip(ClipPlane(n, 0.1)) |> Gradient(:u) |> VonMises(input=:gradient)
    @test count(isfinite, FerriteViz.point_data(pipe, :vonMises)[]) > 0
    # a refined dataset clips (resolve curvature *before* cutting)
    rc = ds |> Refine(1) |> Clip(ClipPlane(n, 0.1))
    @test _kept_volume(rc) ≈ 1.1 * 4 atol = 1e-5
    @test _cap_area(rc, n, 0.1) ≈ 4.0 atol = 1e-5

    # geometry-rebuilding filters reject cut inputs
    qr = QuadratureRule{RefHexahedron}(2)
    vals = [rand(getnquadpoints(qr)) for _ in 1:getncells(grid)]
    clipped = ds |> Clip(ClipPlane(n, 0.1))
    @test_throws ErrorException clipped |> AddQuadraturePointData(qr, vals)
    @test_throws ErrorException clipped |> Refine(1)
    dh1 = DofHandler(grid); add!(dh1, :u, Lagrange{RefHexahedron,1}()); close!(dh1)
    ds1 = FEData(dh1, rand(ndofs(dh1)))
    # ... while whole-cell removal keeps them working
    aligned = ds1 |> Clip(ClipPlane(n, 1 / 3))
    @test (aligned |> AddQuadraturePointData(qr, vals; output=:iv)) isa FerriteViz.FEData
    @test (aligned |> Refine(1)) isa FerriteViz.FEData

    # cache-order independence
    dsA = FEData(dh1, ds1.u[]); dsB = FEData(dh1, ds1.u[])
    FerriteViz.point_data(dsA, :u)
    @test FerriteViz.point_data(dsA |> Clip(ClipPlane(n, 0.3)), :u)[] ≈
          FerriteViz.point_data(dsB |> Clip(ClipPlane(n, 0.3)), :u)[]
end

@testset "Clip: data rules and reactivity" begin
    f(x) = x[1] + 2x[2] - x[3]
    ds, dh, _ = _hexds(2; f)
    n3 = Vec(0.0, 0.0, 1.0)
    # registered point data lerps exactly for data linear in the coordinates
    lin = [p ⋅ Vec(1.0, -1.0, 0.5) for p in _vv3.(ds.coords[])]
    set_point_data!(ds, :lin, lin)
    set_cell_data!(ds, :cd, collect(1.0:getncells(Ferrite.get_grid(dh))))
    clipped = ds |> Clip(ClipPlane(n3, 0.3))
    L = vec(FerriteViz.point_data(clipped, :lin)[])
    pcl = _vv3.(clipped.coords[])
    @test all(isapprox(L[i], pcl[i] ⋅ Vec(1.0, -1.0, 0.5); atol=1e-4) for i in eachindex(pcl))
    @test FerriteViz.cell_data(clipped, :cd)[] == collect(1.0:getncells(Ferrite.get_grid(dh)))
    # Threshold NaNs pass through the lerp
    th = ds |> Threshold(input=:lin, min=0.0) |> Clip(ClipPlane(n3, 0.3))
    @test any(isnan, FerriteViz.point_data(th, :threshold)[])
    @test count(isfinite, FerriteViz.point_data(th, :threshold)[]) > 0

    # frozen topology, reactive positions and colors: with the pure-x warp
    # x' = (1+0.1s)x every output coordinate is an affine combination of the
    # parents, so doubling the solution rescales x by exactly (1.2/1.1)
    grid = generate_grid(Hexahedron, (2, 2, 2))
    dhr = DofHandler(grid); add!(dhr, :u, Lagrange{RefHexahedron,1}()^3); close!(dhr)
    fw(x) = Vec{3}((0.1x[1], 0.0, 0.0))
    ur = zeros(ndofs(dhr)); Ferrite.apply_analytical!(ur, dhr, :u, fw)
    dsr = FEData(dhr, ur)
    w = dsr |> WarpByVector(:u, 1.0)
    c = w |> Clip(ClipPlane(Vec(1.0, 0.0, 0.0), 0.1))
    @test !c.cells_intact
    cb = copy(c.coords[])
    colors = copy(FerriteViz.point_data(c, :u)[])
    FerriteViz.update!(dsr, 2 .* ur)
    scalex(p) = typeof(p)(p[1] * (1.2f0 / 1.1f0), p[2], p[3])
    @test all(isapprox(c.coords[][i], scalex(cb[i]); atol=1e-5) for i in eachindex(cb))
    @test !(FerriteViz.point_data(c, :u)[] ≈ colors)
    # same exactness for extracted isosurfaces under upstream deformation
    dsr2 = FEData(dhr, ur)
    iso = dsr2 |> WarpByVector(:u, 1.0) |> ExtractComponent(1; input=:u, output=:u1) |>
          ExtractIsosurfaces(0.05; input=:u1)
    ib = copy(iso.coords[])
    @test !isempty(ib)
    FerriteViz.update!(dsr2, 2 .* ur)
    @test all(isapprox(iso.coords[][i], scalex(ib[i]); atol=1e-5) for i in eachindex(ib))

    # smoke; cut datasets fall back to the static tessellation (the adaptive
    # base refines whole cells)
    dss, _, _ = _hexds(2; f=x -> x[1]^2 + x[2])
    cs = dss |> Clip(ClipPlane(Vec(0.0, 1.0, 1.0) / sqrt(2.0), 0.1))
    @test FerriteViz._adaptive_capable(dss)
    @test !FerriteViz._adaptive_capable(cs)
    @test solutionplot(cs) isa Makie.FigureAxisPlot
    @test meshplot(cs) isa Makie.FigureAxisPlot
    @test cellplot(cs, rand(8)) isa Makie.FigureAxisPlot
end

@testset "QP Voronoi volume and Clip on quadrature-point data" begin
    _tessvol(t) = sum(_tetv(t.coords[s[1]], t.coords[s[2]], t.coords[s[3]], t.coords[s[4]]) for s in t.simplices; init=0.0)
    for RS in (RefTetrahedron, RefHexahedron, RefPrism, RefPyramid)
        t = FerriteViz.qp_voronoi_tessellation(RS, QuadratureRule{RS}(2))
        surf = FerriteViz.reference_tessellation(RS)
        ctr = sum(surf.coords) / length(surf.coords)
        refvol = sum(_tetv(surf.coords[tri[1]], surf.coords[tri[2]], surf.coords[tri[3]], ctr) for tri in surf.triangles)
        @test _tessvol(t) ≈ refvol atol = 1e-9
        for s in t.simplices   # no tet spans regions
            @test t.vertex_qp[s[1]] == t.vertex_qp[s[2]] == t.vertex_qp[s[3]] == t.vertex_qp[s[4]]
        end
    end
    # duplicate / outside points rejected; boundary points (order-2 prism) work
    dup = QuadratureRule{RefHexahedron}([0.5, 0.5], [Vec(0.1, 0.0, 0.0), Vec(0.1, 0.0, 0.0)])
    @test_throws ErrorException FerriteViz.qp_voronoi_tessellation(RefHexahedron, dup)
    outside = QuadratureRule{RefHexahedron}([0.5, 0.5], [Vec(0.1, 0.0, 0.0), Vec(2.0, 0.0, 0.0)])
    @test_throws ErrorException FerriteViz.qp_voronoi_tessellation(RefHexahedron, outside)
    boundary = QuadratureRule{RefHexahedron}([0.5, 0.5], [Vec(0.1, 0.0, 0.0), Vec(1.0, 0.0, 0.0)])
    @test _tessvol(FerriteViz.qp_voronoi_tessellation(RefHexahedron, boundary)) ≈ 8.0 atol = 1e-9

    # Clip of QP data: exactly region-flat, region volumes and cap areas analytic
    grid = generate_grid(Hexahedron, (1, 1, 1))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefHexahedron,1}()); close!(dh)
    qr = QuadratureRule{RefHexahedron}(2)
    qp = FEData(dh, zeros(ndofs(dh))) |> AddQuadraturePointData(qr, [collect(1.0:8.0)]; output=:iv)
    n = Vec(1.0, 0.0, 0.0)
    clipped = qp |> Clip(ClipPlane(n, 0.2))
    A = vec(FerriteViz.point_data(clipped, :iv)[])
    for tri in clipped.all_triangles   # never blended across a Voronoi wall
        i, j, k = convert(Int, tri[1]), convert(Int, tri[2]), convert(Int, tri[3])
        @test A[i] == A[j] == A[k]
    end
    pc = _vv3.(clipped.coords[])
    vols = Dict{Float64,Float64}()
    for s in clipped.simplices
        @test A[s[1]] == A[s[2]] == A[s[3]] == A[s[4]]
        vols[A[s[1]]] = get(vols, A[s[1]], 0.0) + _tetv(pc[s[1]], pc[s[2]], pc[s[3]], pc[s[4]])
    end
    for (q, ξ) in enumerate(Ferrite.getpoints(qr))
        @test vols[Float64(q)] ≈ (ξ[1] < 0 ? 1.0 : 0.2) atol = 1e-6
    end
    capareas = Dict{Float64,Float64}()
    for tri in clipped.all_triangles
        i, j, k = convert(Int, tri[1]), convert(Int, tri[2]), convert(Int, tri[3])
        all(abs(pc[v][1] - 0.2) <= 1e-6 for v in (i, j, k)) || continue
        capareas[A[i]] = get(capareas, A[i], 0.0) + FerriteViz._tri_area(pc[i], pc[j], pc[k])
    end
    @test length(capareas) == 4 && all(a -> isapprox(a, 1.0; atol=1e-6), values(capareas))
    @test solutionplot(clipped; color=:iv) isa Makie.FigureAxisPlot

    # oblique plane through the Voronoi vertex at the cell center: per-region
    # volumes are analytic — an octant with p positive coordinates of its
    # quadrature point keeps 1, 5/6, 1/6, 0 for p = 0, 1, 2, 3
    ob = qp |> Clip(ClipPlane(Vec(1.0, 1.0, 1.0) / sqrt(3.0), 0.0))
    Ao = vec(FerriteViz.point_data(ob, :iv)[])
    pco = _vv3.(ob.coords[])
    volso = Dict{Float64,Float64}()
    for s in ob.simplices
        @test Ao[s[1]] == Ao[s[2]] == Ao[s[3]] == Ao[s[4]]
        volso[Ao[s[1]]] = get(volso, Ao[s[1]], 0.0) +
                          FerriteViz._signed_tet_volume(pco[s[1]], pco[s[2]], pco[s[3]], pco[s[4]])
    end
    for (q, ξ) in enumerate(Ferrite.getpoints(qr))
        p = count(>(0), ξ)
        expected = (1.0, 5 / 6, 1 / 6, 0.0)[p+1]
        @test get(volso, Float64(q), 0.0) ≈ expected atol = 1e-6
    end
end

@testset "ExtractIsosurfaces" begin
    # 3D: level sets of a linear field are exact planes, through hidden
    # interior cells (no holes), also via a derivation chain
    ds, dh, _ = _hexds(3; f=x -> x[1])
    iso = ds |> ExtractIsosurfaces(0.1)
    @test _surface_area(iso) ≈ 4.0 atol = 1e-6
    @test all(abs(p[1] - 0.1) <= 1e-6 for p in _vv3.(iso.coords[]))
    @test all(isapprox.(vec(FerriteViz.point_data(iso, :u)[]), 0.1; atol=1e-4))
    @test isempty(iso.simplices) && !iso.cells_intact
    multi = ds |> ExtractIsosurfaces([-0.5, 0.1, 0.4])
    @test _surface_area(multi) ≈ 12.0 atol = 1e-5
    @test sort(unique(vec(FerriteViz.point_data(multi, :isovalue)[]))) == [-0.5, 0.1, 0.4]
    @test isempty((ds |> ExtractIsosurfaces(7.0)).all_triangles)            # out of range
    @test _surface_area(ds |> ExtractIsosurfaces(1 / 3)) ≈ 4.0 atol = 1e-5  # nodal level
    gridv = generate_grid(Hexahedron, (3, 3, 3))
    dhv = DofHandler(gridv); add!(dhv, :u, Lagrange{RefHexahedron,1}()^3); close!(dhv)
    uv = zeros(ndofs(dhv)); Ferrite.apply_analytical!(uv, dhv, :u, x -> Vec{3}((x[1], 0.0, 0.0)))
    dchain = FEData(dhv, uv) |> Magnitude(input=:u) |> ExtractIsosurfaces(0.5; input=:magnitude)
    @test _surface_area(dchain) ≈ 8.0 atol = 1e-5   # |x1| = 0.5: two planes

    # composition with Clip in both orders
    ic = ds |> Clip(ClipPlane(Vec(0.0, 0.0, 1.0), 0.25)) |> ExtractIsosurfaces(0.1)
    ci = ds |> ExtractIsosurfaces(0.1) |> Clip(ClipPlane(Vec(0.0, 0.0, 1.0), 0.25))
    @test _surface_area(ic) ≈ 2 * 1.25 atol = 1e-5   # confined to the kept volume
    @test _surface_area(ci) ≈ _surface_area(ic) atol = 1e-6
    @test isempty(ci.simplices)                      # a surface has no volume to cap

    # 2D isolines
    grid2 = generate_grid(Quadrilateral, (4, 4))
    dh2 = DofHandler(grid2); add!(dh2, :u, Lagrange{RefQuadrilateral,1}()); close!(dh2)
    u2 = zeros(ndofs(dh2)); Ferrite.apply_analytical!(u2, dh2, :u, x -> x[1] + x[2])
    ds2 = FEData(dh2, u2)
    iso2 = ds2 |> ExtractIsosurfaces(0.0)
    @test _line_length(iso2) ≈ 2 * sqrt(2.0) atol = 1e-6
    @test isempty(iso2.all_triangles)
    @test all(abs(p[1] + p[2]) <= 1e-6 for p in _vv2.(iso2.coords[]))
    @test solutionplot(iso2; color=:isovalue) isa Makie.FigureAxisPlot
    @test solutionplot(ds2 |> ExtractIsosurfaces(99.0)) isa Makie.FigureAxisPlot  # empty
    @test _line_length(ds2 |> ExtractIsosurfaces([-0.5, 0.5])) ≈ 3 * sqrt(2.0) atol = 1e-5
    # Refine rebuilds whole cells and rejects extracted line geometry
    @test_throws ErrorException iso2 |> Refine(1)

    # input validation
    gridm = generate_grid(Quadrilateral, (2, 2))
    dhm = DofHandler(gridm); add!(dhm, :u, Lagrange{RefQuadrilateral,1}()^2); close!(dhm)
    dsm = FEData(dhm, rand(ndofs(dhm)))
    @test_throws ErrorException dsm |> ExtractIsosurfaces(0.5)             # vector input
    @test_throws ErrorException dsm |> ExtractIsosurfaces(0.5; input=:nope)
    @test_throws ErrorException dsm |> ExtractIsosurfaces(Float64[])
    @test_throws ErrorException dsm |> ExtractIsosurfaces(NaN)
    set_cell_data!(dsm, :cv, rand(getncells(gridm)))
    @test_throws ErrorException dsm |> ExtractIsosurfaces(0.5; input=:cv)  # cell data
    lgrid = generate_grid(Line, (3,))
    ldh = DofHandler(lgrid); add!(ldh, :u, Lagrange{RefLine,1}()); close!(ldh)
    @test_throws ErrorException FEData(ldh, rand(ndofs(ldh))) |> ExtractIsosurfaces(0.5)

    # partial subdomains march only where the field lives (NaN elsewhere)
    dhs = DofHandler(gridm)
    sdh = SubDofHandler(dhs, Set(1:2)); add!(sdh, :p, Lagrange{RefQuadrilateral,1}()); close!(dhs)
    miso = FEData(dhs, collect(1.0:ndofs(dhs))) |> ExtractIsosurfaces(2.5; input=:p)
    @test !any(isnan, reduce(vcat, [collect(p) for p in miso.coords[]]; init=Float64[]))

    # offset/scale robustness of the field-space tolerance
    grids = generate_grid(Quadrilateral, (3, 3))
    dhss = DofHandler(grids); add!(dhss, :u, Lagrange{RefQuadrilateral,1}()); close!(dhss)
    ub = zeros(ndofs(dhss)); Ferrite.apply_analytical!(ub, dhss, :u, x -> 1e6 * (x[1] + x[2]) + 3e8)
    @test _line_length(FEData(dhss, ub) |> ExtractIsosurfaces(3e8)) ≈ 2 * sqrt(2.0) atol = 1e-5
    us = zeros(ndofs(dhss)); Ferrite.apply_analytical!(us, dhss, :u, x -> 1e-7 * (x[1] + x[2]))
    @test _line_length(FEData(dhss, us) |> ExtractIsosurfaces(0.0)) ≈ 2 * sqrt(2.0) atol = 1e-5

    # update! flows into extracted vertices' dof data
    dsr = FEData(dhss, ub)
    riso = dsr |> ExtractIsosurfaces(3e8)
    before = copy(FerriteViz.point_data(riso, :u)[])
    FerriteViz.update!(dsr, 2 .* ub)
    @test FerriteViz.point_data(riso, :u)[] ≈ 2 .* before
end

include("adaptive.jl")

@testset "source hygiene" begin
    src = joinpath(@__DIR__, "..", "src")
    for f in readdir(src; join=true)
        content = read(f, String)
        @test !occursin("@show", content)
        @test !occursin("@info", content)
        # src must not depend on a concrete backend (the backend is the user's
        # choice); detecting the active one by name — as the CairoMakie mesh
        # shim does — is fine, so forbid imports rather than any mention.
        for backend in ("GLMakie", "CairoMakie", "WGLMakie")
            @test !occursin("using $backend", content)
            @test !occursin("import $backend", content)
        end
    end
end
