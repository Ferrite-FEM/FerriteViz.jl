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
    fine = FerriteViz._cell_tessellation(base, 0, 2)
    @test FerriteViz.ntriangles(fine) == FerriteViz.ntriangles(base)
    @test FerriteViz.nedges(fine) == 4 * FerriteViz.nedges(base)
end

@testset "wireframe edges follow the pipeline" begin
    # curved wireframe: a quadratic displacement on a linear grid bends edges
    grid = generate_grid(Quadrilateral, (1,1))
    dh = DofHandler(grid); add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2); close!(dh)
    g_ana(x) = Vec{2}((0.0, x[1]^2))
    u = zeros(ndofs(dh)); Ferrite.apply_analytical!(u, dh, :u, g_ana)
    ds = FEData(dh, u)
    @test length(ds.all_edges) == 4 * 2^3    # auto edge_resolution = 3 for order 2
    warped = ds |> WarpByVector(:u)
    pts = warped.coords[][FerriteViz._visible_edge_indices(warped)]
    # every warped edge vertex satisfies y = y₀ + x² exactly (up to Float32);
    # a straight-chord wireframe would interpolate linearly between corners
    for (p0, p) in zip(ds.coords[][FerriteViz._visible_edge_indices(ds)], pts)
        @test isapprox(p[2], p0[2] + p0[1]^2; atol=1e-5)
    end
    @test meshplot(warped) isa Makie.FigureAxisPlot

    # adaptive=false opts out of the automatic subdivision entirely
    ds0 = FEData(dh, u; adaptive=false)
    @test length(ds0.all_triangles) == 4 && length(ds0.all_edges) == 4
    # ... and the Refine filter reintroduces it with explicit control
    @test length((ds0 |> Refine(1)).all_triangles) == 16
    dsl = FEData(DofHandler(grid), Float64[]; adaptive=false) |> Refine(1) # linear grid, forced
    @test length(dsl.all_triangles) == 16
    @test length((ds0 |> Refine(surface=0, edges=2)).all_edges) == 4 * 4
    # the automatic filter mode reproduces the constructor default exactly
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
