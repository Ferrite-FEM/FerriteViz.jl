# Tests for view-adaptive tessellation (src/adaptive.jl): the FEData↔isubd
# glue and the solutionplot compute-graph wiring. count_tjunctions comes from
# test/isubd.jl (included earlier by runtests.jl).

# scale only the x/y rows of a projective matrix: a uniform scale would cancel
# in the perspective divide
const ZOOM6 = Makie.Mat4f(6, 0, 0, 0, 0, 6, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)

@testset "adaptive solutionplot: graph wiring" begin
    grid = generate_grid(QuadraticQuadrilateral, (4, 4))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = [0.05 * sin(3i) for i in 1:ndofs(dh)]
    ds = FEData(dh, u)

    fig, ax, sp = solutionplot(ds; adaptive=true, px_target=40.0, max_depth=8)
    nbase = 4 * Ferrite.getncells(grid)  # quad center-fan: 4 base triangles/cell
    keys0 = length(sp.subd_keys[])
    @test keys0 >= nbase
    @test length(sp.subd_positions[]) == 3 * keys0 == length(sp.subd_color[])
    @test length(sp.subd_faces[]) == keys0

    # static camera: pulling again must not recompute (decode identity stable)
    p1 = sp.subd_positions[]
    @test sp.subd_positions[] === p1

    # solution update on an unwarped dataset: colors move, geometry untouched
    c1 = copy(sp.subd_color[])
    FerriteViz.update!(ds, 2 .* u)
    @test sp.subd_color[] != c1
    @test sp.subd_positions[] === p1

    # camera zoom-in refines, zoom-out merges back
    cam = Makie.camera(Makie.parent_scene(sp))
    cam.projectionview[] = ZOOM6 * cam.projectionview[]
    keys_zoomed = length(sp.subd_keys[])
    @test keys_zoomed > keys0
    @test length(sp.subd_positions[]) == 3 * keys_zoomed == length(sp.subd_color[])
    cam.projectionview[] = Makie.Mat4f(inv(Makie.Mat4d(ZOOM6)) * Makie.Mat4d(cam.projectionview[]))
    @test length(sp.subd_keys[]) == keys0

    # the depth cap holds under extreme zoom
    cam.projectionview[] = ZOOM6 * ZOOM6 * ZOOM6 * cam.projectionview[]
    @test all(k -> FerriteViz.key_depth(k) <= 8, sp.subd_keys[])

    # plain color passes through
    figp, axp, spp = solutionplot(ds; adaptive=true, color=:red)
    @test !(spp.plots[1].color[] isa AbstractVector)
    # named point-data arrays cannot be resampled at refined vertices
    FerriteViz.set_point_data!(ds, :pd, rand(FerriteViz.num_vertices(ds)))
    @test_throws ErrorException solutionplot(ds; adaptive=true, color=:pd)
end

@testset "adaptive solutionplot: cracks and values on linear geometry" begin
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    # p(x, y) = x + y is reproduced exactly by the P1 field at any ξ
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> x[1] + x[2])
    ds = FEData(dh, u)

    fig, ax, sp = solutionplot(ds; adaptive=true, px_target=30.0, max_depth=6)
    cam = Makie.camera(Makie.parent_scene(sp))
    cam.projectionview[] = ZOOM6 * cam.projectionview[]
    @test length(sp.subd_keys[]) > 4 * Ferrite.getncells(grid)

    mesh = (positions=sp.subd_positions[], faces=sp.subd_faces[])
    @test count_tjunctions(mesh; digits=5) == 0
    # adaptive colors are exact field evaluations at the refined vertices
    pos, col = sp.subd_positions[], sp.subd_color[]
    @test all(isapprox(col[i], pos[i][1] + pos[i][2]; atol=1e-4) for i in eachindex(col))
end

@testset "adaptive solutionplot: warp drives geometry" begin
    grid = generate_grid(QuadraticQuadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = [0.02 * cos(2i) for i in 1:ndofs(dh)]
    ds = FEData(dh, u)
    wds = ds |> WarpByVector(:u, 1.5)

    # provenance: the warp records itself, filters pass it through
    @test length(wds.deformation) == 1
    @test wds.deformation[1][1][] == :u && wds.deformation[1][2][] == 1.5
    @test length((wds |> Gradient(:u)).deformation) == 1     # _rebind path
    @test isempty(ds.deformation)

    figw, axw, spw = solutionplot(wds; adaptive=true, px_target=40.0, max_depth=8)
    p1 = copy(spw.subd_positions[])
    FerriteViz.update!(ds, 3 .* u)
    p2 = spw.subd_positions[]
    @test p1 != p2
    # displaced by scale * field: check one vertex against a direct evaluation
    ev = FerriteViz.FieldEvaluator(dh, :u)
    ξ1 = spw.subd_ξ[][1]
    cell1 = Ferrite.getcells(grid)[1]
    coords1 = Ferrite.getcoordinates(grid, 1)
    base_x = FerriteViz.geometric_map(Ferrite.geometric_interpolation(typeof(cell1)), coords1, ξ1)
    disp = FerriteViz.evaluate_at(ev, 1, coords1, ξ1, wds.u[])
    expected = base_x + 1.5 * disp
    @test isapprox(collect(p2[1]), collect(expected); atol=1e-4)
end
