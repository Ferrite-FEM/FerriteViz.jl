# Tests for error-adaptive tessellation (src/adaptive.jl): the FEData↔isubd
# glue and the solutionplot compute-graph wiring. count_tjunctions comes from
# test/isubd.jl (included earlier by runtests.jl).

@testset "adaptive solutionplot: error-driven graph wiring" begin
    grid = generate_grid(QuadraticQuadrilateral, (4, 4))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> sin(pi * x[1]) * sin(pi * x[2]))
    ds = FEData(dh, u)

    nbase = 4 * Ferrite.getncells(grid)  # quad center-fan: 4 base triangles/cell
    fig, ax, sp = solutionplot(ds; adaptive=true, solution_tol=5e-3, max_depth=8)
    keys0 = length(sp.subd_keys[])
    @test keys0 > nbase                  # the curved field asks for refinement
    @test length(sp.subd_positions[]) == 3 * keys0 == length(sp.subd_color[])
    @test length(sp.subd_faces[]) == keys0

    # no camera in the graph unless px_target is set
    @test !haskey(sp.attributes.outputs, :subd_projectionview)

    # nothing recomputes without an update (decode identity stable)
    p1 = sp.subd_positions[]
    @test sp.subd_positions[] === p1

    # a looser tolerance coarsens, a tighter one refines
    Makie.update!(sp, solution_tol=0.2)
    keys_loose = length(sp.subd_keys[])
    @test keys_loose < keys0
    Makie.update!(sp, solution_tol=5e-4)
    @test length(sp.subd_keys[]) > keys0
    Makie.update!(sp, solution_tol=5e-3)

    # update! drives colors and (for a changed field shape) the refinement
    c1 = copy(sp.subd_color[])
    FerriteViz.update!(ds, 2 .* u)
    @test sp.subd_color[] != c1

    # the depth cap holds
    Makie.update!(sp, solution_tol=1e-9)
    @test all(k -> FerriteViz.key_depth(k) <= 8, sp.subd_keys[])

    # plain color: no solution estimator, flat straight grid ⇒ no refinement
    figp, axp, spp = solutionplot(ds; adaptive=true, color=:red)
    @test !(spp.plots[1].color[] isa AbstractVector)
    @test length(spp.subd_keys[]) == nbase
    # named point-data arrays cannot be resampled at refined vertices
    FerriteViz.set_point_data!(ds, :pd, rand(FerriteViz.num_vertices(ds)))
    @test_throws ErrorException solutionplot(ds; adaptive=true, color=:pd)
end

@testset "adaptive solutionplot: estimator exactness on linear data" begin
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    # p(x, y) = x + y: linear triangles reproduce both geometry and field
    # exactly, so neither estimator asks for anything
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> x[1] + x[2])
    ds = FEData(dh, u)
    fig, ax, sp = solutionplot(ds; adaptive=true, geometry_tol=1e-9, solution_tol=1e-9)
    @test length(sp.subd_keys[]) == 4 * Ferrite.getncells(grid)
end

@testset "adaptive solutionplot: solution error refines, crack-free, exact colors" begin
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())   # Q2 field on straight cells
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> x[1]^2 + 0.5 * x[2]^2)
    ds = FEData(dh, u)

    fig, ax, sp = solutionplot(ds; adaptive=true, solution_tol=1e-3, max_depth=8)
    nbase = 4 * Ferrite.getncells(grid)
    @test length(sp.subd_keys[]) > nbase
    # geometry alone would not have refined the straight grid
    figg, axg, spg = solutionplot(ds; adaptive=true, solution_tol=1.0, geometry_tol=1e-9)
    @test length(spg.subd_keys[]) == nbase

    mesh = (positions=sp.subd_positions[], faces=sp.subd_faces[])
    @test count_tjunctions(mesh; digits=5) == 0
    # colors are exact field evaluations at the refined vertices
    pos, col = sp.subd_positions[], sp.subd_color[]
    @test all(isapprox(col[i], pos[i][1]^2 + 0.5 * pos[i][2]^2; atol=2e-4) for i in eachindex(col))
end

@testset "adaptive solutionplot: geometry error from warps" begin
    grid = generate_grid(QuadraticQuadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u,
                              x -> Ferrite.Vec(0.1 * sin(pi * x[1]) * sin(pi * x[2]), 0.05 * x[1]^2))
    ds = FEData(dh, u)
    wds = ds |> WarpByVector(:u, 1.5)

    # provenance: the warp records itself, filters pass it through
    @test length(wds.deformation) == 1
    @test wds.deformation[1][1][] == :u && wds.deformation[1][2][] == 1.5
    @test length((wds |> Gradient(:u)).deformation) == 1     # _rebind path
    @test isempty(ds.deformation)

    # plain color: refinement here is purely the geometry estimator seeing the warp
    figw, axw, spw = solutionplot(wds; adaptive=true, color=:red, geometry_tol=2e-3, max_depth=8)
    @test length(spw.subd_keys[]) > 4 * Ferrite.getncells(grid)

    # warped geometry follows the solution
    p1 = copy(spw.subd_positions[])
    FerriteViz.update!(ds, 3 .* u)
    p2 = spw.subd_positions[]
    @test p1 != p2
    # displaced by scale * field: check one vertex against a direct evaluation
    ev = FerriteViz.FieldEvaluator(dh, :u)
    ξ1 = spw.subd_ξ[][1]
    k1 = spw.subd_keys[][1]
    cell1_id = 1 + (FerriteViz.key_base(k1) - 1) ÷ 4
    coords1 = Ferrite.getcoordinates(grid, cell1_id)
    gip = Ferrite.geometric_interpolation(typeof(Ferrite.getcells(grid)[cell1_id]))
    base_x = FerriteViz.geometric_map(gip, coords1, ξ1)
    disp = FerriteViz.evaluate_at(ev, cell1_id, coords1, ξ1, wds.u[])
    expected = base_x + 1.5 * disp
    @test isapprox(collect(p2[1]), collect(expected); atol=1e-4)
end

@testset "adaptive solutionplot: watertight conforming refinement" begin
    grid = generate_grid(Quadrilateral, (6, 6))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())
    close!(dh)
    # a sharp peak at the origin: refinement is strongly non-uniform, so
    # refinement-level boundaries (where T-vertices appear) really do occur
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> exp(-40 * (x[1]^2 + x[2]^2)))
    ds = FEData(dh, u)

    # count T-vertices and edges that are drawn by only one triangle; on this
    # straight-sided grid both are exact
    function survey(sp)
        pos, faces = sp.subd_positions[], sp.subd_faces[]
        rk(p) = (round(Float64(p[1]); digits = 7), round(Float64(p[2]); digits = 7))
        verts = Set(rk(p) for p in pos)
        tj = 0
        counts = Dict{Tuple{Any,Any},Int}()
        for f in faces, (i, j) in ((1, 2), (2, 3), (3, 1))
            a, b = pos[f[i]], pos[f[j]]
            m = rk((a + b) / 2)
            (m in verts && m != rk(a) && m != rk(b)) && (tj += 1)
            e = rk(a) <= rk(b) ? (rk(a), rk(b)) : (rk(b), rk(a))
            counts[e] = get(counts, e, 0) + 1
        end
        border(p) = any(t -> isapprox(p[1], t; atol = 1e-6) || isapprox(p[2], t; atol = 1e-6), (-1.0, 1.0))
        holes = count(((e, c),) -> c == 1 && !(border(e[1]) && border(e[2])), counts)
        over = count(==(3), values(counts)) + count(>(3), values(counts))
        return (; tj, holes, over, n = length(faces))
    end

    kw = (; adaptive = true, solution_tol = 2e-3, max_depth = 9)
    _, _, conf = solutionplot(ds; kw..., conforming = true)
    _, _, free = solutionplot(ds; kw..., conforming = false)

    # the refinement really is non-uniform (otherwise the test proves nothing)
    depths = [FerriteViz.key_depth(k) for k in conf.subd_keys[]]
    @test maximum(depths) - minimum(depths) >= 3

    # conforming: no T-vertices, no holes, no overlaps — watertight
    s = survey(conf)
    @test s.tj == 0 && s.holes == 0 && s.over == 0
    # non-conforming: the level boundaries leave T-vertices behind
    @test survey(free).tj > 0
    # and conformity is cheap: forced splits add a modest number of triangles
    @test s.n < 1.5 * survey(free).n

    # the base adjacency is exact and symmetric
    base, _, _ = FerriteViz._isubd_base(ds)
    @test FerriteViz.is_conformable(base)
    for b in 1:length(base.corners)
        n = FerriteViz.diamond_partner(base, FerriteViz.root_key(b))
        n === nothing && continue
        @test FerriteViz.diamond_partner(base, n) == FerriteViz.root_key(b)
    end
    # every quad contributes a four-triangle fan whose split edges are the
    # element edges: interior element edges pair, boundary ones do not
    nbound = count(b -> base.adjacency[b][FerriteViz.EDGE_S][1] == 0, 1:length(base.corners))
    @test nbound == 4 * 6      # the 24 boundary edges of a 6×6 grid

    # coarsening returns to the base tessellation, still conforming
    Makie.update!(conf, solution_tol = 10.0)
    @test length(conf.subd_keys[]) == length(base.corners)
end

@testset "adaptive solutionplot: optional screen-space criterion" begin
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    ds = FEData(dh, rand(ndofs(dh)))
    fig, ax, sp = solutionplot(ds; adaptive=true, px_target=30.0, max_depth=8)
    @test haskey(sp.attributes.outputs, :subd_projectionview)
    keys0 = length(sp.subd_keys[])
    cam = Makie.camera(Makie.parent_scene(sp))
    pv0 = cam.projectionview[]
    zoom6 = Makie.Mat4f(6, 0, 0, 0, 0, 6, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
    cam.projectionview[] = zoom6 * pv0
    @test length(sp.subd_keys[]) > keys0
    # Exact restore: the state coarsens to a fixed point again, but not
    # necessarily to keys0 — the error estimators are not monotone in depth,
    # so detail discovered while zoomed in (whose deviation genuinely exceeds
    # the tolerance) is kept rather than discarded. See the DeviationLoD
    # docstring; the retained state is the more accurate fixed point.
    cam.projectionview[] = pv0
    krestored = sort(copy(sp.subd_keys[]))
    @test length(krestored) >= keys0
    @test sort(copy(sp.subd_keys[])) == krestored          # stable
    # and the retained keys still tile the domain exactly (no overlap, no holes)
    ξs = sp.subd_ξ[]
    area = sum(abs((ξs[3i - 1][1] - ξs[3i - 2][1]) * (ξs[3i][2] - ξs[3i - 2][2]) -
                   (ξs[3i][1] - ξs[3i - 2][1]) * (ξs[3i - 1][2] - ξs[3i - 2][2])) / 2
               for i in 1:length(sp.subd_keys[]))
    @test area ≈ 4.0 * Ferrite.getncells(grid) rtol = 1e-6  # each ref quad has area 4
end
