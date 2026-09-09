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
    ds.adaptivity.solution_tol[] = 5e-3
    ds.adaptivity.max_depth[] = 8
    fig, ax, sp = solutionplot(ds)
    keys0 = length(sp.subd_keys[])
    @test keys0 > nbase                  # the curved field asks for refinement
    @test length(sp.subd_positions[]) == length(sp.subd_color[]) == length(sp.subd_ξ[])
    @test length(sp.subd_positions[]) < 3 * keys0        # vertices are shared within a cell
    @test length(sp.subd_faces[]) == keys0

    # the camera is never an input: refinement is error-driven only
    @test !haskey(sp.attributes.outputs, :subd_projectionview)

    # nothing recomputes without an update (decode identity stable)
    p1 = sp.subd_positions[]
    @test sp.subd_positions[] === p1

    # a looser tolerance coarsens, a tighter one refines
    ds.adaptivity.solution_tol[] = 0.2
    keys_loose = length(sp.subd_keys[])
    @test keys_loose < keys0
    ds.adaptivity.solution_tol[] = 5e-4
    @test length(sp.subd_keys[]) > keys0
    ds.adaptivity.solution_tol[] = 5e-3

    # update! drives colors and (for a changed field shape) the refinement
    c1 = copy(sp.subd_color[])
    FerriteViz.update!(ds, 2 .* u)
    @test sp.subd_color[] != c1

    # the depth cap holds
    ds.adaptivity.solution_tol[] = 1e-9
    @test all(k -> FerriteViz.key_depth(k) <= 8, sp.subd_keys[])

    # plain color: no solution estimator, flat straight grid ⇒ no refinement
    figp, axp, spp = solutionplot(ds; color=:red)
    @test !(spp.plots[1].color[] isa AbstractVector)
    @test length(spp.subd_keys[]) == nbase
    # named point-data arrays cannot be resampled at refined vertices; the
    # plot falls back to the static tessellation instead of erroring
    FerriteViz.set_point_data!(ds, :pd, rand(FerriteViz.num_vertices(ds)))
    _, _, spd = solutionplot(ds; color=:pd)
    @test !haskey(spd.attributes.outputs, :subd_keys)
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
    ds.adaptivity.geometry_tol[] = 1e-9
    ds.adaptivity.solution_tol[] = 1e-9
    fig, ax, sp = solutionplot(ds)
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

    ds.adaptivity.solution_tol[] = 1e-3
    ds.adaptivity.max_depth[] = 8
    fig, ax, sp = solutionplot(ds)
    nbase = 4 * Ferrite.getncells(grid)
    @test length(sp.subd_keys[]) > nbase
    # geometry alone would not have refined the straight grid
    ds.adaptivity.solution_tol[] = 1.0
    ds.adaptivity.geometry_tol[] = 1e-9
    figg, axg, spg = solutionplot(ds)
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

    # provenance: the warp records itself (with the dof handler and solution
    # of the stage it was applied to), filters pass it through
    @test length(wds.deformation) == 1
    @test wds.deformation[1].field[] == :u && wds.deformation[1].scale[] == 1.5
    @test wds.deformation[1].dh === dh && wds.deformation[1].u === ds.u
    @test length((wds |> Gradient(:u)).deformation) == 1     # _rebind path
    @test isempty(ds.deformation)

    # plain color: refinement here is purely the geometry estimator seeing the warp
    wds.adaptivity.geometry_tol[] = 2e-3
    wds.adaptivity.max_depth[] = 8
    figw, axw, spw = solutionplot(wds; color=:red)
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
    disp = FerriteViz.evaluate_at(ev, cell1_id, ξ1, wds.u[])
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
        rk(p) = (round(Float64(p[1]); digits = 7) + 0.0, round(Float64(p[2]); digits = 7) + 0.0)
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

    ds.adaptivity.solution_tol[] = 2e-3
    ds.adaptivity.max_depth[] = 9
    _, _, conf = solutionplot(ds)

    # the refinement really is non-uniform (otherwise the test proves nothing:
    # T-vertices only ever appear at refinement-level boundaries — the
    # unconstrained control lives in test/isubd.jl, at the core level)
    depths = [FerriteViz.key_depth(k) for k in conf.subd_keys[]]
    @test maximum(depths) - minimum(depths) >= 3

    # no T-vertices, no holes, no overlaps — watertight
    s = survey(conf)
    @test s.tj == 0 && s.holes == 0 && s.over == 0

    # the base adjacency is exact and symmetric
    base = FerriteViz._substrate(ds).base
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
    ds.adaptivity.solution_tol[] = 10.0
    @test length(conf.subd_keys[]) == length(base.corners)
end

@testset "adaptive solutionplot: 3D surface is a closed manifold" begin
    grid = generate_grid(Hexahedron, (3, 3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefHexahedron,1}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> exp(-6 * ((x[1] - 1)^2 + x[2]^2 + x[3]^2)))
    ds = FEData(dh, u)
    @test all(ds.solid)                      # an unclipped body is solid everywhere

    # the base is built from the *surface* facets only — six per corner cell of
    # a 3×3×3 block down to one per face-centre cell, fanned into four each
    base = FerriteViz._substrate(ds).base
    nfacets = 6 * 9    # the cube's six sides, nine cells each
    @test length(base.corners) == 4 * nfacets
    @test length(base.corners) < length(ds.all_triangles)   # far fewer than the static path
    # a closed surface: every base triangle has all three neighbours, and every
    # element edge on it carries a diamond
    @test all(t -> all(e -> e[1] != 0, base.adjacency[t]), 1:length(base.corners))
    @test all(t -> FerriteViz.diamond_partner(base, FerriteViz.root_key(t)) !== nothing,
              1:length(base.corners))

    # every drawn edge is shared by exactly two triangles: watertight and closed
    function survey(sp)
        pos, faces = sp.subd_positions[], sp.subd_faces[]
        rk(p) = ntuple(i -> round(Float64(p[i]); digits = 6) + 0.0, 3)
        verts = Set(rk(p) for p in pos)
        tj = 0
        counts = Dict{Any,Int}()
        for f in faces, (i, j) in ((1, 2), (2, 3), (3, 1))
            a, b = pos[f[i]], pos[f[j]]
            m = rk((a + b) / 2)
            (m in verts && m != rk(a) && m != rk(b)) && (tj += 1)
            e = rk(a) <= rk(b) ? (rk(a), rk(b)) : (rk(b), rk(a))
            counts[e] = get(counts, e, 0) + 1
        end
        return (; tj, open = count(==(1), values(counts)), over = count(>(2), values(counts)),
                n = length(faces))
    end

    ds.adaptivity.solution_tol[] = 5e-3
    ds.adaptivity.max_depth[] = 5
    _, _, conf = solutionplot(ds)
    s = survey(conf)
    @test s.tj == 0 && s.open == 0 && s.over == 0
    @test s.n < length(ds.all_triangles)
    # the refinement is non-uniform, so watertightness above is not vacuous
    depths = [FerriteViz.key_depth(k) for k in conf.subd_keys[]]
    @test maximum(depths) > minimum(depths)
end

@testset "adaptive solutionplot: CrinkleClip keeps the surface closed" begin
    grid = generate_grid(Hexahedron, (3, 3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefHexahedron,1}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> x[1] + x[2] + x[3])
    ds = FEData(dh, u)
    clipped = ds |> CrinkleClip(FerriteViz.ClipPlane(Ferrite.Vec(1.0, 1.0, 1.0), 0.3))

    # the clip records the body it kept; filters pass the record on
    @test count(clipped.solid) < Ferrite.getncells(grid)
    @test all(clipped.visible .<= clipped.solid)    # a visible cell is part of the body
    @test isempty(ds.deformation) && all(ds.solid)  # the input is untouched
    @test (clipped |> Gradient(:p)).solid == clipped.solid

    # the cut exposes new surface, and it is closed just like the outer one
    base = FerriteViz._substrate(clipped).base
    @test all(t -> all(e -> e[1] != 0, base.adjacency[t]), 1:length(base.corners))
    clipped.adaptivity.solution_tol[] = 5e-3
    clipped.adaptivity.max_depth[] = 4
    _, _, sp = solutionplot(clipped)
    pos, faces = sp.subd_positions[], sp.subd_faces[]
    rk(p) = ntuple(i -> round(Float64(p[i]); digits = 6) + 0.0, 3)
    counts = Dict{Any,Int}()
    for f in faces, (i, j) in ((1, 2), (2, 3), (3, 1))
        a, b = rk(pos[f[i]]), rk(pos[f[j]])
        e = a <= b ? (a, b) : (b, a)
        counts[e] = get(counts, e, 0) + 1
    end
    @test count(==(1), values(counts)) == 0 && count(>(2), values(counts)) == 0
end

# A custom 2D cell whose reference tessellation lists no element edges — the
# extension contract allows that (meshplot then draws no wireframe, cf. the
# cohesive-cell docs example). No conforming base exists without the edges,
# so such datasets must keep the static path rather than error.
struct RefNoEdgeQuad <: Ferrite.AbstractRefShape{2} end
struct NoEdgeQuadrilateral <: Ferrite.AbstractCell{RefNoEdgeQuad}
    nodes::NTuple{4,Int}
end
Ferrite.geometric_interpolation(::Type{NoEdgeQuadrilateral}) = Lagrange{RefQuadrilateral,1}()
FerriteViz.reference_tessellation(::Type{RefNoEdgeQuad}) =
    FerriteViz.ReferenceTessellation(Ferrite.reference_coordinates(Lagrange{RefQuadrilateral,1}()),
                                     [(1, 2, 3), (1, 3, 4)])

@testset "custom cells without tessellation edges keep the static path" begin
    nodes = [Node((0.0, 0.0)), Node((1.0, 0.0)), Node((1.0, 1.0)), Node((0.0, 1.0)),
             Node((2.0, 0.0)), Node((2.0, 1.0))]
    cells = Ferrite.AbstractCell[Quadrilateral((1, 2, 3, 4)),
                                 NoEdgeQuadrilateral((2, 5, 6, 3))]
    grid = Grid(cells, nodes)
    ds = FEData(DofHandler(grid), Float64[])
    @test ds.adaptivity !== nothing                  # the default config is there
    @test !FerriteViz._adaptive_capable(ds)          # ... but this grid opts out
    _, _, mp = meshplot(ds)
    @test !haskey(mp.attributes.outputs, :subd_keys)
end

# The cohesive-cell docs example: the cell's node numbering (facets (1,2) and
# (3,4)) would make the standard bilinear map fold into a bow-tie, so it
# carries a geometric interpolation that places node i at the right reference
# corner. The reference geometry is then just a quadrilateral — tessellation,
# element edges and the conforming adaptive base included.
struct RefTestCohesive <: Ferrite.AbstractRefShape{2} end
struct TestCohesiveQuad <: Ferrite.AbstractCell{RefTestCohesive}
    nodes::NTuple{4,Int}
end
struct TestCohesiveLagrange <: Ferrite.ScalarInterpolation{RefTestCohesive,1} end
const TCOH_PERM = (1, 2, 4, 3)      # node i sits at reference corner TCOH_PERM[i]
Ferrite.getnbasefunctions(::TestCohesiveLagrange) = 4
Ferrite.reference_shape_value(::TestCohesiveLagrange, ξ, i::Int) =
    Ferrite.reference_shape_value(Lagrange{RefQuadrilateral,1}(), ξ, TCOH_PERM[i])
Ferrite.reference_coordinates(::TestCohesiveLagrange) =
    Ferrite.reference_coordinates(Lagrange{RefQuadrilateral,1}())[collect(TCOH_PERM)]
Ferrite.geometric_interpolation(::Type{TestCohesiveQuad}) = TestCohesiveLagrange()
FerriteViz.reference_tessellation(::Type{RefTestCohesive}) =
    FerriteViz.reference_tessellation(Ferrite.RefQuadrilateral)

@testset "custom cells with a clean parametrization refine adaptively" begin
    # the opened cohesive demo of the docs: two blocks bridged by an interface
    Δ = 0.5
    nodes = [Node((0.0, 0.0)), Node((1.0, 0.0)), Node((1.0, 1.0)), Node((0.0, 1.0)),
             Node((1.0 + Δ, 0.0)), Node((1.0 + Δ, 1.0)), Node((2.0 + Δ, 1.0)), Node((2.0 + Δ, 0.0))]
    cells = Ferrite.AbstractCell[Quadrilateral((1, 2, 3, 4)), Quadrilateral((5, 6, 7, 8)),
                                 TestCohesiveQuad((2, 3, 5, 6))]
    grid = Grid(cells, nodes)
    ds = FEData(DofHandler(grid), Float64[])
    @test FerriteViz._adaptive_capable(ds)
    _, _, mp = meshplot(ds)
    @test haskey(mp.attributes.outputs, :subd_keys)
    # every map is affine here, so nothing refines beyond the base fans
    @test length(mp.subd_keys[]) == 4 * 3

    # the interface facets pair with the neighbouring blocks by node ids, so
    # the base is conforming across the interface: of the 12 fan triangles,
    # two diamonds (4 triangles) sit on the interior facets
    base = FerriteViz._substrate(ds).base
    @test FerriteViz.is_conformable(base)
    nbound = count(b -> base.adjacency[b][FerriteViz.EDGE_S][1] == 0, 1:length(base.corners))
    @test nbound == 8

    # the permuted parametrization maps the reference quad onto the physical
    # interface without folding: reference corner (1,1) carries node 4 of the
    # cell (global node 6), not node 3 — a plain bilinear on the cohesive
    # node order (or coefficients fit against the wrong nodal layout) puts
    # the opposite node there
    b = findfirst(==(3), FerriteViz._substrate(ds).cellmap)
    @test base.mapping(b, Ferrite.Vec(1.0f0, 1.0f0)) ≈ [1.5, 1.0] atol = 1e-5
    @test base.mapping(b, Ferrite.Vec(-1.0f0, -1.0f0)) ≈ [1.0, 0.0] atol = 1e-5
    @test base.mapping(b, Ferrite.Vec(0.0f0, 0.0f0)) ≈ [1.25, 0.5] atol = 1e-5
end

@testset "Adaptivity is a filter" begin
    grid = generate_grid(QuadraticQuadrilateral, (2, 2))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> sin(pi * x[1]) * x[2])

    # keyword sugar: true = defaults, an instance = automatic application
    @test FEData(dh, u).adaptivity isa Adaptivity
    @test FEData(dh, u; adaptivity=false).adaptivity === nothing
    a = Adaptivity(solution_tol=1e-2, max_depth=6)
    @test FEData(dh, u; adaptivity=a).adaptivity === a

    # the filter configures a static dataset without touching the source
    ds = FEData(dh, u; adaptivity=false)
    ads = ds |> Adaptivity(solution_tol=1e-2)
    @test ds.adaptivity === nothing
    @test ads.adaptivity.solution_tol[] == 1e-2
    @test ads.coords === ds.coords               # geometry is shared
    @test ads.subd_cache[] === nothing           # fresh substrate cache
    _, _, sp = solutionplot(ads)
    @test haskey(sp.attributes.outputs, :subd_keys)
    _, _, ss = solutionplot(ds)
    @test !haskey(ss.attributes.outputs, :subd_keys)

    # ... and replaces the settings of an adaptive one
    bds = ads |> Adaptivity(max_depth=4)
    @test bds.adaptivity.max_depth[] == 4
    @test bds.adaptivity !== ads.adaptivity
end

@testset "adaptive substrate: shared per dataset, settings shared too" begin
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> x[1]^2 + 0.5 * x[2]^2)
    ds = FEData(dh, u)

    @test ds.subd_cache[] === nothing            # built lazily
    _, _, sp1 = solutionplot(ds)
    _, _, sp2 = solutionplot(ds)
    sub = FerriteViz._substrate(ds)
    @test ds.subd_cache[] === sub                # both plots hit one cache
    @test length(sub.evaluators) == 1            # ... and share the :p evaluator
    # adaptivity is a dataset property: one setting steers every plot
    n0 = length(sp1.subd_keys[])
    @test length(sp2.subd_keys[]) == n0
    ds.adaptivity.solution_tol[] = 1e-4
    @test length(sp1.subd_keys[]) > n0
    @test length(sp1.subd_keys[]) == length(sp2.subd_keys[])
    # the retained keys of each plot tile the domain exactly (no overlap/holes)
    for sp in (sp1, sp2)
        ξs, fs = sp.subd_ξ[], sp.subd_faces[]
        area = sum(abs((ξs[f[2]][1] - ξs[f[1]][1]) * (ξs[f[3]][2] - ξs[f[1]][2]) -
                       (ξs[f[3]][1] - ξs[f[1]][1]) * (ξs[f[2]][2] - ξs[f[1]][2])) / 2 for f in fs)
        @test area ≈ 4.0 * Ferrite.getncells(grid) rtol = 1e-6  # each ref quad has area 4
    end
    # a filter stage is a new dataset with its own (empty) cache, but the
    # same shared settings object
    gds = ds |> Gradient(:p)
    @test gds.subd_cache[] === nothing
    @test gds.adaptivity === ds.adaptivity
end

@testset "adaptive: warp observables drive the plot" begin
    grid = generate_grid(Quadrilateral, (4, 4))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u, x -> Ferrite.Vec(0.0, 0.2 * sin(pi * (x[1] + 1) / 2)))
    ds = FEData(dh, u)
    scale = Makie.Observable(1.0)
    wds = ds |> WarpByVector(:u, scale)

    wds.adaptivity.geometry_tol[] = 1e-3
    wds.adaptivity.max_depth[] = 8
    _, _, sp = solutionplot(wds; color=:red)
    maxy(ps) = maximum(p -> p[2], ps)
    @test maxy(sp.subd_positions[]) > 1.1        # warped up
    n1 = length(sp.subd_keys[])
    @test n1 > 4 * Ferrite.getncells(grid)       # the warp curves the geometry

    # a slider-driven scale must move the plot, not wait for the next update!
    scale[] = 0.0
    @test maxy(sp.subd_positions[]) ≈ 1.0 atol = 1e-6
    @test length(sp.subd_keys[]) < n1            # flat geometry coarsens
    scale[] = 1.0
    @test maxy(sp.subd_positions[]) > 1.1
    @test length(sp.subd_keys[]) == n1

    # the adaptive wireframe reacts just the same
    wds.adaptivity.geometry_tol[] = 1e-3
    wds.adaptivity.max_depth[] = 8
    _, _, wp = meshplot(wds)
    @test maxy(wp.edge_lines[]) > 1.1
    scale[] = 0.0
    @test maxy(wp.edge_lines[]) ≈ 1.0 atol = 1e-6
    scale[] = 1.0
end

@testset "adaptive: a warp applied before Gradient keeps its dof handler" begin
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u, x -> Ferrite.Vec(0.0, 0.2 * sin(pi * (x[1] + 1) / 2)))
    ds = FEData(dh, u)
    # the gradient dataset's handler has no :u field — the warp must evaluate
    # against the handler and solution captured when it was applied
    gds = ds |> WarpByVector(:u, 1.0) |> Gradient(:u)
    @test :u ∉ Ferrite.getfieldnames(gds.dh)
    gds.adaptivity.geometry_tol[] = 1e-3
    gds.adaptivity.max_depth[] = 8
    _, _, sp = solutionplot(gds; color=:default)
    @test maximum(p -> p[2], sp.subd_positions[]) > 1.1   # the warp is applied
    @test length(sp.subd_keys[]) > 4 * Ferrite.getncells(grid)
end

@testset "adaptive: deviations are memoized on the substrate" begin
    grid = generate_grid(QuadraticQuadrilateral, (4, 4))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> sin(pi * x[1]) * sin(pi * x[2]))
    ds = FEData(dh, u)

    ds.adaptivity.solution_tol[] = 5e-3
    ds.adaptivity.max_depth[] = 8
    _, _, sp1 = solutionplot(ds)
    keys1 = sort(sp1.subd_keys[])
    sub = FerriteViz._substrate(ds)
    # a memo entry is written per evaluated key, so the sizes count evaluations
    ngeo, nsol = length(sub.dev_caches[:geometry]), length(sub.dev_caches[:p])
    @test ngeo > 0 && nsol > 0

    # a second plot decides the same mesh from lookups alone
    _, _, sp2 = solutionplot(ds)
    @test sort(sp2.subd_keys[]) == keys1
    @test length(sub.dev_caches[:geometry]) == ngeo
    @test length(sub.dev_caches[:p]) == nsol

    # loosening the tolerance re-decides without a single new evaluation
    # (every key the merge pass asks about was once a leaf, hence memoized)
    ds.adaptivity.solution_tol[] = 5e-2
    @test length(sp2.subd_keys[]) < length(keys1)
    @test length(sub.dev_caches[:p]) == nsol
    # tightening evaluates only the genuinely new, deeper keys
    ds.adaptivity.solution_tol[] = 5e-4
    @test length(sp2.subd_keys[]) > length(keys1)
    @test length(sub.dev_caches[:p]) > nsol
    # the setting is a dataset property: every plot follows it
    @test sort(sp1.subd_keys[]) == sort(sp2.subd_keys[])

    # a solution update invalidates the memos: cleared, then refilled for the
    # new state — and the plots re-refine against fresh deviations
    n_tight = length(sub.dev_caches[:p])
    FerriteViz.update!(ds, 2 .* u)
    sp1.subd_keys[]; sp2.subd_keys[]
    @test sub.dev_epoch[] == sub.epoch[]
    @test 0 < length(sub.dev_caches[:p]) <= n_tight
end

@testset "adaptive: wireframe and surface share the geometry memo" begin
    grid = generate_grid(Quadrilateral, (4, 4))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u, x -> Ferrite.Vec(0.0, 0.15 * sin(pi * x[1])))
    scale = Makie.Observable(1.0)
    wds = FEData(dh, u) |> WarpByVector(:u, scale)

    wds.adaptivity.geometry_tol[] = 1e-3
    wds.adaptivity.max_depth[] = 8
    _, _, sp = solutionplot(wds; color=:red)
    sp.subd_keys[]
    sub = FerriteViz._substrate(wds)
    ngeo = length(sub.dev_caches[:geometry])
    @test ngeo > 0
    wds.adaptivity.geometry_tol[] = 1e-3
    wds.adaptivity.max_depth[] = 8
    _, _, wp = meshplot(wds)
    @test length(wp.edge_lines[]) > 0
    @test length(sub.dev_caches[:geometry]) == ngeo   # decided from the memo

    # a warp-scale change is a new epoch: the stale geometry deviations are
    # dropped, not reused — the flattened plot re-evaluates and coarsens
    nwarped = length(sp.subd_keys[])
    e0 = sub.epoch[]
    scale[] = 0.0
    @test sub.epoch[] > e0
    @test length(sp.subd_keys[]) < nwarped   # flat geometry coarsens
    @test sub.dev_epoch[] == sub.epoch[]     # the memo was rebuilt for the new state
    scale[] = 1.0
end

@testset "adaptive: derived quantities follow their source fields" begin
    # the elastic chain: Q2 displacement -> Gradient -> Derive(vonmises ∘ σ).
    # The criterion samples the *source* (:gradient); the colors re-apply the
    # recorded closures at the refined vertices.
    grid = generate_grid(Quadrilateral, (4, 4))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    disp(x) = Ferrite.Vec(0.1 * x[1]^2 * x[2], 0.2 * x[1] * x[2]^2)
    σ(∇u) = 2.0 * dev(symmetric(∇u)) + tr(∇u) * one(symmetric(∇u))
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u, disp)
    ds = FEData(dh, u)
    vds = ds |> Gradient(:u) |> Derive(∇u -> FerriteViz.vonmises(σ(∇u)); input=:gradient, output=:σvM)
    @test haskey(vds.point_derivations, :σvM)

    nbase = 4 * Ferrite.getncells(grid)
    vds.adaptivity.solution_tol[] = 1e-3
    vds.adaptivity.max_depth[] = 8
    fig, ax, sp = solutionplot(vds; color=:σvM)
    @test length(sp.subd_keys[]) > nbase          # the source field drove refinement

    # colors are the recorded closures re-applied to the same :gradient field
    # the static path samples: where a refined vertex coincides with a static
    # tessellation vertex, the values must agree (the reference here is the
    # static array, not the analytic gradient — Gradient's projection into its
    # DG space is its own, separately tested concern)
    pos, col = sp.subd_positions[], copy(sp.subd_color[])
    rk(p) = (round(Float64(p[1]); digits=6) + 0.0, round(Float64(p[2]); digits=6) + 0.0)
    static_vals = Dict{Tuple{Float64,Float64},Float64}()
    A = FerriteViz.point_data(vds, :σvM)[]
    for (i, p) in enumerate(vds.coords[])
        static_vals[rk(p)] = A[i, 1]
    end
    matched = 0
    for i in eachindex(col)
        v = get(static_vals, rk(pos[i]), nothing)
        v === nothing && continue
        matched += 1
        @test col[i] ≈ v atol = 1e-5
    end
    @test matched > 50   # plenty of coincident vertices to make that meaningful

    # u -> colors AND tessellation: a pure scaling doubles the (1-homogeneous)
    # colors on the same key set (the solution tolerance is span-relative)...
    k0 = sort(copy(sp.subd_keys[]))
    FerriteViz.update!(ds, 2 .* u)
    @test sort(copy(sp.subd_keys[])) == k0
    @test sp.subd_color[] ≈ 2 .* col rtol = 1e-4
    # ...and a shape change moves the tessellation
    u2 = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u2, dh, :u, x -> Ferrite.Vec(0.0, 0.2 * exp(-30 * sum(abs2, x))))
    FerriteViz.update!(ds, u2)
    @test sort(copy(sp.subd_keys[])) != k0
    FerriteViz.update!(ds, u)

    # Threshold rides the chain: clipped range renders as NaN
    tds = vds |> Threshold(input=:σvM, min=0.5 * maximum(col))
    tds.adaptivity.solution_tol[] = 1e-3
    tds.adaptivity.max_depth[] = 8
    _, _, spt = solutionplot(tds; color=:threshold)
    ct = spt.subd_color[]
    @test any(isnan, ct) && any(!isnan, ct)
    @test all(isnan(ct[i]) || ct[i] >= 0.5 * maximum(col) - 1e-6 for i in eachindex(ct))

    # a tensor-valued derivation cannot color, and says so at plot creation
    dds = ds |> Gradient(:u) |> Deviator(input=:gradient)
    @test_throws ErrorException solutionplot(dds; color=:deviator)

    # a derivation of a raw registered array has no pointwise meaning: no
    # record, so the adaptive path cannot resample it — the plot falls back
    # to the static tessellation
    FerriteViz.set_point_data!(ds, :raw, rand(FerriteViz.num_vertices(ds)))
    rds = ds |> Magnitude(input=:raw)
    @test !haskey(rds.point_derivations, :magnitude)
    _, _, spr = solutionplot(rds; color=:magnitude)
    @test !haskey(spr.attributes.outputs, :subd_keys)
end

@testset "adaptive: solution span comes from the dof values" begin
    # a peak at an edge-midpoint node is invisible to every base-triangle
    # corner (element vertices and fan centres): sampling the span there once
    # collapsed the tolerance reference to ~5e-4, which refined this visually
    # flat field toward the depth cap
    grid = generate_grid(Quadrilateral, (2, 2))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> exp(-30 * (x[1]^2 + (x[2] - 0.5)^2)))
    ds = FEData(dh, u)

    ev = FerriteViz.FieldEvaluator(dh, :p)
    span, mag = FerriteViz._field_span(ev, 1:Ferrite.getncells(grid), u; reduce=false)
    @test span > 0.9
    @test mag >= span
    # relative to the true span, a 50% tolerance asks for almost nothing
    ds.adaptivity.solution_tol[] = 0.5
    ds.adaptivity.max_depth[] = 10
    _, _, sp = solutionplot(ds)
    @test length(sp.subd_keys[]) < 3 * 4 * Ferrite.getncells(grid)
end

@testset "adaptive: cells outside the field's subdomain render as holes" begin
    # every visible cell is tessellated, including cells of a subdomain the
    # color field is not defined on: evaluation there yields `nothing`, which
    # the color transfer maps to NaN (a hole, as in the static path) and the
    # deviation probe to zero (nothing to refine for) — not an error
    grid = generate_grid(Quadrilateral, (4, 2))
    dh = DofHandler(grid)
    sdh1 = SubDofHandler(dh, Set(1:4))
    add!(sdh1, :p, Lagrange{RefQuadrilateral,2}())
    sdh2 = SubDofHandler(dh, Set(5:8))
    add!(sdh2, :q, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    u = collect(range(0.0, 1.0, ndofs(dh)))
    ds = FEData(dh, u)

    _, _, sp = solutionplot(ds; color=:p)
    col = sp.subd_color[]
    @test any(isnan, col)      # the :q half renders as holes
    @test any(!isnan, col)     # the :p half renders values
    ev = FerriteViz._field_evaluator(ds.subd_cache[], dh, :p)
    ξ = Ferrite.Vec(0.25, -0.5)
    @test FerriteViz.evaluate_at(ev, 5, ξ, u) === nothing
    @test FerriteViz.evaluate_at(ev, 1, ξ, u) isa Real
end

@testset "adaptive: the pipeline samples in the render number type" begin
    grid = generate_grid(QuadraticQuadrilateral, (2, 2))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,2}())
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :p, x -> sin(pi * x[1]) * x[2])
    ds = FEData(dh, u)

    _, _, sp = solutionplot(ds)
    sub = ds.subd_cache[]
    @test ds.adaptivity.sample_type === Float32
    @test FerriteViz._sample_type(sub) === Float32
    @test eltype(eltype(eltype(sub.base.corners))) === Float32
    @test sub.base.mapping(1, sub.base.corners[1][1]) isa Tensors.Vec{2,Float32}

    # Float64 stays available as an opt-in — a dataset property, carried
    # through filters
    ds64 = FEData(dh, u; adaptivity=Adaptivity(sample_type=Float64))
    @test (ds64 |> FerriteViz.Refine(1)).adaptivity.sample_type === Float64
    solutionplot(ds64)
    @test FerriteViz._sample_type(ds64.subd_cache[]) === Float64

    # the eps-scaled tolerance floor: a large constant field samples with
    # Float32 noise proportional to its value, which must never read as
    # deviation — however tight the requested tolerance
    uc = fill(5.0, ndofs(dh))
    dsc = FEData(dh, uc)
    nbase = 4 * Ferrite.getncells(grid)
    dsc.adaptivity.solution_tol[] = 1e-9
    dsc.adaptivity.max_depth[] = 8
    _, _, spc = solutionplot(dsc)
    @test length(spc.subd_keys[]) == nbase
end

@testset "per-cell coefficients agree with the shape-function sum" begin
    # every standard Lagrange space is representable in a monomial basis;
    # anything else must fall back rather than fit something wrong
    for (ip, supported) in ((Lagrange{RefTriangle,1}(), true),
                            (Lagrange{RefTriangle,2}(), true),
                            (Lagrange{RefQuadrilateral,1}(), true),
                            (Lagrange{RefQuadrilateral,2}(), true),
                            (Lagrange{RefHexahedron,1}(), true),
                            (Lagrange{RefTetrahedron,2}(), true),
                            (Lagrange{RefPrism,1}(), false))
        @test (FerriteViz.PolyBasis(ip) !== nothing) == supported
    end
    # a vectorized interpolation reuses its scalar basis
    @test FerriteViz.PolyBasis(Lagrange{RefQuadrilateral,2}()^2) !== nothing

    # and the coefficients reproduce the shape-function sum — exactly in
    # Float64, to the sample type's resolution in the Float32 default
    for (celltype, ip) in ((Quadrilateral, Lagrange{RefQuadrilateral,1}()),
                           (QuadraticQuadrilateral, Lagrange{RefQuadrilateral,2}()),
                           (Triangle, Lagrange{RefTriangle,1}()))
        grid = generate_grid(celltype, (3, 3))
        dh = DofHandler(grid); add!(dh, :f, ip); close!(dh)
        u = zeros(ndofs(dh))
        Ferrite.apply_analytical!(u, dh, :f, x -> sin(pi * x[1]) + x[2]^2)
        ev = FerriteViz.FieldEvaluator(dh, :f, Float64)
        ev32 = FerriteViz.FieldEvaluator(dh, :f)
        @test ev.poly !== nothing
        FerriteViz.prepare!(ev, 1:Ferrite.getncells(grid), u)
        FerriteViz.prepare!(ev32, 1:Ferrite.getncells(grid), u)
        inside = Ferrite.getrefshape(ip) === RefTriangle ?
                 (Ferrite.Vec(0.2, 0.3), Ferrite.Vec(0.1, 0.6), Ferrite.Vec(0.0, 0.0)) :
                 (Ferrite.Vec(0.13, -0.21), Ferrite.Vec(-0.7, 0.4), Ferrite.Vec(0.9, 0.9))
        for cell in (1, 5, 9), ξ in inside
            direct = FerriteViz._sum_shape_values(ev.ips[1], ev.celldofs_field[cell], ξ, u)
            @test FerriteViz.evaluate(ev.poly, cell, ξ) ≈ direct atol = 1e-12
            @test FerriteViz.evaluate_at(ev, cell, ξ, u) ≈ direct atol = 1e-12
            @test FerriteViz.evaluate_at(ev32, cell, ξ, u) ≈ direct atol = 1e-5
        end
    end
end

@testset "adaptive QP data: higher-order geometry, crisp regions" begin
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u, x -> Ferrite.Vec(0.0, 0.2 * sin(pi * (x[1] + 1) / 2)))
    qr = Ferrite.QuadratureRule{RefQuadrilateral}(2)         # 2x2 Gauss: 4 regions/cell
    ncells = Ferrite.getncells(grid)
    states = Makie.Observable([Float64[10i + j for j in 1:4] for i in 1:ncells])
    qds = FEData(dh, u) |> AddQuadraturePointData(qr, states) |> WarpByVector(:u, 1.0)
    @test qds.qp_partition !== nothing
    @test length(qds.deformation) == 1

    qds.adaptivity.geometry_tol[] = 5e-4
    qds.adaptivity.max_depth[] = 8
    _, _, sp = solutionplot(qds; color=:qpdata)
    sub = FerriteViz._substrate(qds)
    @test length(unique(sub.groupmap)) == 4 * ncells         # one sharing group per region
    @test length(sp.subd_keys[]) > length(sub.base.corners)  # the warp curves the regions

    # every vertex carries exactly its region's value — nothing is blended
    pos, col = sp.subd_positions[], copy(sp.subd_color[])
    allvals = Float32.(reduce(vcat, states[]))
    @test all(c -> c in allvals, col)
    # crisp jumps: region-boundary positions appear duplicated with different values
    rk(p) = (round(Float64(p[1]); digits=5) + 0.0, round(Float64(p[2]); digits=5) + 0.0)
    bycoord = Dict{Any,Set{Float32}}()
    for i in eachindex(pos)
        push!(get!(Set{Float32}, bycoord, rk(pos[i])), col[i])
    end
    @test count(s -> length(s) >= 2, values(bycoord)) > 10

    # watertight and exactly tiling: no T-vertices, and the reference areas of
    # the drawn triangles sum to the full domain (no holes, no overlaps)
    faces = sp.subd_faces[]
    @test count_tjunctions((positions=pos, faces=faces); digits=5) == 0
    ξs = sp.subd_ξ[]
    area = sum(abs((ξs[f[2]][1] - ξs[f[1]][1]) * (ξs[f[3]][2] - ξs[f[1]][2]) -
                   (ξs[f[3]][1] - ξs[f[1]][1]) * (ξs[f[2]][2] - ξs[f[1]][2])) / 2 for f in faces)
    @test area ≈ 4.0 * ncells rtol = 1e-6

    # updating the internal variables recolors; the mesh has no reason to move
    k0 = length(sp.subd_keys[])
    states[] = [2 .* s for s in states[]]
    @test sp.subd_color[] ≈ 2 .* col rtol = 1e-6
    @test length(sp.subd_keys[]) == k0

    # the geometry answers to the tolerance (a dataset property: every plot
    # of the dataset re-refines), and the surface follows the warp
    qds.adaptivity.geometry_tol[] = 5e-5
    qds.adaptivity.max_depth[] = 10
    @test length(sp.subd_keys[]) > k0
    @test maximum(p -> p[2], pos) > 1.1

    # non-scalar entries need a reducing extract, said at plot creation
    tstates = [[Ferrite.Vec(1.0, 2.0) for _ in 1:4] for _ in 1:ncells]
    tqds = FEData(dh, u) |> AddQuadraturePointData(qr, tstates)
    @test_throws ErrorException solutionplot(tqds; color=:qpdata)
end

@testset "adaptive QP data: wireframe stays on element edges" begin
    # unwarped linear grid: element edges are the grid lines, Voronoi rims of
    # a 2x2 Gauss rule are the cell midlines — the wireframe must draw the
    # former and mask the latter
    grid = generate_grid(Quadrilateral, (3, 3))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefQuadrilateral,1}())
    close!(dh)
    qr = Ferrite.QuadratureRule{RefQuadrilateral}(2)
    states = [rand(4) for _ in 1:Ferrite.getncells(grid)]
    qds = FEData(dh, zeros(ndofs(dh))) |> AddQuadraturePointData(qr, states)
    qds.adaptivity.geometry_tol[] = 1e-3
    qds.adaptivity.max_depth[] = 6
    _, _, wp = meshplot(qds)
    segs = wp.edge_lines[]
    nseg = length(segs) ÷ 2
    @test nseg > 0
    lines = (-1.0, -1 / 3, 1 / 3, 1.0)
    @test all(1:nseg) do i
        a, b = segs[2i - 1], segs[2i]
        (isapprox(a[1], b[1]; atol=1e-6) && any(l -> isapprox(a[1], l; atol=1e-5), lines)) ||
            (isapprox(a[2], b[2]; atol=1e-6) && any(l -> isapprox(a[2], l; atol=1e-5), lines))
    end
end

@testset "adaptive QP data: clipped 3D body renders internal variables" begin
    grid = generate_grid(Hexahedron, (2, 2, 2))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefHexahedron,1}())
    close!(dh)
    qr = Ferrite.QuadratureRule{RefHexahedron}(2)            # 8 QPs
    states = [Float64[10i + j for j in 1:8] for i in 1:Ferrite.getncells(grid)]
    qds = FEData(dh, zeros(ndofs(dh))) |> AddQuadraturePointData(qr, states) |>
          CrinkleClip(FerriteViz.ClipPlane(Ferrite.Vec(1.0, 1.0, 1.0), 0.0))
    # the registered static array survives the clip (the vertex layout is
    # unchanged), and so does the partition record the adaptive path draws from
    @test haskey(qds.point_data, :qpdata)
    @test qds.qp_partition !== nothing
    qds.adaptivity.geometry_tol[] = 1e-3
    qds.adaptivity.max_depth[] = 4
    _, _, sp = solutionplot(qds; color=:qpdata)
    pos, col, faces = sp.subd_positions[], sp.subd_color[], sp.subd_faces[]
    @test length(faces) > 0
    @test all(c -> c in Float32.(reduce(vcat, states)), col)
    # the clipped surface is a closed manifold: every drawn edge shared twice
    rk3(p) = ntuple(i -> round(Float64(p[i]); digits=5) + 0.0, 3)
    counts = Dict{Any,Int}()
    for f in faces, (i, j) in ((1, 2), (2, 3), (3, 1))
        a, b = rk3(pos[f[i]]), rk3(pos[f[j]])
        e = a <= b ? (a, b) : (b, a)
        counts[e] = get(counts, e, 0) + 1
    end
    @test count(==(1), values(counts)) == 0 && count(>(2), values(counts)) == 0
end

@testset "adaptive meshplot: the wireframe is the surface's own edges" begin
    grid = generate_grid(QuadraticQuadrilateral, (4, 4))
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
    close!(dh)
    u = zeros(ndofs(dh))
    Ferrite.apply_analytical!(u, dh, :u,
        x -> Ferrite.Vec(0.18 * sin(pi * x[1]) * sin(pi * x[2]),
                         0.22 * sin(pi * x[1]) * cos(0.5pi * x[2])))
    wds = FEData(dh, u) |> WarpByVector(:u, 1.0)

    # the fixed subdivision does not react to the tolerance; the adaptive one
    # does (the static reference opts out and subdivides uniformly)
    sds = FEData(dh, u; adaptivity=false) |> Refine() |> WarpByVector(:u, 1.0)
    _, _, static = meshplot(sds)
    nstatic = length(static.edge_lines[]) ÷ 2
    counts = map((1.0e-3, 1.0e-4, 1.0e-5)) do tol
        wds.adaptivity.geometry_tol[] = tol
        wds.adaptivity.max_depth[] = 12
        _, _, wp = meshplot(wds)
        length(wp.edge_lines[]) ÷ 2
    end
    @test issorted(counts)
    @test counts[1] < nstatic < counts[3]

    # at one geometry tolerance and without solution refinement, wireframe and
    # surface come from the same key set — so every drawn segment must be an
    # edge of the drawn surface, exactly
    wds.adaptivity.geometry_tol[] = 1.0e-4
    wds.adaptivity.max_depth[] = 12
    _, _, sp = solutionplot(wds; color = :red)
    _, _, wp = meshplot(wds)
    key(p) = (round(Float64(p[1]); digits = 6) + 0.0, round(Float64(p[2]); digits = 6) + 0.0)
    pos, faces = sp.subd_positions[], sp.subd_faces[]
    surface_edges = Set{Any}()
    for f in faces, (i, j) in ((1, 2), (2, 3), (3, 1))
        a, b = key(pos[f[i]]), key(pos[f[j]])
        push!(surface_edges, a <= b ? (a, b) : (b, a))
    end
    segs = wp.edge_lines[]
    nseg = length(segs) ÷ 2
    @test nseg > 0
    @test all(1:nseg) do i
        a, b = key(segs[2i - 1]), key(segs[2i])
        (a <= b ? (a, b) : (b, a)) in surface_edges
    end

    # element edges are drawn once, not once per adjacent cell
    @test nseg < 2 * length(unique(map(i -> key(segs[2i - 1]), 1:nseg)))

    # and it follows the solution like everything else
    before = copy(wp.edge_lines[])
    FerriteViz.update!(wds, 2 .* u)
    @test wp.edge_lines[] != before
end
