# The plot owns a live cut of the continuous domain; FEData's arrays are the
# coarse snapshot used by static filters and data inspection.
_render_area(m) = sum(FerriteViz._tri_area(_vv3(m.positions[convert(Int,f[1])]),
                                          _vv3(m.positions[convert(Int,f[2])]),
                                          _vv3(m.positions[convert(Int,f[3])])) for f in m.faces; init=0.0)

function _render_closed_volume(m)
    edges = Dict{Tuple{NTuple{3,Float64},NTuple{3,Float64}},Int}()
    volume = 0.0
    key(p) = ntuple(i -> round(Float64(p[i]); digits=5) + 0.0, 3)
    for face in m.faces
        ps = ntuple(i -> _vv3(m.positions[convert(Int,face[i])]), 3)
        volume += ps[1] ⋅ (ps[2] × ps[3]) / 6
        for (i,j) in ((1,2),(2,3),(3,1))
            a,b = key(ps[i]),key(ps[j])
            edge = isless(a,b) ? (a,b) : (b,a)
            edges[edge] = get(edges,edge,0)+1
        end
    end
    return all(==(2),values(edges)), volume
end

@testset "cut rendering: planes and pruning" begin
    ds,_,_ = _hexds(3; f=x -> x[1]+2x[2]-x[3])
    for offset in (0.1, 1/3, -1.0, 1.0, 4.0)
        cut = ds |> Clip(ClipPlane(Vec(1.,0.,0.),offset))
        st = FerriteViz.CutRenderState(cut,:u)
        m = FerriteViz._render_cut(st,nothing)
        expected = offset <= -1 ? 0.0 : offset >= 1 ? 24.0 : 8+8*(offset+1)
        @test _render_area(m) ≈ expected atol=2e-5
        @test all(p -> p[1] <= offset+1e-6,m.positions)
    end
    slab = ds |> Clip(ClipPlane(Vec(1.,0.,0.),0.4)) |>
                 Clip(ClipPlane(Vec(-1.,0.,0.),0.0))
    m = FerriteViz._render_cut(FerriteViz.CutRenderState(slab,:u),nothing)
    @test _render_area(m) ≈ 11.2 atol=2e-5
    same = slab |> Clip(ClipPlane(Vec(2.,0.,0.),0.8))
    @test _render_area(FerriteViz._render_cut(FerriteViz.CutRenderState(same,:u),nothing)) ≈ 11.2 atol=2e-5

    source = ds |> Refine(1)
    refined_cut = source |> Clip(ClipPlane(Vec(1.,0.,0.),0.1))
    coarse_cut = ds |> Clip(ClipPlane(Vec(1.,0.,0.),0.1))
    mref = FerriteViz._render_cut(FerriteViz.CutRenderState(refined_cut,:u),nothing)
    mbase = FerriteViz._render_cut(FerriteViz.CutRenderState(coarse_cut,:u),nothing)
    @test length(mref.faces) > length(mbase.faces)
    @test _render_area(mref) ≈ _render_area(mbase) atol=2e-5
    @test _render_closed_volume(mref)[1]

    quadratic,_,_ = _hexds(2; ip=Lagrange{RefHexahedron,2}(), f=x -> x[2]^2)
    curved_colors = quadratic |> Clip(ClipPlane(Vec(1.,0.,0.),0.1))
    cfg = Adaptivity(solution_tol=0.01,max_depth=6)
    refined = FerriteViz._render_cut(FerriteViz.CutRenderState(curved_colors,:u),cfg)
    closed, volume = _render_closed_volume(refined)
    @test closed
    @test volume ≈ 4.4 atol=2e-5

    large,_,_ = _hexds(12)
    cut = large |> Clip(ClipPlane(Vec(1.,0.,0.),0.1))
    m = FerriteViz._render_cut(FerriteViz.CutRenderState(cut,:u),nothing)
    @test m.candidates < getncells(Ferrite.get_grid(large.dh)) ÷ 3
end

@testset "cut rendering: nonlinear interior and compute graph" begin
    grid = generate_grid(Hexahedron,(1,1,1))
    dh = DofHandler(grid)
    ip = Lagrange{RefHexahedron,2}()
    add!(dh,:u,ip^3); add!(dh,:p,ip); close!(dh)
    u = zeros(ndofs(dh))
    warp(x) = Vec(0.3*(1-x[1]^2)*(1-x[2]^2)*(1-x[3]^2),0.0,0.0)
    affine(x) = x ⋅ Vec(1.,2.,3.)
    Ferrite.apply_analytical!(u,dh,:u,warp)
    Ferrite.apply_analytical!(u,dh,:p,x -> affine(x+warp(x)))
    cfg = Adaptivity(geometry_tol=0.02,solution_tol=0.02,max_depth=6)
    ds = FEData(dh,u;adaptivity=cfg)
    scale = Makie.Observable(1.0)
    cut = ds |> WarpByVector(:u,scale) |> Clip(ClipPlane(Vec(1.,0.,0.),0.1))
    st = FerriteViz.CutRenderState(cut,:p)
    coarse = FerriteViz._render_cut(st,nothing)
    fine = FerriteViz._render_cut(st,cfg)
    @test fine.depth > 0
    @test length(fine.faces) > length(coarse.faces)
    @test maximum(abs(fine.colors[i]-affine(_vv3(fine.positions[i]))) for i in eachindex(fine.colors)) < 2e-6
    @test all(p -> p[1] <= 0.100001,fine.positions)
    @test _render_area(fine) ≈ 16.8 atol=2e-5
    closed, volume = _render_closed_volume(fine)
    @test closed
    @test volume ≈ 4.4 atol=2e-5

    _,_,sp = solutionplot(cut;color=:p)
    _,_,wf = meshplot(cut)
    before = copy(sp.cut_positions[])
    for amplitude in (1.5,0.0,1.0)
        scale[] = amplitude
        frame = copy(u)
        Ferrite.apply_analytical!(frame,dh,:p,x -> affine(x + amplitude*warp(x)))
        FerriteViz.update!(ds,frame)
        @test maximum(abs(sp.cut_colors[][i]-affine(_vv3(sp.cut_positions[][i])))
                      for i in eachindex(sp.cut_colors[])) < 2e-6
        @test all(p -> p[1] <= 0.100001,sp.cut_positions[])
        @test length(sp.cut_colors[]) == length(sp.cut_positions[])
        @test all(p -> p[1] <= 0.100001,wf.edge_lines[])
    end
    FerriteViz.update!(ds,2u)
    @test sp.cut_positions[] != before
    @test all(p -> p[1] <= 0.100001,sp.cut_positions[])
    cfg.max_depth[] = 0
    @test length(sp.cut_faces[]) == length(coarse.faces)
    cfg.max_depth[] = 6
    @test length(sp.cut_faces[]) > length(coarse.faces)

    # A cell can enter the cut later even if the apply-time snapshot is empty.
    translated = zeros(ndofs(dh))
    Ferrite.apply_analytical!(translated,dh,:u,x -> Vec(4.,0.,0.))
    source = FEData(dh,translated;adaptivity=false)
    moving = source |> WarpByVector(:u) |> Clip(ClipPlane(Vec(1.,0.,0.),0.0))
    @test isempty(moving.all_triangles)
    _,_,plot = solutionplot(moving;color=:p)
    @test isempty(plot.cut_faces[])
    FerriteViz.update!(source,zeros(ndofs(dh)))
    @test !isempty(plot.cut_faces[])
    @test all(p -> p[1] <= 1e-6,plot.cut_positions[])
end

@testset "cut bounds include nonnodal extrema" begin
    # x'=x+2(1-y²) overshoots the original cell AABB. The right-hand
    # portion is still found when all geometric corner nodes miss the plane.
    grid = generate_grid(Hexahedron,(1,1,1))
    dh = DofHandler(grid);add!(dh,:u,Lagrange{RefHexahedron,2}()^3);close!(dh)
    u=zeros(ndofs(dh))
    Ferrite.apply_analytical!(u,dh,:u,x -> Vec(2*(1-x[2]^2),0.,0.))
    ds=FEData(dh,u;adaptivity=Adaptivity(geometry_tol=0.03,max_depth=6))
    cut=ds |> WarpByVector(:u) |> Clip(ClipPlane(Vec(-1.,0.,0.),-2.5))
    m=FerriteViz._render_cut(FerriteViz.CutRenderState(cut,:black),ds.adaptivity)
    @test !isempty(m.faces)
    @test all(p -> p[1] >= 2.5-1e-6,m.positions)
end
