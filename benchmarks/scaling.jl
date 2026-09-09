# Scaling comparison across FerriteViz revisions: FEData construction,
# solutionplot and meshplot time-to-first-plot (through `colorbuffer`, so the
# whole pipeline down to the rendered image counts), and steady-state
# `update!`, over mesh sizes in 2D and 3D.
#
#   julia scaling.jl <path-to-FerriteViz> <label> <out.json> [--adaptive]
#
# Run once per revision (e.g. against a master worktree and a branch checkout)
# with identical Makie/GLMakie versions and compare the JSONs. `--adaptive`
# additionally measures the error-adaptive path (branch only): solutionplot
# with the default tolerances and the paired adaptive meshplot.
#
# The 2D scene is a Q2 grid with a Q2 field, the 3D scene Q1 hexes with a Q2
# field — both render curved, so the automatic static subdivision and the
# adaptive estimators have real work to do. Code is warmed up on a tiny grid
# first; the per-size numbers are warm-code/cold-data.

import Pkg
fviz, label, outfile = ARGS[1], ARGS[2], ARGS[3]
adaptive = "--adaptive" in ARGS
Pkg.activate(; temp=true, io=devnull)
Pkg.develop(path=fviz; io=devnull)
Pkg.add([Pkg.PackageSpec(name="GLMakie", version="0.13.13"),
         Pkg.PackageSpec(name="Ferrite"), Pkg.PackageSpec(name="JSON")]; io=devnull)

using GLMakie, Ferrite, FerriteViz, JSON
GLMakie.activate!()

function scene(dim::Int, N::Int)
    if dim == 2
        grid = generate_grid(QuadraticQuadrilateral, (N, N))
        dh = DofHandler(grid)
        add!(dh, :p, Lagrange{RefQuadrilateral,2}())
        close!(dh)
        u = zeros(ndofs(dh))
        Ferrite.apply_analytical!(u, dh, :p, x -> sin(pi * x[1]) * sin(pi * x[2]))
    else
        grid = generate_grid(Hexahedron, (N, N, N))
        dh = DofHandler(grid)
        add!(dh, :p, Lagrange{RefHexahedron,2}())
        close!(dh)
        u = zeros(ndofs(dh))
        # varies on the boundary faces (a product of sines vanishes there,
        # leaving the drawn surface identically zero — nothing to measure)
        Ferrite.apply_analytical!(u, dh, :p, x -> sin(pi * (x[1] + x[2] + x[3]) / 2))
    end
    return dh, u
end

function measure(dim::Int, N::Int; adaptive::Bool)
    dh, u = scene(dim, N)
    t_fed = @elapsed ds = FEData(dh, u)
    m_fed = Base.summarysize(ds)

    sp_kw = adaptive ? (; adaptive=true, solution_tol=5.0e-3, geometry_tol=1.0e-3, max_depth=8) : (;)
    mp_kw = adaptive ? (; adaptive=true, geometry_tol=1.0e-3, max_depth=8) : (;)
    fig = Figure(size=(600, 500))
    ax = dim == 3 ? Axis3(fig[1, 1]) : Axis(fig[1, 1])
    local sp
    t_sp = @elapsed begin
        sp = solutionplot!(ax, ds; sp_kw...)
        colorbuffer(fig)
    end
    t_mp = @elapsed begin
        meshplot!(ax, ds; mp_kw...)
        colorbuffer(fig)
    end

    u2 = 1.05 .* u
    flip = Ref(false)
    step!() = (flip[] = !flip[]; FerriteViz.update!(ds, flip[] ? u2 : u); colorbuffer(fig); nothing)
    step!(); step!()
    t_upd = minimum((@elapsed step!()) for _ in 1:6)

    # drawn triangles (the static path builds every cell's but draws the
    # visible subset — in 3D the boundary cells); built separately
    ntri = adaptive ? length(sp.subd_keys[]) :
           count(i -> ds.visible[ds.triangle_cell_map[i]], eachindex(ds.triangle_cell_map))
    return (; dim, N, label=(adaptive ? "adaptive" : "static"), ndofs=ndofs(dh),
            ncells=Ferrite.getncells(Ferrite.get_grid(dh)), ntri,
            ntri_built=length(ds.all_triangles), t_fed, m_fed, t_sp, t_mp, t_upd)
end

fmt(t) = lpad(round(1000t; digits=1), 9)
rows = Any[]
modes = adaptive ? (false, true) : (false,)
for a in modes            # warm up compile on a tiny scene per mode
    measure(2, 2; adaptive=a)
    measure(3, 2; adaptive=a)
end
for (dim, sizes) in ((2, (16, 32, 64, 128)), (3, (6, 12, 24)))
    for N in sizes, a in modes
        r = measure(dim, N; adaptive=a)
        push!(rows, r)
        println(dim, "D N=", lpad(N, 3), " ", rpad(r.label, 8), ": tri=", lpad(r.ntri, 8),
                "  FEData", fmt(r.t_fed), "ms  solutionplot", fmt(r.t_sp),
                "ms  meshplot", fmt(r.t_mp), "ms  update", fmt(r.t_upd), "ms")
        GC.gc()
    end
end
open(outfile, "w") do io
    JSON.print(io, Dict("label" => label, "rows" => rows))
end
println("wrote ", outfile)
