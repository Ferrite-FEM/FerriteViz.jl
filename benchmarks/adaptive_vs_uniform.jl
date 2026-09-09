# Accuracy-matched comparison: uniform `Refine(n)` against error-adaptive
# tessellation, on a curved+warped nonlinear example.
#
#   julia adaptive_vs_uniform.jl <path-to-FerriteViz> static  out.json
#   julia adaptive_vs_uniform.jl <path-to-FerriteViz> matched out.json <static.json>
#
# `static` sweeps uniform refinement levels (run it against a checkout without
# the adaptive path to compare releases); `matched` tunes the two adaptive
# tolerances until both measured errors meet each level's, and reports the
# cheapest mesh that does. SCENARIO=smooth|local picks a globally smooth field
# or a localized bump.
#
# Both modes measure the SAME thing, independent of how the mesh was made: the
# drawn surface is piecewise linear over the triangles, so its error against
# the exact geometry / exact field is sampled inside every drawn triangle by
# comparing the linear interpolation of the triangle's vertex data against the
# exact value at the corresponding reference point.

import Pkg
fviz_path, mode, outfile = ARGS[1], ARGS[2], ARGS[3]
Pkg.activate(; temp=true, io=devnull)
Pkg.develop(path=fviz_path; io=devnull)
Pkg.add([Pkg.PackageSpec(name="Makie", version="0.24.13"),
         Pkg.PackageSpec(name="Ferrite"), Pkg.PackageSpec(name="JSON")]; io=devnull)

using Makie, Ferrite, FerriteViz, LinearAlgebra, JSON
import FerriteViz as FV

############
# The case #
############
# Curved + warped: a Q2 grid displaced by an analytic sine field (the geometry
# nonlinearity) and coloured by a second analytic field (the solution
# nonlinearity). Both are resolved exactly by the Q2 spaces, so the *only*
# error left is the piecewise-linear rendering — which is what we measure.
const N = 6
const WARP = 1.0

grid = generate_grid(QuadraticQuadrilateral, (N, N))
dh = DofHandler(grid)
add!(dh, :u, Lagrange{RefQuadrilateral,2}()^2)
add!(dh, :s, Lagrange{RefQuadrilateral,2}())
close!(dh)
const SCENARIO = get(ENV, "SCENARIO", "smooth")
u = zeros(ndofs(dh))
if SCENARIO == "smooth"
    # globally nonlinear: sine waves of the same amplitude everywhere, so the
    # rendering error is spread evenly and uniform refinement is near-optimal
    Ferrite.apply_analytical!(u, dh, :u,
        x -> Ferrite.Vec(0.12 * sin(pi * x[1]) * sin(pi * x[2]), 0.10 * sin(2pi * x[1]) * x[2]))
    Ferrite.apply_analytical!(u, dh, :s, x -> sin(pi * x[1]) * cos(pi * x[2]))
else
    # localized: one bump in the geometry and one in the field, flat elsewhere
    Ferrite.apply_analytical!(u, dh, :u,
        x -> Ferrite.Vec(0.0, 0.45 * exp(-25 * ((x[1] - 0.2)^2 + (x[2] + 0.1)^2))))
    Ferrite.apply_analytical!(u, dh, :s, x -> exp(-25 * ((x[1] - 0.2)^2 + (x[2] + 0.1)^2)))
end
u2 = 1.3 .* u      # a second state, for the update-path timing

#####################
# Exact evaluators  #
#####################
struct Probe{PV}
    pv::PV
    dofs::Vector{Vector{Int}}
end
function Probe(dh, field::Symbol)
    sdh = first(dh.subdofhandlers)
    ip = Ferrite.getfieldinterpolation(sdh, field)
    gip = Ferrite.geometric_interpolation(Ferrite.getcelltype(sdh))
    pv = Ferrite.PointValues(ip, gip; update_gradients=false)
    rng = Ferrite.dof_range(sdh, field)
    dofs = [Ferrite.celldofs(dh, c)[rng] for c in 1:Ferrite.getncells(Ferrite.get_grid(dh))]
    return Probe(pv, dofs)
end
function (p::Probe)(cell::Int, ξ, coords, uvec)
    Ferrite.reinit!(p.pv, coords, ξ)
    return Ferrite.function_value(p.pv, 1, @views(uvec[p.dofs[cell]]))
end

const grid_ = Ferrite.get_grid(dh)
const cellcoords = [Ferrite.getcoordinates(grid_, c) for c in 1:Ferrite.getncells(grid_)]
const gips = [Ferrite.geometric_interpolation(typeof(c)) for c in Ferrite.getcells(grid_)]
const probe_u = Probe(dh, :u)
const probe_s = Probe(dh, :s)

# exact drawn geometry: the cell's geometric map plus the warp displacement
exact_x(cell, ξ, uvec) =
    FV.geometric_map(gips[cell], cellcoords[cell], ξ) + WARP * probe_u(cell, ξ, cellcoords[cell], uvec)
exact_s(cell, ξ, uvec) = probe_s(cell, ξ, cellcoords[cell], uvec)

# Barycentric sample points strictly inside a triangle, where linear
# interpolation is worst.
const SAMPLES = [(1/3, 1/3, 1/3), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5),
                 (2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 1/6, 2/3),
                 (0.25, 0.25, 0.5), (0.5, 0.25, 0.25), (0.25, 0.5, 0.25)]

function measure(positions, refcoords, faces, cell_of, colors, uvec)
    eg = es = 0.0
    for (t, f) in enumerate(faces)
        c = cell_of[t]
        ξ1, ξ2, ξ3 = refcoords[f[1]], refcoords[f[2]], refcoords[f[3]]
        x1, x2, x3 = positions[f[1]], positions[f[2]], positions[f[3]]
        s1, s2, s3 = colors[f[1]], colors[f[2]], colors[f[3]]
        for (a, b, c3) in SAMPLES
            ξ = a * ξ1 + b * ξ2 + c3 * ξ3
            xlin = a * Float64.(x1) + b * Float64.(x2) + c3 * Float64.(x3)
            slin = a * s1 + b * s2 + c3 * s3
            eg = max(eg, norm(collect(exact_x(c, ξ, uvec)) - collect(xlin)))
            es = max(es, abs(exact_s(c, ξ, uvec) - slin))
        end
    end
    return eg, es
end

# minimum over repeats: the timings are milliseconds, so one stray GC pause or
# a late specialization would otherwise dominate
function best(f, n=5)
    t = Inf
    for _ in 1:n; t = min(t, @elapsed f()); end
    return t
end

payload(positions, faces, colors) =
    Base.summarysize(positions) + Base.summarysize(faces) + Base.summarysize(colors)

#####################
# The two pipelines #
#####################
Vec2 = Ferrite.Vec{2,Float64}

function run_static(n::Int)
    build() = begin
        ds = FEData(dh, copy(u); adaptive=false) |> Refine(n) |> WarpByVector(:u, WARP)
        fig = solutionplot(ds; color=:s)
        (ds, fig)
    end
    build()                                   # warm up
    t_build = best(build, 3)
    (ds, fig) = build()
    flip = Ref(false)
    step!() = (flip[] = !flip[]; FerriteViz.update!(ds, flip[] ? u2 : u); ds.coords[]; FV.point_data(ds, :s)[])
    step!(); step!()
    t_update = best(step!)
    m_update = @allocated step!()
    while ds.u[] != u; step!(); end

    positions = ds.coords[]
    faces = ds.all_triangles
    refcoords = [Vec2((ds.reference_coords[v, 1], ds.reference_coords[v, 2])) for v in 1:size(ds.reference_coords, 1)]
    colors = vec(FV.point_data(ds, :s)[])
    eg, es = measure(positions, refcoords, faces, ds.triangle_cell_map, colors, u)
    return Dict("label" => "Refine($n)", "param" => n, "ntri" => length(faces),
                "nvert" => length(positions), "geom_err" => eg, "sol_err" => es,
                "t_build" => t_build, "t_update" => t_update, "mem_update" => m_update,
                "payload" => payload(positions, faces, colors))
end

function run_adaptive(gtol::Float64, stol::Float64=gtol)
    build() = begin
        ds = FEData(dh, copy(u); adaptive=false)
        wds = ds |> WarpByVector(:u, WARP)
        fig, ax, sp = solutionplot(wds; color=:s, adaptive=true,
                                   geometry_tol=gtol, solution_tol=stol, max_depth=12)
        (ds, sp)
    end
    build()
    t_build = best(build, 3)
    (ds, sp) = build()
    flip = Ref(false)
    step!() = (flip[] = !flip[]; FerriteViz.update!(ds, flip[] ? u2 : u); sp.subd_positions[]; sp.subd_color[])
    step!(); step!()
    t_update = best(step!)
    m_update = @allocated step!()
    while ds.u[] != u; step!(); end

    cellmap = FV._substrate(ds |> WarpByVector(:u, WARP)).cellmap
    keys = sp.subd_keys[]
    cell_of = [cellmap[FV.key_base(k)] for k in keys]
    positions, faces, colors = sp.subd_positions[], sp.subd_faces[], sp.subd_color[]
    refcoords = sp.subd_ξ[]
    eg, es = measure(positions, refcoords, faces, cell_of, colors, u)
    return Dict("label" => "adaptive g=$gtol s=$stol", "geometry_tol" => gtol, "solution_tol" => stol, "ntri" => length(faces),
                "nvert" => length(positions), "geom_err" => eg, "sol_err" => es,
                "t_build" => t_build, "t_update" => t_update, "mem_update" => m_update,
                "payload" => payload(positions, faces, colors),
                "max_depth_used" => maximum(FV.key_depth, keys))
end

results = Dict{String,Any}[]
if mode == "static"
    for n in 0:4
        r = run_static(n)
        println(r["label"], ": tri=", r["ntri"], " eg=", r["geom_err"], " es=", r["sol_err"],
                " build=", round(1000r["t_build"], digits=1), "ms upd=", round(1000r["t_update"], digits=2), "ms")
        push!(results, r)
    end
elseif mode == "adaptive"
    for tol in (3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
        r = run_adaptive(tol)
        println(r["label"], ": tri=", r["ntri"], " eg=", r["geom_err"], " es=", r["sol_err"],
                " build=", round(1000r["t_build"], digits=1), "ms upd=", round(1000r["t_update"], digits=2), "ms")
        push!(results, r)
    end
else
    # Accuracy-matched: for every uniform level, tune the two tolerances
    # independently until both measured errors meet that level's, and keep the
    # cheapest mesh that does. The estimators bound the error nearly linearly
    # in their tolerance, so scaling by target/measured converges in a few
    # rounds; the clamp keeps a bad first guess from overshooting.
    targets = JSON.parsefile(ARGS[4])["results"]
    for t in targets
        egt, est = t["geom_err"], t["sol_err"]
        gtol, stol = egt / 3, est / 2          # first guess: tolerances are relative
        best = nothing
        for iter in 1:6
            r = run_adaptive(gtol, stol)
            if r["geom_err"] <= egt && r["sol_err"] <= est &&
               (best === nothing || r["ntri"] < best["ntri"])
                best = r
            end
            gtol *= clamp(egt / max(r["geom_err"], 1e-15), 0.3, 3.0)
            stol *= clamp(est / max(r["sol_err"], 1e-15), 0.3, 3.0)
        end
        if best === nothing
            println("target ", t["label"], ": NOT MATCHED")
        else
            best["target"] = t["label"]
            best["target_ntri"] = t["ntri"]
            println("target ", t["label"], " (", t["ntri"], " tri): adaptive ", best["ntri"],
                    " tri (", round(best["ntri"] / t["ntri"], digits=2), "x)  eg=", best["geom_err"],
                    " es=", best["sol_err"], " build=", round(1000best["t_build"], digits=1),
                    "ms upd=", round(1000best["t_update"], digits=2), "ms")
            push!(results, best)
        end
    end
end
open(outfile, "w") do io
    JSON.print(io, Dict("mode" => mode, "results" => results))
end
println("saved ", outfile)
