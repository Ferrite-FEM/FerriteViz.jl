# Does the adaptive path hold up at "real 3D" scale — ~100k FE cells — and
# what do the substrate's shared deviation memos cost and save there?
#
#   julia substrate_scale.jl <path-to-FerriteViz> [ncells_per_side]
#
# The scene: a hex block (default 46³ ≈ 97k cells) with a Q2 scalar field
# carrying a localized feature on one face, drawn adaptively. Only the
# boundary surface enters the adaptive base (the interior is extracted away),
# so the numbers scale with the drawn surface, not the cell count — which is
# the point to verify. Reported:
#
#   * substrate build (base fans, adjacency, geometry coefficients),
#   * the first plot's refinement (fills the deviation memos),
#   * a second identical plot's refinement (should be lookups),
#   * tolerance-only re-decisions (loosen: pure lookups; tighten: only new keys),
#   * a solution update (memos invalidated, full re-evaluation),
#   * memo entry counts and their measured memory.
#
# The dataset is built with `adaptive=false` (flat static tessellation): the
# static Refine() default would tessellate every *volume* cell 4× subdivided,
# which an adaptive-only session never draws — at this scale that is the
# difference between ~2.4M and ~9.6M static triangles kept alive.

import Pkg
fviz_path = ARGS[1]
const NSIDE = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 46
Pkg.activate(; temp=true, io=devnull)
Pkg.develop(path=fviz_path; io=devnull)
Pkg.add([Pkg.PackageSpec(name="Makie", version="0.24.13"),
         Pkg.PackageSpec(name="Ferrite")]; io=devnull)

using Makie, Ferrite, FerriteViz
import FerriteViz as FV

fmt(t) = string(round(1000t; digits=1), " ms")
mib(b) = string(round(b / 2^20; digits=1), " MiB")

println("grid: ", NSIDE, "³ = ", NSIDE^3, " hexes")
t = @elapsed begin
    grid = generate_grid(Hexahedron, (NSIDE, NSIDE, NSIDE))
    dh = DofHandler(grid)
    add!(dh, :p, Lagrange{RefHexahedron,2}())
    close!(dh)
end
println("grid + dofs (", ndofs(dh), " dofs): ", fmt(t))

u = zeros(ndofs(dh))
Ferrite.apply_analytical!(u, dh, :p,
                          x -> exp(-60 * ((x[1] - 0.3)^2 + (x[2] + 0.2)^2 + (x[3] - 1.0)^2)))
u2 = 2 .* u

t = @elapsed ds = FEData(dh, u; adaptive=false)
println("FEData (static tessellation, flat): ", fmt(t))
t = @elapsed sub = FV._substrate(ds)
println("substrate (", length(sub.base.corners), " base triangles): ", fmt(t))

kw = (; adaptive=true, solution_tol=5e-3, max_depth=8)
decide(sp) = (sp.subd_keys[]; sp.subd_positions[]; sp.subd_color[]; length(sp.subd_keys[]))

t1 = @elapsed begin
    _, _, sp1 = solutionplot(ds; kw...)
    n = decide(sp1)
end
println("plot 1, first decision (", n, " triangles): ", fmt(t1))
t2 = @elapsed begin
    _, _, sp2 = solutionplot(ds; kw...)
    decide(sp2)
end
println("plot 2, same criterion:              ", fmt(t2), "   (", round(t1 / t2; digits=1), "x)")

tl = @elapsed (Makie.update!(sp2, solution_tol=5e-2); decide(sp2))
println("tolerance loosened (", length(sp2.subd_keys[]), " triangles):  ", fmt(tl))
tt = @elapsed (Makie.update!(sp2, solution_tol=5e-4); decide(sp2))
println("tolerance tightened (", length(sp2.subd_keys[]), " triangles): ", fmt(tt))
Makie.update!(sp2, solution_tol=5e-3); decide(sp2)

tu = @elapsed (FerriteViz.update!(ds, u2); decide(sp1); decide(sp2))
println("solution update, both plots:         ", fmt(tu))

for (term, c) in sub.dev_caches
    println("memo :", term, ": ", length(c), " entries, ", mib(Base.summarysize(c)))
end
println("substrate total: ", mib(Base.summarysize(sub)))
