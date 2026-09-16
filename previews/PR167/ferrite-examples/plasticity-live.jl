# Live-plotting variant of the plasticity example: the viewer pops up before
# time stepping begins and updates after every converged step via
# `FerriteViz.update!` (see the live plotting section of the docs).
include("plasticity.jl")

u, dh, u_history, σ, κ = solve(true)
