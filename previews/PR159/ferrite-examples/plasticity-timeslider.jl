# Time-slider variant of the plasticity example: solve first, then step
# through the stored solution history interactively.
include("plasticity.jl")

u, dh, u_history, σ, κ = solve()
ds = FEData(dh, u)
fig = ferriteviewer(ds, u_history)
display(fig)
