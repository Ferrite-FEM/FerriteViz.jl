# Advanced Topics

## Live plotting

Plotting while a computational heavy simulation is performed can be easily achieved with FerriteVis.jl.
Every plotter object of type `MakiePlotter` holds a property called `u` which is a so called `Observable`.
If an `Observable` changes, all its dependencies are triggered to change as well. So, all we need to do is to update
the observable `plotter.u`. For this purpose the function [`FerriteVis.update!`](@ref) is provided. It takes a `plotter:MakiePlotter`
and a new solutiuon vector `u_new` and updates `plotter.u`, thereby all open plots called with `plotter` are updated.
