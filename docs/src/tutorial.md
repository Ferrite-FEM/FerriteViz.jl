# Tutorial

## Solve a Boundary Value Problem

Start with solving a boundary value problem as you would usually do with Ferrite. It is crucial that you safe your used DofHandler
and solution vector because we need to pass those objects to `MakiePlotter`.


## Plot your results

!!! tip "Plotting Functions"
    Currently, [`FerriteViz.solutionplot`](@ref), [`FerriteViz.wireframe`](@ref), [`FerriteViz.surface`](@ref), [`FerriteViz.arrows`](@ref) and their mutating analogues with `!` are defined for `MakiePlotter`.
    Due to the nature of the documentation we need `WGLMakie`, however, you can simply exchange any `WGLMakie` call by `GLMakie`.

```@example 1
import JSServe # hide
JSServe.Page(exportable=true, offline=true) # hide
```

You can start by plotting your mesh

```@example 1
import FerriteViz
using Ferrite
import WGLMakie #activating the backend, switch to GLMakie or CairoMakie (for 2D) locally
WGLMakie.set_theme!(resolution=(800, 400)) # hide

grid = generate_grid(Hexahedron,(3,3,3))
FerriteViz.wireframe(grid,markersize=50,strokewidth=2)
```

FerriteViz.jl also supports showing labels for `Ferrite.AbstractGrid` entities, such as node- and celllabels, as well as plotting cellsets.

```@example 1
grid = generate_grid(Quadrilateral,(3,3))
addcellset!(grid,"s1",Set((1,4,7)))
addcellset!(grid,"s2",Set((2,5,8)))
addcellset!(grid,"s3",Set((3,6,9)))
FerriteViz.wireframe(grid,markersize=5,strokewidth=1,nodelabels=true,celllabels=true,cellsets=true)
```

If you solve some boundary value problem with Ferrite.jl keep in mind to safe your `dh::DofHandler` and solution vector `u::Vector{T}` in some variable.
With them, we create the `MakiePlotter` struct that dispatches on the plotting functions.

```@example 1
include("ferrite-examples/incompressible-elasticity.jl") #defines variables dh and u

plotter = FerriteViz.MakiePlotter(dh,u)
FerriteViz.arrows(plotter)
```

Per default, all plotting functions grab the first field in the `DofHandler`, but of course you can plot a different field as well.
The next plot will show the pressure instead of the displacement

```@example 1
FerriteViz.solutionplot(plotter,field=:p)
```

For certain 2D problems it makes sense to visualize the result as a `surface` plot. To showcase the combination with the mutating versions of the plotting functions,
the `solutionplot` function is plotted below the `surface` plot

```@example 1
FerriteViz.surface(plotter)
FerriteViz.solutionplot!(plotter,colormap=:magma)
WGLMakie.current_figure()
```

However, in structural mechanics we often would like to see the deformed configuration,
which can be achieved by providing a `deformation_field::Symbol` as a keyword argument.

```@example 1
include("ferrite-examples/plasticity.jl") #only defines solving function
u, dh, uhistory, ??, ?? = solve()
plotter = FerriteViz.MakiePlotter(dh,u)

FerriteViz.solutionplot(plotter,colormap=:thermal,deformation_field=:u)
FerriteViz.wireframe!(plotter,deformation_field=:u,markersize=25,strokewidth=1)
WGLMakie.current_figure()
```

FerriteViz.jl also supports to plot cell data, such as the **averaged** von-Mises stress or the drag stress of the plasticity example.
```@example 1
FerriteViz.cellplot(plotter,??,colormap=:thermal,deformation_field=:u,deformation_scale=2.0)
FerriteViz.wireframe!(plotter,deformation_field=:u,markersize=25,strokewidth=1,deformation_scale=2.0)
WGLMakie.current_figure()
```

Further, this package provides an interactive viewer that you can call with `ferriteviewer(plotter)` and
`ferriteviewer(plotter,u_history)` for time dependent views, respectively.
If you want to live plot your solution while solving some finite element system, consider to take a look at the advanced topics page.
