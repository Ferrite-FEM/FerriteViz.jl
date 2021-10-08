# Tutorial

## Solve a Boundary Value Problem

Start with solving a boundary value problem as you would usually do with Ferrite. It is crucial that you safe your used DofHandler
and solution vector because we need to pass those objects to `MakiePlotter`.


## Plot your results

Currently, `surface`, `mesh`, `arrows` and their mutating analogues with `!` are defined for `MakiePlotter`.
Due to the nature of the documentation we need `WGLMakie`, however, you can simply exchange any `WGLMakie` call by `GLMakie`.

```@example 1
import JSServe # hide
JSServe.Page(exportable=true, offline=true) # hide
```

```@example 1
include("ferrite-examples/incompressible-elasticity.jl")

using FerriteVis
import WGLMakie
WGLMakie.set_theme!(resolution=(800, 400)) # hide
plotter = MakiePlotter(dh,u)

WGLMakie.arrows(plotter)
```

Per default, all plotting functions grab the first field in the `DofHandler`, but of course you can plot a different field as well.
The next plot will show the pressure instead of the displacement

```@example 1
WGLMakie.mesh(plotter,field=2)
```

For certain 2D problems it makes sense to visualize the result as a `surface` plot. To showcase the combination with the mutating versions of the plotting functions,
the `mesh` function is plotted below the `surface` plot

```@example 1
WGLMakie.surface(plotter)
WGLMakie.mesh!(plotter)
WGLMakie.current_figure()
```

However, in structural mechanics we often would like to see the deformed configuration,
so there is a custom function `warp_by_vector`, which does the same as `warp by glyph` in Paraview.

```@example 1
warp_by_vector(plotter)
```
