# Tutorial

## Solve a Boundary Value Problem

Start with solving a boundary value problem as you would usually do with Ferrite. It is crucial that you safe your used DofHandler
and solution vector because we need to pass those objects to `MakiePlotter`.


## Plot your results

Currently, `surface`, `mesh`, `arrows` and their mutating analogues with `!` are defined for `MakiePlotter`

```@example
include("ferrite-examples/incompressible-elasticity.jl")

using FerriteVis
import GLMakie
plotter = MakiePlotter(dh,u)

GLMakie.arrows(plotter)
```

Per default, all plotting functions grab the first field in the `DofHandler`, but of course you can plot a different field as well.
The next plot will show the pressure instead of the displacement

```@example 
include("ferrite-examples/incompressible-elasticity.jl")

using FerriteVis
import GLMakie
plotter = MakiePlotter(dh,u)

GLMakie.mesh(plotter,field=2)
```

For certain 2D problems it makes sense to visualize the result as a `surface` plot. To showcase the combination with the mutating versions of the plotting functions,
the `mesh` function is plotted below the `surface` plot

```@example 
include("ferrite-examples/incompressible-elasticity.jl")

using FerriteVis
import GLMakie
plotter = MakiePlotter(dh,u)

GLMakie.surface(plotter)
GLMakie.mesh!(plotter)
GLMakie.current_figure()
```

However, in structural mechanics we often would like to see the deformed configuration,
so there is a custom function `warp_by_vector`, which does the same as `warp by glyph` in Paraview.

```@example
include("ferrite-examples/incompressible-elasticity.jl")

using FerriteVis
import GLMakie
plotter = MakiePlotter(dh,u)

warp_by_vector(plotter)
```
