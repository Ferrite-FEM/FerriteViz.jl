# Tutorial

## Solve a Boundary Value Problem

Start with solving a boundary value problem as you would usually do with Ferrite. It is crucial that you safe your used DofHandler
and solution vector because we need to pass those objects to `MakiePlotter`.


## Plot your results

!!! tip "Plotting Functions"
    Currently, [`FerriteVis.solutionplot`](@ref), [`FerriteVis.wireframe`](@ref), [`FerriteVis.surface`](@ref), [`FerriteVis.arrows`](@ref) and their mutating analogues with `!` are defined for `MakiePlotter`.
    Due to the nature of the documentation we need `WGLMakie`, however, you can simply exchange any `WGLMakie` call by `GLMakie`.

```@example 1
import JSServe # hide
JSServe.Page(exportable=true, offline=true) # hide
```

You can start by plotting your mesh

```@example 1
import FerriteVis
using Ferrite
import WGLMakie #activating the backend, switch to GLMakie or CairoMakie (for 2D) locally
WGLMakie.set_theme!(resolution=(800, 400)) # hide

grid = generate_grid(Hexahedron,(3,3,3))
FerriteVis.wireframe(grid,markersize=50,strokewidth=2)
```

FerriteVis.jl also supports showing labels for `Ferrite.AbstractGrid` entities, such as node- and celllabels.

```@example 1
grid = generate_grid(Quadrilateral,(3,3))
FerriteVis.wireframe(grid,markersize=5,strokewidth=1,nodelabels=true,celllabels=true)
```

If you solve some boundary value problem with Ferrite.jl keep in mind to safe your `dh::DofHandler` and solution vector `u::Vector{T}` in some variable.
With them, we create the `MakiePlotter` struct that dispatches on the plotting functions.

```@example 1
include("ferrite-examples/incompressible-elasticity.jl") #defines variables dh and u

plotter = FerriteVis.MakiePlotter(dh,u)
FerriteVis.arrows(plotter)
```

Per default, all plotting functions grab the first field in the `DofHandler`, but of course you can plot a different field as well.
The next plot will show the pressure instead of the displacement

```@example 1
FerriteVis.solutionplot(plotter,field=:p)
```

For certain 2D problems it makes sense to visualize the result as a `surface` plot. To showcase the combination with the mutating versions of the plotting functions,
the `solutionplot` function is plotted below the `surface` plot

```@example 1
FerriteVis.surface(plotter)
FerriteVis.solutionplot!(plotter,colormap=:magma)
WGLMakie.current_figure()
```

However, in structural mechanics we often would like to see the deformed configuration,
which can be achieved by providing a `deformation_field::Symbol` as a keyword argument.

```@example 1
include("ferrite-examples/plasticity.jl") #only defines solving function
u, dh, uhistory, σ, κ = solve()
plotter = FerriteVis.MakiePlotter(dh,u)

FerriteVis.solutionplot(plotter,colormap=:thermal,deformation_field=:u)
FerriteVis.wireframe!(plotter,deformation_field=:u,markersize=25,strokewidth=1)
WGLMakie.current_figure()
```

FerriteVis.jl also supports to plot quadrature data, such as the von-Mises stress or the drag stress of the plasticity example.
```@example 1
FerriteVis.quadratureplot(plotter,σ,colormap=:thermal,deformation_field=:u,deformation_scale=2.0)
FerriteVis.wireframe!(plotter,deformation_field=:u,markersize=25,strokewidth=1)
WGLMakie.current_figure()
```

Further, this package provides an interactive viewer that you can call with `ferriteviewer(plotter)` and
`ferriteviewer(plotter,u_history)` for time dependent views, respectively.
If you want to live plot your solution while solving some finite element system, consider to take a look at the advanced topics page.
