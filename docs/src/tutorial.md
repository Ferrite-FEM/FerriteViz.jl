# Tutorial

## Solve a Boundary Value Problem

Start with solving a boundary value problem as you would usually do with Ferrite. It is crucial that you safe your used DofHandler
and solution vector because we need to pass those objects to `MakiePlotter`.


## Basics

!!! tip "Plotting Functions"
    Currently, [`FerriteViz.solutionplot`](@ref), [`FerriteViz.wireframe`](@ref), [`FerriteViz.surface`](@ref), [`FerriteViz.arrows`](@ref) and their mutating analogues with `!` are defined for `MakiePlotter`.
    Due to the nature of the documentation we need `WGLMakie`, however, you can simply exchange any `WGLMakie` call by `GLMakie`.

### Mesh utilities

```@example 1
import JSServe, WGLMakie # hide
JSServe.Page(exportable=true, offline=true) # hide
WGLMakie.activate!() # hide
WGLMakie.Makie.inline!(true) # hide
```

You can start by plotting your mesh

```@example 1
import FerriteViz
using Ferrite
import WGLMakie #activating the backend, switch to GLMakie or CairoMakie (for 2D) locally
WGLMakie.set_theme!(resolution=(800, 400)) # hide

grid = generate_grid(Hexahedron,(3,3,3))
FerriteViz.wireframe(grid,markersize=10,strokewidth=2)
```

FerriteViz.jl also supports showing labels for `Ferrite.AbstractGrid` entities, such as node- and celllabels, as well as plotting cellsets.

```@example 1
grid = generate_grid(Quadrilateral,(3,3))
addcellset!(grid,"s1",Set((1,4,7)))
addcellset!(grid,"s2",Set((2,5,8)))
addcellset!(grid,"s3",Set((3,6,9)))
FerriteViz.wireframe(grid,markersize=10,strokewidth=1,nodelabels=true,celllabels=true,cellsets=true)
```

### Solution field of a boundary value problem

If you solve some boundary value problem with Ferrite.jl keep in mind to safe your `dh::DofHandler` and solution vector `u::Vector{T}` in some variable.
With them, we create the `MakiePlotter` struct that dispatches on the plotting functions.

```@example 1
include("ferrite-examples/incompressible-elasticity.jl") #defines variables dh_quadratic and u_quadratic

plotter = FerriteViz.MakiePlotter(dh_quadratic,u_quadratic)
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

### Deformed mesh for mechanical boundary value problem

However, in structural mechanics we often would like to see the deformed configuration,
which can be achieved by providing a `deformation_field::Symbol` as a keyword argument.

```@example 1
include("ferrite-examples/plasticity.jl") #only defines solving function
u, dh, uhistory, σ, κ = solve()
plotter = FerriteViz.MakiePlotter(dh,u)

FerriteViz.solutionplot(plotter,colormap=:thermal,deformation_field=:u)
FerriteViz.wireframe!(plotter,deformation_field=:u,markersize=10,strokewidth=1)
WGLMakie.current_figure()
```

### Showing per-cell data

FerriteViz.jl also supports to plot cell data, such as the **averaged** von-Mises stress or the drag stress of the plasticity example.
```@example 1
u, dh, uhistory, σ, κ = solve()
FerriteViz.cellplot(plotter,σ,colormap=:thermal,deformation_field=:u,deformation_scale=2.0)
FerriteViz.wireframe!(plotter,deformation_field=:u,markersize=10,strokewidth=1,deformation_scale=2.0)
WGLMakie.current_figure()
```
For a more granular investigation of the stress field consult the advanced tutorial.

### Interior of a 3D domain

For 3D problems we can also inspect the interior of the domain. Currenly we only have crinkle clipping
implemented and it can be used as follows:
```@example 1
clip_plane = FerriteViz.ClipPlane(Vec((0.0,0.5,0.5)), 0.7)
clipped_plotter = FerriteViz.crinkle_clip(plotter, clip_plane)
FerriteViz.solutionplot(clipped_plotter,deformation_field=:u,colormap=:thermal,deformation_scale=2.0)
WGLMakie.current_figure()
```
Note that we can replace the plane withs some other object or a decision function. Such a function takes
the grid and a cell index as input and returns a boolean which decides whether a cell is visible or not.

### What's next?

Further, this package provides an interactive viewer that you can call with `ferriteviewer(plotter)` and
`ferriteviewer(plotter,u_history)` for time dependent views, respectively.
If you want to live plot your solution while solving some finite element system, consider to take a look at the advanced topics page.
