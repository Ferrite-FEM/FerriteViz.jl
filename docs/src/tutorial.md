# Tutorial

## Solve a Boundary Value Problem

Start with solving a boundary value problem as you would usually do with Ferrite. It is crucial that you save your used DofHandler
and solution vector because we need to pass those objects to `FEData`.

## Basics

!!! tip "Plotting Functions"
    Currently, [`FerriteViz.solutionplot`](@ref), [`FerriteViz.meshplot`](@ref), [`FerriteViz.surfaceplot`](@ref), [`FerriteViz.arrowplot`](@ref), [`FerriteViz.cellplot`](@ref) and their mutating analogues with `!` are defined for `FEData`.
    Due to the nature of the documentation we need `WGLMakie`, however, you can simply exchange any `WGLMakie` call by `GLMakie`.

### Mesh utilities

```@example 1
import WGLMakie, Bonito # hide
Bonito.Page() # hide
WGLMakie.activate!() # hide
WGLMakie.Makie.inline!(true) # hide
```

You can start by plotting your mesh

```@example 1
import FerriteViz
using FerriteViz: FEData, WarpByVector, Gradient, VonMises, CrinkleClip, ClipPlane
using Ferrite
import WGLMakie #activating the backend, switch to GLMakie or CairoMakie (for 2D) locally
WGLMakie.set_theme!(size=(800, 400)) # hide

grid = generate_grid(Hexahedron,(3,3,3))
FerriteViz.meshplot(grid,markersize=10,linewidth=2)
```

FerriteViz.jl also supports showing labels for `Ferrite.AbstractGrid` entities, such as node- and celllabels, as well as plotting cellsets.

```@example 1
grid = generate_grid(Quadrilateral,(3,3))
addcellset!(grid,"s1",Set((1,4,7)))
addcellset!(grid,"s2",Set((2,5,8)))
addcellset!(grid,"s3",Set((3,6,9)))
FerriteViz.meshplot(grid,markersize=10,linewidth=1,nodelabels=true,celllabels=true,cellsets=true)
```

### Solution field of a boundary value problem

If you solve some boundary value problem with Ferrite.jl keep in mind to save your `dh::DofHandler` and solution vector `u::Vector{T}` in some variable.
With them, we create the [`FEData`](@ref) source that all plotting functions and filters operate on.

```@example 1
include("ferrite-examples/incompressible-elasticity.jl") #defines variables dh_quadratic and u_quadratic

ds = FEData(dh_quadratic,u_quadratic)
FerriteViz.arrowplot(ds)
```

Per default, all plotting functions grab the first field in the `DofHandler` (reduced to its magnitude if vector-valued).
Color by a different field by naming it — the next plot shows the pressure instead of the displacement:

```@example 1
FerriteViz.solutionplot(ds,color=:p)
```

For certain 2D problems it makes sense to visualize the result as a `surfaceplot` plot. To showcase the combination with the mutating versions of the plotting functions,
the `solutionplot` function is plotted below the `surfaceplot` plot

```@example 1
FerriteViz.surfaceplot(ds)
FerriteViz.solutionplot!(ds,colormap=:magma)
WGLMakie.current_figure()
```

### Deformed mesh for mechanical boundary value problems

In structural mechanics we often would like to see the deformed configuration.
Deformation is a *filter*: [`WarpByVector`](@ref) displaces the geometry by a
vector field, and everything plotted downstream of it is deformed. Filters are
applied by piping a dataset into them:

```@example 1
include("ferrite-examples/plasticity.jl") #only defines solving function
u, dh, u_history, σ, κ = solve()
ds = FEData(dh,u)

warped = ds |> WarpByVector(:u)
FerriteViz.solutionplot(warped,colormap=:thermal)
FerriteViz.meshplot!(warped,markersize=10,linewidth=1)
WGLMakie.current_figure()
```

### Showing per-cell data

FerriteViz.jl also supports plotting cell data, such as the **averaged** von-Mises stress or the drag stress of the plasticity example.
```@example 1
warped2 = ds |> WarpByVector(:u, 2.0)
FerriteViz.cellplot(warped2,σ,colormap=:thermal)
FerriteViz.meshplot!(warped2,markersize=10,linewidth=1)
WGLMakie.current_figure()
```
For a more granular investigation of the stress field consult the advanced tutorial.

### Interior of a 3D domain

For 3D problems we can also inspect the interior of the domain. Currently only crinkle clipping
is implemented, as the [`CrinkleClip`](@ref) filter:
```@example 1
clipped = ds |> CrinkleClip(ClipPlane(Vec((0.0,0.5,0.5)), 0.7)) |> WarpByVector(:u, 2.0)
FerriteViz.solutionplot(clipped,colormap=:thermal)
WGLMakie.current_figure()
```
Note that we can replace the plane with some other object or a decision function. Such a function takes
the grid and a cell index as input and returns a boolean which decides whether a cell is visible or not.

### What's next?

Further, this package provides an interactive viewer that you can call with `ferriteviewer(ds)` and
`ferriteviewer(ds,u_history)` for time dependent views, respectively. The viewer is composable —
its layout and controls are built from `Makie.SpecApi` and can be fully customized (see the
[composable viewer](atopics.md#Composable-viewer) section).
If you want to live plot your solution while solving some finite element system, consider to take a look at the advanced topics page.
