# Tutorial

## Solve a Boundary Value Problem

Start with solving a boundary value problem as you would usually do with Ferrite. It is crucial that you safe your used DofHandler
and solution vector because we need to pass those objects to `MakiePlotter`.

```@example
include("ferrite-examples/incompressible-elasticity.jl")

import FerriteVis
import GLMakie
plotter = FerriteVis.MakiePlotter(dh,u)

warp_by_vector(plotter)
```
