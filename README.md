# FerriteViz.jl

[![Build Status](https://github.com/ferrite-fem/FerriteViz.jl/workflows/CI/badge.svg)](https://github.com/ferrite-fem/FerriteViz.jl/actions)
[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg

[docs-dev-url]: http://ferrite-fem.github.io/FerriteViz.jl/dev/

Small package to visualize your [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) results. Currently supports only Makie,
but the goal is to extend it for different plotting packages.

The package is highly experimental and breaking changes about the internal machinery are about to come.
Likely, only a small fraction of the interface will change over time.

## Installation

```julia
pkg> add FerriteViz
```

## Usage

Simply grab your solution vector and the corresponding dof handler to create a data source,
then pipe it through filters into the plots — ParaView's Source → Filter → Representation model:

```julia
import FerriteViz, GLMakie
using FerriteViz
dh, u = solve_problem()
ds = FEData(dh, u)
FerriteViz.solutionplot(ds)
# or composed, e.g. the von Mises invariant of the displacement gradient on
# the deformed mesh (apply a constitutive law with Derive to get a stress):
FerriteViz.solutionplot(ds |> WarpByVector(:u) |> Gradient(:u) |> VonMises(), color=:vonMises)
```

For a guide check out [the tutorial section](https://ferrite-fem.github.io/FerriteViz.jl/dev/tutorial.html) - or just enjoy the gallery below!

## Features

- composable, reactive filters: `WarpByVector`, `Gradient`, `CrinkleClip`, `Refine`,
  `FirstOrderRefinement`, `ExtractComponent`, `Magnitude`, `Norm1`, `VonMises`, `Deviator`,
  `Threshold`, `Derive`, ...
- `AddQuadraturePointData` renders internal variables — data known only at the quadrature points —
  by partitioning every cell into the exact Voronoi regions of its quadrature points, without
  averaging over the cell or smoothing onto a nodal field
- `solutionplot` FE solution contour plot on arbitrary finite element mesh (in Makie called `mesh` plots)
- `ferriteviewer` viewer with toggles and menus that update the plot, composable through a
  `layout(ds, state)` hook and pluggable `Control`s
- `meshplot` plots the finite element mesh and optionally labels nodes and cells
- `arrowplot` - also called `quiver` plots, in paraview `glyph` filter
- `surfaceplot` 2D solutions in 3D space as surface, in paraview `warp by scalar` filter
- custom cell types via a single `reference_tessellation` method
- `Wedge` and `Pyramid` cells out of the box; curved (higher-order geometry) cells render
  curved — surfaces *and* the `meshplot` wireframe are subdivided in reference space and
  mapped through the geometric interpolation, automatically whenever the geometry or a
  field is nonlinear (so a quadratic displacement warp bends edges too) and tunable via
  `FEData`'s `resolution`/`edge_resolution` keywords
- synchronous plotting while your simulation runs with any of the above listed options
- mutating versions of the above listed functions (except for the viewer)
- full integration into the Makie ecosystem, e.g. themes, layouts etc. 
- GPU powered plotting with GLMakie.jl, jupyter/pluto notebook plotting with WGLMakie.jl and vector graphics with CairoMakie.jl

## Missing Features

- visualization of boundary conditions
- subdomain entity plotting, e.g. facesets, edgesets and so on
- ...

For a detailed list of planned features take a look into the [issue tracker](https://github.com/Ferrite-FEM/FerriteViz.jl/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement).
Helping hands are always welcome.
Just join the discussion in the corresponding issues.

## Gallery

Pulling the Ferrite.jl logo with a [cohesive zone material model](https://github.com/kimauth/FerriteCohesiveZones.jl).

![](https://media2.giphy.com/media/v1.Y2lkPTc5MGI3NjExbHg2NXdkdnhqc2xxZ3BpZm96d3dqdjhqeHN0ODRvNmlic2hzM2RwNCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/Z8ziYf5Gju7s6FzEOA/giphy.gif)
![](https://media.giphy.com/media/v1.Y2lkPTc5MGI3NjExeDdrazJ4YmVxeDAzaGNpM2sxdGx5MGFwZm8zenUzbmZweDEzZjZpdyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/tECNrtIbG8QX0LHY7b/giphy.gif)

Credits to [Kim Auth](https://github.com/kimauth/)
