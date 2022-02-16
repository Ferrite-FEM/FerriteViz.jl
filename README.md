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


## Features

- `solutionplot` FE solution contour plot on arbitrary finite element mesh (in Makie called `mesh` plots)
- `ferriteviewer` viewer with toggles and menus that update the plot
- `wireframe` plots the finite element mesh and optionally labels nodes and cells
- `arrows` - also called `quiver` plots, in paraview `glyph` filter
- `surface` 2D solutions in 3D space as surface, in paraview `warp by scalar` filter
- synchronous plotting while your simulation runs with any of the above listed options
- mutating versions of the above listed functions (except for the viewer)
- deformed plots available for `solutionplot` and `wireframe`
- full integration into the Makie ecosystem, e.g. themes, layouts etc. 
- GPU powered plotting with GLMakie.jl, jupyter/pluto notebook plotting with WGLMakie.jl and vector graphics with CairoMakie.jl

## Missing Features

- correct visualization of nonlinear geometry faces/edges
- visualization of boundary conditions
- subdomain entity plotting, e.g. cellsets, facesets, edgesets and so on
- labels for face, edge, vertex and dof numbering (technically there is no face and edge numbering in Ferrite.jl)
- plotting gradient, divergence or other fields that can be derived from the FE solution
