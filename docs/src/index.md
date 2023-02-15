# FerriteViz.jl

FerriteViz.jl is a small package to visualize your Ferrite.jl results. Currently all Makie backends are supported and thus,
you can visualize your results in a GLMakie window, inside Pluto/Jupyter notebooks via WGLMakie and produce nice vector graphics with
CairoMakie.

In the future this package tries to adapt also other plotting packages, such as Plots.jl and PGFPlotsX.jl. Contributions are highly welcome.

## Getting Started

Install FerriteViz.jl with the in-built package manager of Julia

```julia
pkg> add FerriteViz
```

Do your computation with Ferrite.jl and save the used `DofHandler` and solution vector into a variable. Pass those two variables into
the `MakiePlotter` constructor

```julia
plotter = MakiePlotter(dh,u)
```

Now, you can use `solutionplot`, `wireframe`, `arrows`, `surface` or the viewer via `ferriteviewer`. 
Note that the mutating `solutionplot!`, `wireframe!`, `arrows!` and `surface!` are available as well.

## Unique features

This package offers a set of unique features that are not easily reproducible with other export options of Ferrite.jl:

- [`FerriteViz.solutionplot`](@ref) FE solution contour plot on arbitrary finite element mesh (in Makie called `mesh` plots)
- [`FerriteViz.ferriteviewer`](@ref) viewer with toggles and menus that update the plot
- [`FerriteViz.wireframe`](@ref) plots the finite element mesh and optionally labels nodes and cells
- [`FerriteViz.arrows`](@ref) - also called `quiver` plots, in paraview `glyph` filter
- [`FerriteViz.surface`](@ref) 2D solutions in 3D space as surface, in paraview `warp by scalar` filter
- synchronous plotting while your simulation runs with any of the above listed options
- mutating versions of the above listed functions (except for the viewer)
- deformed plots available for `solutionplot` and `wireframe`
- full integration into the Makie ecosystem, e.g. themes, layouts etc. 
- GPU powered plotting with GLMakie.jl, jupyter/pluto notebook plotting with WGLMakie.jl and vector graphics with CairoMakie.jl
- visualization of high order solutions via first order refinement
- visualization of non-conforming solutions, e.g. for Crouzeix-Raviart ansatz
