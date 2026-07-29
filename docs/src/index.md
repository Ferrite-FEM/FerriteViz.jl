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
the `FEData` constructor

```julia
ds = FEData(dh,u)
```

Besides the dof fields, an `FEData` carries named point- and cell-data arrays
([`FerriteViz.set_point_data!`](@ref), [`FerriteViz.set_cell_data!`](@ref)); filters read and
write them by name, and every plot can be colored by any of them.

Now, you can use `solutionplot`, `meshplot`, `arrowplot`, `surfaceplot` or the viewer via `ferriteviewer` —
and compose transformations ParaView-style by piping the data through filters:

```julia
solutionplot(ds |> WarpByVector(:u, 2.0) |> Gradient(:u) |> VonMises(); color=:vonMises)
```

Note that the mutating `solutionplot!`, `meshplot!`, `arrowplot!` and `surfaceplot!` are available as well.

## Unique features

This package offers a set of unique features that are not easily reproducible with other export options of Ferrite.jl:

- a composable, reactive filter pipeline (ParaView's Source → Filter → Representation model):
  [`FerriteViz.WarpByVector`](@ref), [`FerriteViz.Gradient`](@ref), [`FerriteViz.CrinkleClip`](@ref),
  [`FerriteViz.Refine`](@ref), [`FerriteViz.FirstOrderRefinement`](@ref), [`FerriteViz.ExtractComponent`](@ref),
  [`FerriteViz.Magnitude`](@ref), [`FerriteViz.Norm1`](@ref), [`FerriteViz.VonMises`](@ref),
  [`FerriteViz.Deviator`](@ref), [`FerriteViz.Threshold`](@ref), [`FerriteViz.Derive`](@ref), ...
- [`FerriteViz.AddQuadraturePointData`](@ref) renders internal variables — data known only at the
  quadrature points, with no interpolation defining it elsewhere — by partitioning every cell into
  the exact Voronoi regions of its quadrature points, so the values are neither averaged over the
  cell nor smoothed onto a nodal field
- [`FerriteViz.solutionplot`](@ref) FE solution contour plot on arbitrary finite element mesh (in Makie called `mesh` plots)
- [`FerriteViz.ferriteviewer`](@ref) viewer with toggles and menus that update the plot, composable
  through a `layout(ds, state)` hook and pluggable [`FerriteViz.Control`](@ref)s
- [`FerriteViz.meshplot`](@ref) plots the finite element mesh and optionally labels nodes and cells
- [`FerriteViz.arrowplot`](@ref) - also called `quiver` plots, in paraview `glyph` filter
- [`FerriteViz.surfaceplot`](@ref) 2D solutions in 3D space as surface, in paraview `warp by scalar` filter
- synchronous plotting while your simulation runs with any of the above listed options
- mutating versions of the above listed functions (except for the viewer)
- deformed plots for any representation via the [`FerriteViz.WarpByVector`](@ref) filter
- support for custom cell types by implementing a single [`FerriteViz.reference_tessellation`](@ref)
  method (walked through on the [custom cells](cohesive.md) page)
- `Wedge` and `Pyramid` cells out of the box; curved (higher-order geometry) cells render
  curved — surfaces *and* the [`FerriteViz.meshplot`](@ref) wireframe are subdivided in
  reference space and mapped through the geometric interpolation, automatically whenever
  the geometry or a field is nonlinear (so a quadratic displacement warp bends edges too)
  and tunable via [`FEData`](@ref)'s `resolution`/`edge_resolution` keywords
- full integration into the Makie ecosystem, e.g. themes, layouts etc. 
- GPU powered plotting with GLMakie.jl, jupyter/pluto notebook plotting with WGLMakie.jl and vector graphics with CairoMakie.jl
- visualization of high order solutions via first order refinement
- visualization of non-conforming solutions, e.g. for Crouzeix-Raviart ansatz

## Viewing the docs locally

To view the docs locally use the provided live server:

```julia
include("docs/liveserver.jl")
```

Opening the html files in the browser directly might fail with a CORS error, manifesting itself figures which don't render correctly.
