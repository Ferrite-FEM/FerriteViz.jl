# Developer Documentation

Note that these functions could be removed or change in behavior between minor version changes! Use and dispatch on these with care!

## Architecture

FerriteViz is structured in three layers, following the ParaView model:

1. **Tessellation** (`src/tessellation.jl`): every reference shape describes its
   surface triangulation with a single [`FerriteViz.ReferenceTessellation`](@ref).
   [`FEData`](@ref) lays this out per cell with *duplicated* vertices, so
   discontinuous (L2) fields render with their inter-element jumps intact, and
   maps reference coordinates through the cell's geometric interpolation (curved
   cells tessellate correctly).
2. **Data pipeline** (`src/dataset.jl`, `src/filters.jl`): [`FEData`](@ref)
   holds the solution as an `Observable` plus named point-/cell-data arrays;
   filters derive new datasets while sharing the source observable, so
   [`FerriteViz.update!`](@ref) propagates through the entire pipeline. The
   coordinates and triangles live in `ShaderAbstractions.Buffer`s shared into a
   `GeometryBasics.Mesh` — updates mutate GPU data in place without rebuilding.
3. **Representations** (`src/representations.jl`): thin Makie recipes that take
   an `FEData` and the *name* of the array to color by.

## Adding support for a custom cell type

Implement one method. For a 3D reference shape,
[`FerriteViz.facet_based_tessellation`](@ref) usually is all you need — e.g. if
pyramids were not already supported, this would make `FEData` and all
representations work for them:

```julia
FerriteViz.reference_tessellation(::Type{Ferrite.RefPyramid}) =
    FerriteViz.facet_based_tessellation(Ferrite.RefPyramid)
```

For 2D shapes, construct the [`FerriteViz.ReferenceTessellation`](@ref)
directly (coordinates in reference space, triangles indexing into them; shared
coordinates are fine — per-cell duplication is `FEData`'s job).

To additionally support [`FirstOrderRefinement`](@ref) for a high-order
interpolation, provide its [`FerriteViz.first_order_subcells`](@ref) table (and
[`FerriteViz.linear_celltype`](@ref) for a new reference shape).

## Data layout

Point-data arrays are `Matrix{Float64}` (nvertices × ncomponents) with tensor
components in Tensors.jl linear (column-major) order; scalars have one column.
Cell-data arrays are per-cell `Vector`s of arbitrary element type.

## Reference

```@docs
FerriteViz.ReferenceTessellation
FerriteViz.reference_tessellation
FerriteViz.facet_based_tessellation
FerriteViz.first_order_subcells
FerriteViz.linear_celltype
FerriteViz.ntriangles
FerriteViz.num_vertices
FerriteViz.transfer_solution
FerriteViz.transfer_scalar_celldata
FerriteViz.interpolate_gradient_field
FerriteViz._tensorsjl_gradient_accessor
```
