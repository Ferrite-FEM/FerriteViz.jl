# Developer Documentation

Note that these functions could be removed or change in behavior between minor version changes! Use and dispatch on these with care!

## Implicit adaptive subdivision (isubd)

`src/isubd.jl` holds the CPU core of view-adaptive tessellation via implicit
longest-edge bisection, in the spirit of
[demo-isubd-terrain](https://github.com/jdupuy/opengl-framework/tree/master/demo-isubd-terrain):
a persistent buffer of `UInt64` subdivision keys (base triangle id · bisection
path), a split/merge/keep streaming pass ([`FerriteViz.update_keys!`](@ref))
driven by a level-of-detail criterion, and a triangle-soup emission pass
([`FerriteViz.decode_keys!`](@ref)). It is Makie-free and unexported; the
adaptive recipes wire it into the plots' compute graphs. Every hot function is
an element-wise pass over flat buffers, so a GPU port (KernelAbstractions
kernel / Mantle compute pass) is a lowering, not a rewrite.

Pointwise field evaluation — the hot path, since the estimators query it far
more often than the rendering does — goes through per-cell monomial
coefficients ([`FerriteViz.PolyBasis`](@ref),
[`FerriteViz.PolyField`](@ref), `src/polyeval.jl`): on a fixed cell the field
is a polynomial in the reference coordinate, so it is rewritten once per
solution update and then evaluated with a handful of multiply-adds. This is
also the representation the planned fragment-shader evaluation needs. The
basis verifies itself against the shape functions before use, and
interpolations it cannot represent fall back to summing shape functions.

Refinement is driven by interpolation-error estimators
([`FerriteViz.DeviationLoD`](@ref), combined with
[`FerriteViz.CombinedLoD`](@ref)) rather than by the camera, and it is
*conforming*: [`FerriteViz.force_split!`](@ref) splits a triangle together
with the leaf across its split edge — its
[`FerriteViz.diamond_partner`](@ref), located by pure key algebra
([`FerriteViz.key_neighbour`](@ref)) over a base-adjacency table — so the
drawn surface stays watertight. The base table is built from exact global
vertex ids, and from a fan over the cells' element edges (2D) or the surface
facets' element edges (3D), which is what makes every split edge an element
edge shared by exactly two base triangles. The 3D surface facets are those
whose neighbouring cell is missing or not part of the body (`FEData.solid`),
so the adaptive path draws a closed manifold rather than every facet of every
visible cell.

That split-edge identity is also what gives [`meshplot`](@ref) an adaptive
wireframe: a point lies on an element edge exactly when its barycentric weight
for the base triangle's apex vanishes — which the key's transform gives
exactly, the bisection weights being dyadic — so the drawn segments fall out
of the same key set and coincide with the surface's own edges rather than
approximating them.

## Architecture

FerriteViz is structured in three layers, following the ParaView model:

1. **Tessellation** (`src/tessellation.jl`): every reference shape describes its
   surface triangulation *and* its wireframe edge segments with a single
   [`FerriteViz.ReferenceTessellation`](@ref) — the edges are separate from the
   triangles because the triangulation contains interior diagonals that are not
   finite element edges. [`FEData`](@ref) lays this out per cell with
   *duplicated* vertices, so discontinuous (L2) fields render with their
   inter-element jumps intact, and maps reference coordinates through the cell's
   geometric interpolation. High-order cells (or high-order fields) get their
   reference tessellation subdivided first ([`FerriteViz.subdivide`](@ref),
   driven by the [`Refine`](@ref) filter whose automatic mode `FEData`
   applies unless constructed with `adaptive=false`), so curved geometry
   and deformation render curved. Since the wireframe's vertices are ordinary
   tessellation vertices, [`meshplot`](@ref) inherits warping, clipping and
   refinement from the pipeline without any special-casing.
   `src/qptessellation.jl` adds a second, quadrature
   rule dependent reference geometry: the Voronoi partition of a reference shape
   induced by its quadrature points, which [`AddQuadraturePointData`](@ref) uses to
   render internal variables piecewise constant.
2. **Data pipeline** (`src/dataset.jl`, `src/filters.jl`): [`FEData`](@ref)
   holds the solution as an `Observable` plus named point-/cell-data arrays;
   filters derive new datasets while sharing the source observable, so
   [`FerriteViz.update!`](@ref) propagates through the entire pipeline. The
   coordinates and triangles live in `ShaderAbstractions.Buffer`s shared into a
   `GeometryBasics.Mesh` — updates mutate GPU data in place without rebuilding.
3. **Representations** (`src/representations.jl`): thin Makie recipes that take
   an `FEData` and the *name* of the array to color by. The recipes are
   new-style (declared attribute blocks) and compute derived values in the
   plot's `ComputeGraph`: `FEData`'s Observables enter the graph via
   `ComputePipeline.add_input!`, transformations are `map!` edges, and the
   child plots draw from graph nodes. Following a *named* data array when the
   name attribute changes stays Observable-side (see `resolve_color`) — which
   array a plot listens to is a structural change, and a graph edge's
   dependencies are fixed at registration.

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
directly (coordinates in reference space, triangles and wireframe edge
segments indexing into them; shared coordinates are fine — per-cell
duplication is `FEData`'s job). Edges may be omitted, in which case
[`meshplot`](@ref) draws no wireframe for cells of that shape.

## Data layout

Point-data arrays are `Matrix{Float64}` (nvertices × ncomponents) with tensor
components in Tensors.jl linear (column-major) order; scalars have one column.
Cell-data arrays are per-cell `Vector`s of arbitrary element type.

## Reference

```@docs
FerriteViz.IsubdBase
FerriteViz.IsubdMesh
FerriteViz.update_keys!
FerriteViz.refine_keys!
FerriteViz.decode_keys!
FerriteViz.decode_topology!
FerriteViz.decode_positions!
FerriteViz.UniformLoD
FerriteViz.DeviationLoD
FerriteViz.deviation
FerriteViz.CombinedLoD
FerriteViz.CachedLoD
FerriteViz.excess_levels
FerriteViz.PolyBasis
FerriteViz.PolyField
FerriteViz.key_neighbour
FerriteViz.diamond_partner
FerriteViz.force_split!
FerriteViz.conforming_update!
FerriteViz.leb_order
FerriteViz.ReferenceTessellation
FerriteViz.reference_tessellation
FerriteViz.facet_based_tessellation
FerriteViz.subdivide
FerriteViz.QPTessellation
FerriteViz.qp_voronoi_tessellation
FerriteViz.ntriangles
FerriteViz.num_vertices
FerriteViz.transfer_solution
FerriteViz.transfer_scalar_celldata
FerriteViz.interpolate_gradient_field
FerriteViz._tensorsjl_gradient_accessor
```
