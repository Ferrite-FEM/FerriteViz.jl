# FerriteViz.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
 - The Makie recipes were migrated to the new-style `@recipe` with declared,
   documented attribute blocks, and compute derived values in the plot's
   `ComputeGraph` (Makie ≥ 0.24 / ComputePipeline) instead of Observable
   lift chains: `FEData`'s Observables enter the graph via `add_input!`,
   transformations are `map!` edges, and child plots draw from graph nodes.
   User-facing API is unchanged; recipe defaults are now visible to Makie
   (e.g. a spec-linked `Colorbar` resolves the recipe's colormap by itself),
   the plots' attribute docstrings are auto-generated, and common Makie
   attributes (`alpha`, `colorscale`, `lowclip`/`highclip`, `transparency`,
   `visible`, …) now forward to the drawn primitives. Passing an *unknown*
   keyword to a recipe is now an error instead of being silently ignored.

### Added
 - Experimental error-adaptive tessellation for `solutionplot` (#161):
   `solutionplot(ds; adaptive=true)` re-tessellates the visible cells by
   longest-edge bisection driven by two interpolation-error estimators —
   how badly the flat triangles approximate the exact geometry (dofhandler
   interpolation, including warps; `geometry_tol`, relative to the grid's
   bounding-box diagonal) and how badly the linear vertex colors approximate
   the exact field polynomial (`solution_tol`, relative to the field's value
   span) — refining where either asks. Refinement is *conforming*: a triangle
   is always split together with the leaf across its split edge (forced down
   first when it is coarser), so every drawn edge is a full edge of the
   triangle on the other side and the rendered surface is watertight — no
   hanging nodes, gaps or color seams at refinement-level boundaries. The
   estimators additionally measure every triangle edge, so the residual
   deviation of any drawn edge is bounded by the tolerance. Pass
   `conforming=false` for the cheaper independent-per-triangle refinement.
   In 3D the adaptive path extracts the actual surface — the facets whose
   neighbour is missing or was removed by a [`CrinkleClip`](@ref) — instead of
   tessellating every facet of every visible cell, so it draws a closed
   manifold with several times fewer triangles than the static path.
   [`meshplot`](@ref) takes the same `adaptive` flag: the wireframe is then
   the element edges as the refinement subdivided them, rather than a fixed
   number of segments per edge. Paired with an adaptive `solutionplot` at the
   same `geometry_tol`, every wireframe segment is an edge of the drawn
   surface exactly, so the lines cannot drift off the surface they trace.
   The whole chain (solution →
   subdivision keys → decoded triangles → per-vertex field evaluation) lives
   in the plot's ComputeGraph and follows `FerriteViz.update!`; the camera is
   not an input unless the optional screen-space criterion is enabled with
   `px_target`. Geometry, connectivity and colors are emitted by a single
   graph edge, so rapid event bursts can never render them against different
   refinement states. The color must be a dof field name or a plain color
   (registered point-data arrays live on the static tessellation and cannot
   be resampled). `FEData` records upstream `WarpByVector` applications in a
   new `deformation` field so the displaced geometry can be evaluated at
   arbitrary reference coordinates, and a `solid` field recording which cells
   make up the body (as opposed to `visible`, the cells contributing surface)
   so the surface facets can be identified after a clip.
   Benchmarked against uniform `Refine` at matched geometry *and* solution
   error (`benchmarks/adaptive_vs_uniform.jl`): a localized feature needs
   3–10× fewer triangles adaptively, while a globally smooth field (where
   uniform refinement is near optimal) saves only about 12 %. Updates cost
   3–9× less than the first working version — field values are summed straight
   from the reference shape functions instead of through `PointValues` (3×,
   and again that much for vectorized interpolations), the estimators sample
   four points instead of sixteen, decoded vertices are shared within a cell
   (5× fewer), connectivity is a separate graph node from the values it
   carries, the nodes hand out their buffers instead of copies, and every
   field is rewritten once per update into per-cell monomial coefficients so a
   sample is a few multiply-adds. What remains is dominated by pointwise
   evaluation inside the estimators: deciding a mesh costs about 30 field
   evaluations per triangle against the ~1 that filling a fixed one needs, so
   the adaptive path is not yet faster in wall clock even where it draws far
   fewer triangles.
 - `mantle_mwe/`: a self-contained example of the adaptive pipeline with no
   FerriteViz, Ferrite or Makie dependency — two hard-coded curved cells, the
   key-update and decode passes as KernelAbstractions kernels, per-cell
   coefficient buffers, per-fragment field evaluation, and a software
   rasterizer producing the reference image. Written for prototyping the
   pipeline against a GPU backend; its README states exactly which primitives
   the backend has to provide.
 - Internal (unexported) CPU core for view-adaptive tessellation via implicit
   longest-edge bisection (`src/isubd.jl`, in the spirit of jdupuy's
   demo-isubd-terrain, #161): `UInt64` subdivision keys, a split/merge/keep
   streaming pass with an isotropic screen-space LoD criterion, and a
   buffer-reusing triangle-soup decode carrying per-vertex reference
   coordinates. Not yet wired into any recipe.
 - `ReferenceTessellation` now carries the wireframe edge segments of its
   reference shape (from `Ferrite.reference_edges`) next to the surface
   triangles, and `FEData` instantiates them per cell like the triangles:
   `all_edges`/`edge_cell_map`/`cell_edge_offsets` index into the same
   tessellation vertices as the surface. Custom shapes without edges simply
   render no wireframe.
 - `FEData` grew an `adaptive::Bool=true` keyword: by default it tessellates
   with the reworked `Refine` filter's automatic mode (see below), so curved
   and high-order-deformed cells render curved out of the box;
   `adaptive=false` opts out and keeps the flat base tessellation, e.g. to
   subdivide only one branch of a pipeline with an explicit `ds |> Refine(n)`.
   The tessellation the flag picks may change in a future release; such a
   change is breaking. `adaptive=false` and explicit `Refine(n)` counts are
   stable.

### Changed
 - `Refine` is reworked (breaking): instead of relatively subdividing every
   triangle of the current tessellation into 4 dedicated-vertex triangles, it
   re-tessellates every cell from its reference shape with reference-space
   subdivision (`FerriteViz.subdivide`, shared vertices — several times less
   memory for the same picture) and now takes two counts, `surface` rounds
   for the triangles and `edges` rounds for the wireframe
   (`Refine(n; edges=n)` / `Refine(surface=..., edges=...)`). Counts are
   *absolute*, not relative to the input's tessellation, and a count of
   `nothing` (the default) is chosen per cell type: no subdivision when
   geometry and all fields are (multi-)linear, otherwise 1 surface and 3 edge
   rounds — the automatic mode `FEData` applies by default. Note `Refine()`
   now means this automatic mode, not one relative round, and applying it
   after `AddQuadraturePointData` no longer refines the quadrature-point
   partition (which was moot anyway: the partition renders piecewise-constant
   data exactly, and the refined dataset dropped the quadrature-point arrays).
 - `meshplot` draws the wireframe from the dataset's tessellation edges
   instead of connecting grid nodes with straight lines. The wireframe now
   respects the whole pipeline: it follows `WarpByVector` (including
   high-order, discontinuous and non-dof-field warps), is hidden with cells
   removed by `CrinkleClip` (previously clipped cells kept their edges), is
   refined by `Refine`, and bends along curved cell edges. Node markers and
   labels are restricted to nodes of visible cells; in 3D, edges of interior
   cells are no longer drawn.
 - `AddQuadraturePointData` datasets keep the finite element cell edges in
   their rebuilt vertex layout (valued by the nearest quadrature point), so
   `meshplot` keeps working downstream of it.
 - datasets with high-order geometry or fields tessellate finer by default
   (see above), costing those cell types about 4× the triangles; pass
   `adaptive=false` to `FEData` to restore the previous flat tessellation and
   its memory footprint.

### Removed
 - the `FirstOrderRefinement` filter and its extension-point registry
   (`FerriteViz.first_order_subcells`/`FerriteViz.linear_celltype`)
   (breaking): the reworked `Refine` covers resolving high-order fields —
   without the flattening artifacts of the first-order re-discretization —
   and `FEData`'s adaptive default already renders them curved out of the
   box.

## [0.3.0] - 2026-07-28

Full rewrite of the internals around a ParaView-style pipeline:
Source (`FEData`) → Filters → Representations. Breaking release.

### Added
 - `FEData` data source with named point-/cell-data arrays
   (`point_data`/`cell_data`/`set_point_data!`/`set_cell_data!`)
 - composable, reactive filters (applied with `|>`): `WarpByVector`,
   `Gradient`, `CrinkleClip`, `Refine`, `FirstOrderRefinement`, `ExtractComponent`,
   `Magnitude`, `Norm1`, `VonMises`, `Deviator`, `Threshold`, `Derive`
 - `AddQuadraturePointData` filter for internal variables (data known only at the
   quadrature points, with no interpolation defining it elsewhere): every cell is
   partitioned into the exact Voronoi regions of its quadrature points and each
   region is filled with that point's value, so the data is neither averaged over
   the cell nor smoothed onto a nodal field. Accepts `values[cell][qp]`,
   `values[cell, qp]` or an `Observable` of either, an `extract` function for
   material-state structs, and one quadrature rule per reference shape for mixed
   grids. The result is point data, so `VonMises`, `Deviator`, `Derive`, ...
   compose with it; two `AddQuadraturePointData` filters sharing a quadrature rule
   reproduce the same vertex layout and keep each other's arrays, which is what
   allows several quadrature point quantities to be combined in one `Derive`.
 - `Derive` accepts several inputs: `input` may be a vector of names and the
   function then takes one argument per name, taken from the same tessellation
   vertex (point data) or cell (cell data), e.g.
   `Derive((σ, εᵖ) -> σ ⊡ εᵖ; input = [:σ, :εᵖ])`. A single `Symbol` keeps
   working. All inputs must be of the same kind (all point or all cell data).
 - composable, `Makie.SpecApi`-based `ferriteviewer`: a `layout(ds, state)` hook
   returns a `GridLayoutSpec`, pluggable `Control`s (`FieldMenu`, `ProcessMenu`,
   `ColormapMenu`, `LabelsToggle`, `DeformationToggle`, `TimeSlider`) feed the
   view state, and spec helpers (`panelspec`, `solutionplotspec`, …) build the
   panels. The defaults reproduce the previous single-panel view; everything is
   overridable. `panelspec` resolves the array a panel colors by *name*, so a
   linked `Colorbar` shows the data range instead of Makie's `(0, 1)` fallback.
 - one-method extension point for custom cell types:
   `reference_tessellation(::Type{<:AbstractRefShape})`
   (with `facet_based_tessellation` for 3D shapes)
 - out-of-the-box support for `Wedge` (`RefPrism`) and `Pyramid`
   (`RefPyramid`) cells
 - curved (higher-order geometry) cells now tessellate through the geometric
   interpolation; `Refine` improves geometry resolution, not just the solution
 - `vonmises(σ)` helper (previously only in doc examples; not exported, use
   `FerriteViz.vonmises` or import it explicitly)

### Modified (breaking)
 - `MakiePlotter(dh, u)` → `FEData(dh, u)` (which copies `u`; `update!`
   no longer mutates the caller's vector)
 - recipe `field`/`process`/`deformation_field`/`deformation_scale` keyword
   arguments are replaced by upstream filters and the `color` attribute:
   `solutionplot(plotter, field=:p)` → `solutionplot(ds, color=:p)`,
   `solutionplot(plotter, deformation_field=:u)` →
   `solutionplot(ds |> WarpByVector(:u))`,
   `solutionplot(plotter, field=:gradient, process=f)` →
   `solutionplot(ds |> Gradient(:u) |> Derive(f, output=:name), color=:name)`
 - `arrows` → `arrowplot` (Makie 0.24 owns the `Arrows` recipe name);
   its `normalize` attribute works now
 - `crinkle_clip(!)` → `CrinkleClip` filter, `uniform_refinement` → `Refine`,
   `for_discretization` → `FirstOrderRefinement` (which now also supports
   vector fields)
 - `wireframe` alias removed; the recipe is `meshplot` (now exported, along
   with the other recipe functions and the filters)
 - `for_nodes`/`for_base_geometry_type`/`for_interpolation` →
   `first_order_subcells`/`linear_celltype`

### Changed
 - `ferriteviewer` internals fully rewritten around `PlotSpec`. The
   `ferriteviewer(ds)` / `ferriteviewer(ds, u_history)` signatures are preserved.
 - the documentation is restructured around the pipeline: the tutorial covers the
   plotting recipes first, then chaining filters, with `Gradient |> Derive` shown
   on the mixed displacement/pressure formulation and `AddQuadraturePointData |>
   Derive` on the plastic work density of the plasticity example. The recommended
   practices page (formerly "advanced topics") no longer repeats that material and
   focuses on the reasoning behind the discontinuous gradient and the quadrature
   point partition.

### Removed
 - `transfer_solution`'s `process` keyword, `postprocess`, `x₁`/`x₂`/`x₃`,
   `l1`/`l2` (use the data-derivation filters)
 - unused `StaticArrays` dependency

### Fixed
 - CairoMakie renders again ([#118][github-118], [#146][github-146]): its software
   mesh path expects plain `Vector`s and does not accept a `ShaderAbstractions.Buffer`
   for the faces, so the representations unwrap the shared buffers to their
   underlying arrays when CairoMakie is the active backend (GL/WGLMakie keep the
   buffer-backed, live-updating path).
 - point data holding a tensor whose dimension differs from the grid's — a shell
   or plane-strain problem carrying 3D stresses on a 2D grid — no longer dies
   with a `MethodError` from inside Tensors.jl. Rows of 4, 6 and 9 components are
   interpreted as `Tensor{2,2}`, `SymmetricTensor{2,3}` and `Tensor{2,3}`, so
   `VonMises`, `Deviator` and friends work on them; widths that cannot be
   interpreted now report the component count instead.
 - the constitutive law in the gradient example used `ones` (a tensor of ones)
   instead of `one` (the identity) for the volumetric term.

## [0.2.3] - 2026-06-22
### Added
 - `colorrange` attribute for `FerriteViz.surface` ([#122][github-122])
 - logo to the docs ([#124][github-124])

### Modified
 - Ferrite 1.0 and Makie 0.24 compatibility ([#103][github-103])
 - use `Makie.automatic` for `colorrange` instead of a manual min/max check ([#120][github-120])
 - modernize CI ([#137][github-137])
 - bump `julia-actions/setup-julia` to v2 ([#123][github-123])
 - bump `actions/checkout` to v5 ([#138][github-138]) and v6 ([#139][github-139])
 - update tutorial docs ([#126][github-126])

## [0.2.2] - 2023-11-10
### Added
 - uniform refinement for high-order solutions ([#97][github-97])
 - dependabot for GitHub actions ([#101][github-101])
 - attempt to increase internal machinery test coverage ([#104][github-104])

### Modified
 - `README.md` improvements with example gifs ([#96][github-96])
 - CI trigger only for PRs and master ([#105][github-105])
 - update docs to Documenter v1 ([#106][github-106])
 - update Makie in docs to v0.19.12 ([#109][github-109])

### Fixed
 - 0 `ntriangles` for empty domains ([#92][github-92])
 - correct link for plasticity example ([#93][github-93])
 - colorbar for 0 values and `ferriteviewer` deformation default changed to false ([#95][github-95])

## [0.2.1] - 2023-05-24
### Added
 - Basic culling where all faces of all boundary elements are rendered ([#56][github-56]).
 - Citation file ([#65](github-65))
 - Support for MixedDofHandler ([#70][github-70])
### Modified
 - Removed unnecessary extra dispatches for three-dimensional case ([#56][github-56]).
 - function barrier for `transfer_solution` such that its closer to type groundedness ([#68][github-68]).
 - `MakiePlotter` holds now `ShaderAbstractions.Buffer`s ([#69][github-69])
    - triangles are now stored in `Buffer`s with Observables
    - triangle coords are now `Buffers`s with Observables
- replace overcomplicated ternary operators by `begin end` expressions ([#69][github-69])
- remove unused functions ([#69][github-69])
- default linear rendering of high order triangles ([#83][github-83])
- keyword argument `copy_fields` added to `interpolate_gradient_field` ([#83][github-83])
### Fixed
 - Renamed `Crincle` to `Crinkle` ([#56][github-56]).
 - wireframe plot could not selectively disable the plotting of the nodes ([#83][github-83])
 - let CI error if example block errors ([#71][github-71])
 - removed bug in `transfer_solution` from ([#70][github-70]) in ([#89][github-89])
 - fix JSServe documentation issue ([#85][github-85])

## [0.2.0] - 2023-03-06
### Added
 - Functionality to obtain a first-order refined mesh and the corresponding
   dof handler and solution to approximately visualize high order solutions ([#57][github-57]).
 - Subtitles for the tutorial to find useful stuff faster ([#57][github-57]).
 - Crincle clip in 3D ([#56][github-56]), which basically removes all elements above some surface,
   which can be described by a function, from visualization.
 - `ClipPlane` to describe planar surfaces in 3D ([#56][github-56]).
 - Docs and helper for gradient field visualization based on interpolation ([#51][github-51]).
   Currently only useful in 2D, because we have no clip plane feature to introspect the interior
   of 3D problems.
 - Manufactured heat problem to test correctness of gradient field computation and as a
   helper to generate scalar-valued solutions with different ansatz ([#51][github-51]).

### Modified
 - Incompressible elasticity solver now takes the Ansatz functions and the actual material
   parameters instead of the poisson number the Ansatz functions ([#51][github-51]).

### Fixed
 - Visualization of non-conforming solution fields in 3D ([#59][github-59]).
 - An unknown bug has been fixed, which computes the colorbar `(min,max)` wrong. Now the `max` is
   set to be `1.01` of `min` guaranteeing that the value is larger than `min` if close to zero ([#51][github-51]).
 - Update Makie dependencies to fix some visualization bugs ([#51][github-51]).

[github-51]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/51
[github-56]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/56
[github-57]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/57
[github-59]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/59
[github-65]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/65
[github-63]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/63
[github-68]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/68
[github-69]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/69
[github-70]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/70
[github-71]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/71
[github-83]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/83
[github-85]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/85
[github-89]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/89
[github-92]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/92
[github-93]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/93
[github-95]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/95
[github-96]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/96
[github-97]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/97
[github-101]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/101
[github-104]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/104
[github-105]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/105
[github-106]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/106
[github-109]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/109
[github-103]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/103
[github-120]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/120
[github-122]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/122
[github-123]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/123
[github-124]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/124
[github-126]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/126
[github-137]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/137
[github-138]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/138
[github-139]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/139
[github-118]: https://github.com/Ferrite-FEM/FerriteViz.jl/issues/118
[github-146]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/146

[Unreleased]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.3...v0.3.0
[0.2.3]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.1...v0.2.0
[0.2.0]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.0...v0.1.4
