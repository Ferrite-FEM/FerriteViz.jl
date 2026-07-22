# FerriteViz.jl Redesign Plan — Pipeline Architecture + Fable Brief

## Context

FerriteViz.jl (v0.2.3, ~1,674 LOC in 4 files) visualizes Ferrite FEM data with
Makie, backend-agnostically. Its core idea: turn a Ferrite mesh (arbitrary,
possibly mixed cell types) into a **"discontinuous" (L2) triangle mesh** where
each cell owns its own duplicated vertices and its own reference-space
triangulation (a quad → 4 triangles through a center vertex, a hex → 6 faces ×
4). Solution fields are evaluated per duplicated vertex from the owning
element's DOFs, so element-to-element jumps are preserved. GPU acceleration
comes from storing coordinates and triangle indices in `ShaderAbstractions.Buffer`s
wrapped around Observables, shared into a `GeometryBasics.Mesh` so in-place
buffer mutation updates the GPU without rebuilding.

The code carries significant historic debt. Three problems motivate this work:

1. **Hard to extend to custom cell types.** Adding a new Ferrite cell/refshape
   requires editing 5–8 scattered sites across two files with *inconsistent
   dispatch keys* (cell-type aliases, refshapes, and interpolations):
   `ntriangles` (`utils.jl:193`), three `decompose!` methods (`utils.jl:213/254/308`),
   `linear_face_cell` (`utils.jl:19`), `midpoint` (`utils.jl:327`), plus the LOR
   tables `for_nodes`/`for_base_geometry_type`/`for_interpolation` (`lor_tools.jl:3/117/119`).
   There is no single "here is how you tessellate a cell" interface.

2. **No composable data/filter abstraction (ParaView-style).** Data transforms
   exist but are inconsistent in shape: `transfer_solution` returns an array,
   `interpolate_gradient_field` returns a new `(dh,u)`, `crinkle_clip`/
   `uniform_refinement`/`for_discretization` return new plotters/grids, and the
   `process` scalar reductions (`postprocess`, `x₁/x₂/x₃`, `l1/l2`) are ad-hoc
   kwargs. Von Mises / deviatoric stress are not in the package at all — only in
   doc examples as user closures. Users cannot compose new scenes from existing
   pieces the way ParaView chains Source → Filter → Representation.

3. **Code duplication.** The deformation/warp `@lift` block is copy-pasted across
   `SolutionPlot` (`makieplotting.jl:51-67`), `CellPlot` (`102-117`) and
   `MeshPlot` (`218-231`); the `field==:default → first field` block is repeated
   in `SolutionPlot`/`SurfacePlot`/`Arrows`; cellset coloring is duplicated
   between the two `MeshPlot` methods; the two `Elementinfo` methods
   (`AbstractCell` vs `Interpolation`, `438-498` vs `500-560`) are near
   line-for-line identical; shared `Attributes` blocks repeat. Debug leftovers
   remain: `@show` (`makieplotting.jl:247`), `@info WF` (`253`); `Arrows.normalize`
   is `#TODO: broken`; `transfer_solution` is self-labeled "peak inefficiency"
   (`utils.jl:362`); `const wireframe = meshplot` is `# FIXME` (`672`).

**Decisions (from the user):** full pipeline rewrite (new core, not additive);
breaking changes allowed, shipped as **v0.3.0**; deliverable is a single
comprehensive Fable brief (reproduced verbatim at the end of this file).

**Intended outcome:** a three-layer architecture — (1) a single cell
*tessellation interface*, (2) a reactive *DataSet + Filter graph* (ParaView
model), (3) *thin representations* (recipes) — that makes custom cell types a
one-method extension, makes scenes composable from named filters, and removes
the duplicated recipe internals. GPU-buffer sharing, Observable reactivity, live
`update!`, and Makie-backend-agnosticism are preserved; numeric correctness is
guarded by the existing analytic-field tests.

---

## Target architecture

### Layer 1 — Tessellation interface (extensibility)

One extension point replaces the scattered geometric dispatch. A cell's
reference-space tessellation is described by a single value:

```julia
struct ReferenceTessellation{refdim,T}
    coords::Vector{Ferrite.Vec{refdim,T}}   # tessellation vertex coords in REFERENCE space
    triangles::Vector{NTuple{3,Int}}        # surface triangles indexing into coords
end

# THE extension point. Implement one method to support a new cell:
reference_tessellation(::Type{<:Ferrite.AbstractRefShape}, order::Int=1)::ReferenceTessellation
```

Everything geometric is *derived* from this, not re-dispatched:
- `ntriangles(cell)` = `length(reference_tessellation(refshape, order).triangles)`.
- Physical coords: push `coords` through the cell's **geometric interpolation**
  (`Ferrite.spatial_coordinate` / `PointValues` on the geo ip). This fixes the
  "TODO use geometric interpolation here" (`utils.jl:622`) and "deal with
  nonlinear geometries" (`utils.jl:326`) issues for free — curved elements
  tessellate correctly.
- `reference_coords` per vertex = the `coords` entries (already what
  `transfer_solution` consumes).
- `midpoint` = centroid of the tessellation mapped through the geo ip.
- 3D cells: either keep the generic face-based construction as one method, or
  encode the boundary triangulation directly in `reference_tessellation`.
- **Unify LOR with tessellation.** `for_nodes` (high-order → first-order
  sub-cells) is the *same concept* at higher order. Fold `lor_tools.jl`'s tables
  into one registry keyed by `(refshape, order)`; `for_discretization` and
  `uniform_refinement` become consumers of that registry.

Provide default implementations for RefLine/RefTriangle/RefQuadrilateral/
RefTetrahedron/RefHexahedron (+ RefPrism/RefPyramid if cheap). Net effect:
adding a custom cell type = implement `reference_tessellation` (+ its geometric
interpolation, which Ferrite already needs) and nothing else.

### Layer 2 — DataSet + reactive Filter graph (ParaView model)

Replace `MakiePlotter` as the core with a **DataSet** source and a graph of
**Filters**.

- **`FEData` (Source):** wraps `dh`, `u::Observable`, the tessellation, and
  named **point-data** and **cell-data** arrays (a dictionary of derived arrays,
  like ParaView data arrays). Owns the `ShaderAbstractions.Buffer`s and the
  shared `GeometryBasics.Mesh` (GPU path unchanged). Built once; live-updatable
  via `update!(ds, u)`.

- **`AbstractFilter`:** a reactive node with `apply(f, ds::FEData)::FEData`.
  Filters hold Observables so the chain stays live: an `update!` on the source
  `u` propagates through every filter to the GPU buffers. Composable with `|>`.
  Port existing operations into uniform filters, and add the ones users
  currently hand-roll:
  - `WarpByVector(field, scale)` — was baked into recipe kwargs.
  - `Gradient(field)` — wraps `interpolate_gradient_field` (+ its
    `MatrixizedInterpolation` machinery), lifting the single-subdofhandler
    restriction where feasible.
  - `CrinkleClip(plane)` / `Clip` — was `crinkle_clip`.
  - `Refine(n)` — was `uniform_refinement`.
  - `FirstOrderRefinement()` — was `for_discretization`.
  - Point/cell-data derivations replacing `process` closures, each producing a
    **named** array: `Component(i)` (x₁/x₂/x₃), `Magnitude` (l2), `Norm1` (l1),
    `VonMises`, `Deviator`, `Threshold`. Von Mises + deviator become
    first-class (currently only in doc examples).

- **Representations (thin recipes):** `solutionplot`, `surfaceplot`, `arrows`,
  `wireframe`/`meshplot`, `cellplot` become thin: take a `FEData` and a **named
  array** to color by. No embedded field-resolution or deformation logic — those
  are filters upstream. Example:
  ```julia
  ds   = FEData(dh, u)
  pipe = ds |> WarpByVector(:u, 2.0) |> Gradient(:u) |> VonMises()
  solutionplot(pipe; color = :vonMises)
  ```

### Layer 3 — Representation dedup

With logic pushed into filters, extract the residual shared helpers and delete
dead code:
- `_resolve_array(ds, name)` for default/named array selection.
- One `_warp_buffer!(ds, coords, scale)` GPU-buffer update helper (replaces the
  three duplicated deformation blocks).
- `cellset_data(grid)` shared by the mesh/wireframe representations.
- Merge the two `Elementinfo` methods via a shared `_draw_reference_element!`
  taking the geometric interpolation + refshape.
- A shared base `Attributes` set for colormap/colorrange/nan_color/shading.
- Remove `@show`/`@info` debug lines; fix `Arrows.normalize`; drop the
  `const wireframe` FIXME alias in favor of a real exported name.

---

## Critical files

- `src/utils.jl` — becomes/soaks into Layer 1 (tessellation) + Layer 2 (`FEData`,
  filters). `MakiePlotter` (`:24`), `decompose!` (`:213/254/308`), `ntriangles`
  (`:193`), `transfer_solution` (`:364`), `interpolate_gradient_field` (`:483`),
  `crinkle_clip` (`:137/164`), `uniform_refinement` (`:587`),
  `MatrixizedInterpolation` (`:684-767`).
- `src/lor_tools.jl` — `for_nodes`/`for_*` tables folded into the tessellation
  registry; `for_discretization` (`:138`) becomes a filter.
- `src/makieplotting.jl` — recipes rewritten as thin representations; dedup
  helpers; `ferriteviewer` (`:576`) rewired onto filters/named arrays; remove
  debug leftovers (`:247/253`), FIXME alias (`:672`).
- `src/FerriteViz.jl` — new module structure + exports (add `FEData`, filter
  constructors, `apply`); consider splitting into `tessellation.jl`,
  `dataset.jl`, `filters.jl`, `representations.jl`.
- `test/runtests.jl` — keep the analytic-field numeric checks (`:31-177`, the
  `f_ana` correctness assertions at various orders/shapes); update to the new
  API surface.
- `docs/src/*.md` (`api.md`, `atopics.md`, `devdocs.md`, `tutorial.md`) +
  `docs/src/ferrite-examples/*.jl` — update to the pipeline API; move von
  Mises/deviator out of examples into documented filters.
- `Project.toml` — bump to `0.3.0`; add `[compat]` entries for any new deps
  (aim for none beyond the current set).

---

## Reuse (do not rewrite from scratch)

- The **L2 discontinuous triangulation** semantics and per-duplicated-vertex
  evaluation (`_transfer_solution!`, `utils.jl:397`) are correct — keep the
  algorithm, re-home it behind the tessellation interface.
- The **GPU pattern** (`ShaderAbstractions.Buffer` + shared `GeometryBasics.Mesh`
  + in-place `buffer[1:end] = ...`, `utils.jl:94-97`) is the right approach —
  preserve it inside `FEData`.
- `MatrixizedInterpolation` (`utils.jl:684-767`) and
  `_tensorsjl_gradient_accessor` (`:473`) are needed by the `Gradient` filter —
  reuse as-is.
- Quad 4-triangle and hex face decompositions are deliberate (they render the
  linear solution mode correctly) — encode them as the default
  `reference_tessellation`s, don't "simplify" to 2-triangle splits.

---

## Verification

1. **Unit/numeric:** run `test/runtests.jl` headless (`Pkg.test`), keeping the
   analytic-field assertions green at each order/shape (Tri/Quad/Tet/Hex,
   orders 2–3). Backend-agnostic: tests must not require a GPU; use `CairoMakie`
   or no backend for CI.
2. **Extensibility check:** add a throwaway test that registers
   `reference_tessellation` for one additional refshape and confirms
   `ntriangles`/`FEData`/`solutionplot` work with *no other edits* — this is the
   acceptance criterion for goal #1.
3. **Pipeline check:** a test composing `FEData |> WarpByVector |> Gradient |>
   VonMises` and asserting the named array matches a hand-computed value — the
   acceptance criterion for goal #2.
4. **Dedup check:** grep the representations for the removed duplicated blocks
   and debug macros; confirm the three deformation blocks collapse to one helper.
5. **Docs build:** `docs/make.jl` builds; `atopics.md` gradient/high-order/clip
   examples and `tutorial.md` run under the new API.
6. **Backend-agnostic guard:** ensure `src/` imports only `Makie` (no concrete
   backend), matching the current contract.

---

## THE FABLE BRIEF (hand this to Fable verbatim)

> ## Task: Redesign FerriteViz.jl around a tessellation interface + a reactive ParaView-style data/filter pipeline (v0.3.0, breaking changes allowed)
>
> ### Background
> FerriteViz.jl visualizes Ferrite.jl FEM data with Makie, **backend-agnostically**
> (`src/` may import only `Makie`, never GLMakie/CairoMakie/WGLMakie). The core
> idea: a Ferrite mesh (arbitrary, possibly mixed cell types) is turned into a
> **"discontinuous" (L2) triangle mesh** — every cell owns duplicated vertices and
> its own reference-space triangulation (quad → 4 triangles via a center vertex;
> hex → 6 faces × 4). Solution fields are evaluated per duplicated vertex from the
> owning element's DOFs so inter-element jumps are preserved. GPU acceleration:
> coordinates and triangle indices live in `ShaderAbstractions.Buffer`s wrapping
> Observables, shared into a `GeometryBasics.Mesh`, so in-place buffer mutation
> updates the GPU without rebuilding. `plotter.u` is an Observable and
> `update!(plotter, u)` drives live plotting.
>
> Current code is ~1,674 LOC in `src/FerriteViz.jl`, `src/utils.jl`,
> `src/makieplotting.jl`, `src/lor_tools.jl`. It works but has three structural
> problems you will fix.
>
> ### Goals
> 1. **One extension point for custom cell types.** Today, supporting a new
>    Ferrite refshape means editing 5–8 scattered sites with inconsistent dispatch
>    keys: `ntriangles`, three `decompose!` methods, `linear_face_cell`,
>    `midpoint` (all in `utils.jl`), and `for_nodes`/`for_base_geometry_type`/
>    `for_interpolation` (`lor_tools.jl`). Replace all of it with a single
>    interface:
>    ```julia
>    struct ReferenceTessellation{refdim,T}
>        coords::Vector{Ferrite.Vec{refdim,T}}   # reference-space vertices
>        triangles::Vector{NTuple{3,Int}}        # surface triangles into coords
>    end
>    reference_tessellation(::Type{<:Ferrite.AbstractRefShape}, order::Int=1)
>    ```
>    Derive `ntriangles`, the physical-coordinate tessellation (map reference
>    coords through the cell's **geometric interpolation** so curved/nonlinear
>    elements work — fixes existing `utils.jl:326,622` TODOs), per-vertex
>    reference coords, and `midpoint` generically from it. Provide defaults for
>    RefLine/Triangle/Quadrilateral/Tetrahedron/Hexahedron (Prism/Pyramid if
>    cheap), reproducing the **existing** quad-4-triangle and hex-face
>    decompositions exactly (they render the linear solution mode correctly — do
>    not simplify them). Fold the `lor_tools.jl` `for_nodes` tables into the same
>    registry keyed by `(refshape, order)`; `for_discretization` and
>    `uniform_refinement` consume it. **Acceptance:** registering
>    `reference_tessellation` for one new refshape makes `ntriangles`, the data
>    source, and `solutionplot` work with no other edits.
>
> 2. **A composable ParaView-style pipeline: Source → Filter → Representation.**
>    - **Source `FEData`:** wraps `dh`, `u::Observable`, the tessellation, and
>      named **point-data**/**cell-data** arrays (a dict of derived arrays). Owns
>      the `ShaderAbstractions.Buffer`s + shared `GeometryBasics.Mesh` (keep the
>      exact GPU pattern). `update!(ds, u)` stays live.
>    - **`AbstractFilter`:** reactive node, `apply(f, ds)::FEData`, composable via
>      `|>`, Observable-carrying so `update!` propagates through the whole chain to
>      the GPU buffers. Port existing ops to uniform filters and add the ones users
>      currently hand-roll: `WarpByVector(field,scale)` (was a recipe kwarg),
>      `Gradient(field)` (wrap `interpolate_gradient_field` +
>      `MatrixizedInterpolation`), `CrinkleClip(plane)`, `Refine(n)`,
>      `FirstOrderRefinement()`, and named data derivations replacing the `process`
>      closures: `Component(i)`, `Magnitude`, `Norm1`, `Threshold`, and — newly
>      first-class — `VonMises`, `Deviator` (currently only in doc examples).
>    - **Representations:** rewrite `solutionplot`/`surfaceplot`/`arrows`/
>      `wireframe`(`meshplot`)/`cellplot` as **thin** recipes that take a `FEData`
>      and color by a **named array** — no embedded field-resolution or deformation
>      logic (those are upstream filters). Target usage:
>      ```julia
>      ds = FEData(dh, u)
>      solutionplot(ds |> WarpByVector(:u,2.0) |> Gradient(:u) |> VonMises(); color=:vonMises)
>      ```
>    **Acceptance:** the pipeline above renders, and a test asserts the `:vonMises`
>    array matches a hand-computed value; live `update!` still updates every
>    representation through the filter chain.
>
> 3. **Remove duplication and dead code.** Collapse the three duplicated
>    deformation `@lift` blocks (`makieplotting.jl:51-67,102-117,218-231`) into one
>    `_warp_buffer!` helper; the repeated `field==:default` block into
>    `_resolve_array`; the two duplicated cellset-coloring blocks into
>    `cellset_data(grid)`; the two near-identical `Elementinfo` methods
>    (`438-498` vs `500-560`) into one shared `_draw_reference_element!`. Delete
>    debug leftovers `@show` (`:247`) and `@info WF` (`:253`), fix the `#TODO:
>    broken` `Arrows.normalize`, and replace the `# FIXME const wireframe =
>    meshplot` alias (`:672`) with a properly exported name.
>
> ### Hard constraints
> - Julia + Ferrite 1.0 API (`getspatialdim`, `geometric_interpolation`,
>   `PointValues`, `facet_to_element_transformation`, `subdofhandlers`) and
>   Makie 0.24 / GeometryBasics 0.5 / ShaderAbstractions 0.3–0.5 (see
>   `Project.toml`); add **no** new dependencies if avoidable.
> - `src/` imports only `Makie` — stay backend-agnostic.
> - Preserve the GPU buffer sharing + in-place mutation pattern and the
>   Observable-driven reactivity; `update!` must live-update.
> - Preserve the L2 discontinuous-evaluation numerics (`_transfer_solution!`) —
>   re-home it, don't change what it computes.
> - Where reasonable, lift the current single-subdofhandler / single-field
>   restrictions in `interpolate_gradient_field` and `for_discretization` to
>   support mixed grids.
> - Breaking API changes are allowed; bump `Project.toml` to **0.3.0** and update
>   all of `test/`, `docs/src/*.md`, and `docs/src/ferrite-examples/*.jl` to the
>   new API.
>
> ### Suggested internal phasing (ship each green)
> 1. Tessellation interface + registry; re-home `MakiePlotter`/`decompose!` onto
>    it; keep tests green.
> 2. `FEData` source + `update!` + GPU buffers on top of the tessellation.
> 3. Filter graph (`apply`, `|>`) + port/create the filters listed above.
> 4. Thin representations + dedup helpers; delete dead code; rewire
>    `ferriteviewer`.
> 5. Update tests, docs, examples; bump version.
>
> ### Verification (must all pass)
> - `Pkg.test` green headless (CairoMakie or no backend; no GPU required),
>   keeping the analytic-field correctness assertions in `test/runtests.jl`
>   (`f_ana` at Tri/Quad/Tet/Hex, orders 2–3).
> - New extensibility test: one extra `reference_tessellation` method →
>   `solutionplot` works with no other edits.
> - New pipeline test: `FEData |> WarpByVector |> Gradient |> VonMises` renders and
>   the named array is numerically correct.
> - `docs/make.jl` builds; `atopics.md` + `tutorial.md` examples run under the new
>   API.
> - Grep confirms the duplicated blocks and debug macros are gone and `src/`
>   imports only `Makie`.
>
> ### Additional notes on the existing example scripts
> - **Two example scripts are stale (pre-1.0 Ferrite API) — port them, don't
>   assume they run.** `docs/src/ferrite-examples/plasticity-live.jl` and
>   `plasticity-timeslider.jl` still use the pre-1.0 Ferrite API
>   (`CellVectorValues`, `FaceVectorValues`, `QuadratureRule{3,RefTetrahedron}`,
>   `getfaceset`, `nfaces`), whereas `plasticity.jl` was already updated to
>   Ferrite 1.0 (`CellValues`, `FacetValues`, `getfacetset`, `nfacets`). Update
>   these two to the Ferrite 1.0 API as part of the docs/examples pass, and
>   collapse the ~3× duplication across the three `plasticity*.jl` scripts
>   (they are ~10 KB each and nearly identical).
> - **Reuse the existing von Mises / deviator code when building the `VonMises`
>   and `Deviator` filters.** The only von-Mises/deviator implementations in the
>   repo live in those plasticity example scripts (via Tensors.jl `dev`/`tr`,
>   e.g. `plasticity.jl:50/63/65`) — not in `src/`. Lift that logic into the new
>   first-class filters rather than re-deriving it, and delete the ad-hoc copies
>   from the examples once the filters exist.
