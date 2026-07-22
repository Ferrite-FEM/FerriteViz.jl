# FerriteViz v0.3.0 implementation plan (Fable) — rev 3

Derived from `rewrite-plan.md` ("the brief"). Deviations from the brief are
marked **[DEVIATION]**. Changes since rev 1 respond to the codex review; the
resolution of each round-1 finding is listed at the end.

## File layout

```
src/FerriteViz.jl        module, includes, exports
src/tessellation.jl      Layer 1: tessellation interface + LOR subcell registry
src/dataset.jl           Layer 2a: FEData source, update!, data transfer, named data
src/gradient.jl          interpolate_gradient_field, MatrixizedInterpolation
src/filters.jl           Layer 2b: AbstractFilter + filters, ClipPlane
src/representations.jl   Layer 3: recipes, elementinfo, ferriteviewer
```
`utils.jl`, `makieplotting.jl`, `lor_tools.jl` are deleted (content re-homed).

## Layer 1 — tessellation.jl

```julia
struct ReferenceTessellation{refdim,T}
    coords::Vector{Ferrite.Vec{refdim,T}}   # tessellation vertices, reference space (T = Float64)
    triangles::Vector{NTuple{3,Int}}        # surface triangles indexing into coords
end
reference_tessellation(::Type{<:Ferrite.AbstractRefShape})  # THE extension point
```

**[DEVIATION]** No `order::Int=1` parameter (the brief's signature has one).
Nothing consumes an order>1 *surface* tessellation: render resolution is
`Refine(n)`'s job (reference-space subdivision re-mapped through the geometric
interpolation), and the LOR tables are interpolation-keyed (below). A vestigial
parameter would suggest an unimplemented capability; the extension contract
stays "one zero-argument method per refshape". If a per-order surface
tessellation is ever wanted, adding the argument back with a default is
non-breaking.

**Invariant (explicit):** `ReferenceTessellation` is a genuinely *indexed* mesh —
coords may be shared between triangles; triangles index into `coords`. There is
no three-coords-per-triangle expansion requirement. `FEData` construction lays
out, per cell, `length(tess.coords)` duplicated vertices (this is where the L2
inter-cell discontinuity lives — duplication is *per cell*, not per triangle)
and records both `cell_triangle_offsets` and `cell_vertex_offsets`, so
preallocation and per-cell iteration work for any user tessellation.
`transfer_*` iterates each cell's *vertex range* directly (not via triangles),
which is simpler than today and writes each vertex exactly once.

Default methods, reproducing the existing *triangle structure* exactly (same
triangles, same reference positions — buffer layout may share vertices within a
cell, which renders identically since the field is continuous inside a cell):
- `RefLine`: 2 coords, **0 triangles**. `FEData`/`meshplot` work on line grids;
  `solutionplot` renders nothing (lines have no surface) — documented.
  **[DEVIATION — interpretation]** the brief lists RefLine as a default; an
  empty triangle list is the only meaning consistent with a surface-triangle
  interface.
- `RefTriangle`: 3 coords, 1 triangle.
- `RefQuadrilateral`: 5 coords (4 corners + center), the deliberate 4-triangle
  fan (1,2,5),(2,3,5),(3,4,5),(4,1,5) — order and orientation as today.
- 3D via a generic exported helper `facet_based_tessellation(RefShape)`: for
  each face in `Ferrite.reference_faces(RefShape)`, pick the 2D tessellation by
  vertex count (3 → RefTriangle, 4 → RefQuadrilateral) and map its coords
  through `Ferrite.facet_to_element_transformation(ξ, RefShape, face_idx)`.
  Shipped: `RefTetrahedron` (4 tris), `RefHexahedron` (6×4 = 24 tris,
  reproducing today's per-face fans), `RefPrism` (Wedge), **`RefPyramid`**
  (both come free from the helper; codex round 1 rightly rejected withholding
  pyramid as artificial).

LOR registry (same file, second — optional — extension point):
```julia
first_order_subcells(ip::Ferrite.Interpolation) -> NTuple{N,NTuple{M,Int}}
linear_celltype(::Type{RefShape})               -> e.g. Triangle
```
**[DEVIATION]** The brief wants the `for_nodes` tables rekeyed to
`(refshape, order)`; codex round 1 correctly noted that node count/ordering
belongs to the *interpolation* (serendipity vs Lagrange can share
shape+order), so the tables keep interpolation keys
(`Union{Lagrange{shape,order},DiscontinuousLagrange{shape,order}}` — as today,
minus the redundant concrete-cell-type aliases), consolidated next to the
tessellation registry. `for_base_geometry_type` → `linear_celltype(refshape)`;
`for_interpolation(ip)` → `Lagrange{getrefshape(ip),1}()` inline. Rendering
extensibility needs only `reference_tessellation`; first-order refinement of a
high-order field additionally needs `first_order_subcells`.

Derived generically (all scattered per-refshape dispatch deleted):
- `ntriangles(cell) = length(reference_tessellation(getrefshape(cell)).triangles)`,
  `nvertices_tess(cell) = length(...coords)`.
- Physical coords: every tessellation vertex mapped through the cell's full
  **geometric interpolation**: `x(ξ) = Σᵢ reference_shape_value(ip_geo, ξ, i) xᵢ`
  over *all* geometry nodes (not just corners) — exact for curved cells; for
  linear cells bit-compatible with today's construction (quad center = corner
  mean, hex face center = face-corner mean). Fixes `utils.jl:326,622` TODOs.
- `midpoint(cell, points)` (cell-label placement) = the reference centroid
  `ξc = mean(reference_coordinates(ip_geo))` mapped through the geometric
  interpolation with the (possibly deformed) `points[cell.nodes]` — honors
  curved geometry, generic over refshapes.
- `linear_face_cell` deleted (`facet_based_tessellation` + `reference_faces`
  subsume it, including in `Elementinfo`).

## Layer 2a — dataset.jl: `FEData`

`MakiePlotter` becomes `FEData(dh, u; topology=...)`:

```julia
struct FEData{dim,DH,T1,TOP,M,TRI} <: AbstractPlotter
    dh::DH
    u::Makie.Observable{Vector{T1}}        # this dataset's dof vector (maybe lifted)
    source_u::Makie.Observable             # ROOT solution observable — update! target
    topology::TOP
    visible::Vector{Bool}                  # IMMUTABLE after construction (see below)
    gridnodes::Makie.Observable{Vector{GeometryBasics.Point{dim,Float32}}}  # warpable
    coords::Makie.Observable{Vector{GeometryBasics.Point{dim,Float32}}}     # tess coords (base or lifted/warped)
    coords_buffer::ShaderAbstractions.Buffer{...}                           # = Buffer(coords) → GPU
    all_triangles::Vector{TRI}
    vis_triangles::ShaderAbstractions.Buffer{TRI,...}
    triangle_cell_map::Vector{Int}
    cell_triangle_offsets::Vector{Int}
    cell_vertex_offsets::Vector{Int}       # NEW: per-cell vertex ranges
    reference_coords::Matrix{Float64}      # nverts × refdim (Float64, independent of T1)
    mesh::M                                # GeometryBasics.Mesh(coords_buffer, vis_triangles)
    point_data::Dict{Symbol,Makie.Observable}   # named per-vertex arrays (lazy cache + user data)
    cell_data::Dict{Symbol,Makie.Observable}    # named PER-CELL arrays (length ncells)
end
```

- GPU pattern unchanged: `ShaderAbstractions.Buffer` + shared
  `GeometryBasics.Mesh`; `Buffer(obs)` registers `on(obs) → update!(buffer)`
  (verified in ShaderAbstractions/src/types.jl:164), so buffers over *lifted*
  observables stay live with no manual plumbing. Storing the wrapped `coords`
  Observable (not just the buffer) gives representations like `surfaceplot` a
  liftable handle.
- **Visibility is immutable per dataset.** The mutating `crinkle_clip!` API is
  gone; `CrinkleClip` returns a new `FEData` (defensively `copy`ing the input's
  mask before modifying — masks are never mutated after construction).
  Consequently a cached point-data array can never go stale within its
  dataset — this keeps the visible-cell-skip optimization in
  `_transfer_solution!` sound (interior cells of large 3D grids are not
  evaluated until a clip exposes them, matching today's performance) while
  fixing the round-1 invalidation finding: every visibility/geometry-changing
  filter output starts with a **fresh point-data cache** recomputed against its
  own `visible`.
- **Named data:**
  - `point_data(ds, name)::Observable{Matrix{Float64}}` — cached; if `name` is
    a field of `ds.dh`, created as `lift(u -> transfer_solution(ds, u; field_name=name), ds.u)`.
    `:default` → first field of `ds.dh`. Arrays are `Matrix` (nvertices ×
    ncomponents), scalar = 1 column; tensor components use Tensors.jl linear
    (column-major) order — the documented layout contract (this is what
    `transfer_solution(...; process=identity)` produces today and what the
    tests reconstruct with `Tensor{2,dim}(row)`).
  - `cell_data(ds, name)::Observable` — **per-cell** (length ncells) arrays;
    expansion to vertices happens at use-time via the re-homed
    `transfer_scalar_celldata` (now driven by `cell_vertex_offsets`).
  - User registration: `set_point_data!(ds, name, array_or_observable)` and
    `set_cell_data!(ds, name, array_or_observable)`. Inputs are **normalized to
    the canonical layout at the boundary**: a plain `Vector` becomes an n×1
    `Matrix` (point data) / stays a per-cell `Vector` (cell data), a plain
    `Matrix` is size-validated; an `Observable` is adopted through a
    canonicalizing `lift` (so callers can keep driving their own observable
    while every stored array satisfies the row-wise contract the derivation
    filters rely on). Name collisions: overwrite, documented.
  - **Branch ownership:** every `FEData` owns its `Dict` objects. Filters never
    share a `Dict`; where entries remain valid they are carried over by
    *shallow-copying the dict* (sharing the entry Observables), so an insertion
    or overwrite on one pipeline branch can never leak into a sibling branch.
    Entry validity: geometry-compatible filters (`WarpByVector`) carry over
    both dicts; cell-count-preserving geometry filters (`CrinkleClip`,
    `Refine`) carry over `cell_data` only (vertex count/visibility changed ⇒
    `point_data` starts fresh); dh/grid-changing filters (`Gradient`,
    `FirstOrderRefinement`) start both fresh. Documented rule of thumb:
    register user arrays at the end of the pipeline.
- `FerriteViz.update!(ds, u_new)` (unexported; see Public surface) →
  writes/notifies `ds.source_u`; propagates through
  lifted `u`s, data arrays, warp buffers to the GPU. Calling it on any dataset
  in a pipeline updates the root (single source of truth).
- Re-homed internals: `transfer_solution`, `_transfer_solution!` (per-vertex
  numerics untouched; iteration switches from triangle-indices to the cell's
  vertex range — each vertex written once), `transfer_scalar_celldata`,
  `getsubdofhandlers`. `dof_to_node` is **replaced** by Ferrite's maintained
  `evaluate_at_grid_nodes` for the nodal (wireframe) path.

## gradient.jl

`interpolate_gradient_field`, `get_gradient_interpolation`,
`MatrixizedInterpolation`, `_tensorsjl_gradient_accessor` moved essentially
verbatim. **[DEVIATION — scope, codex-endorsed]** The single-subdofhandler
restriction stays as a clear, tested error and a documented limitation.

## Layer 2b — filters.jl

```julia
abstract type AbstractFilter end
apply(f::AbstractFilter, ds::FEData)::FEData
(f::AbstractFilter)(ds::FEData) = apply(f, ds)   # ds |> WarpByVector(:u, 2.0) |> ...
```

Every `apply` returns a new `FEData` (fresh buffers/mesh where geometry
changes), sharing `source_u`. Pipelines fork without clobbering each other —
unlike today, where two recipes with different `deformation_field`s fight over
one shared buffer.

Geometry filters:
- `WarpByVector(field=:default, scale=1.0)` (positional, per the brief's
  example; `scale` may be a number or an Observable). **Warp baseline = the
  input dataset's `coords` observable**, lifted directly:
  `coords = lift((base, pd, s) -> base .+ s .* to_points(pd), ds.coords,
  point_data(ds, field), scale)` — so warps compose (warp-of-warp, or warp
  after `Refine`) with no extra base-coordinate field; the displacement comes
  from the **named point-data array** (validated: ncomponents == spatial dim),
  and the result is converted to `Point{dim,Float32}` so the GPU buffer stays
  Float32. Also warps `gridnodes` the same way via `evaluate_at_grid_nodes` so
  `meshplot` draws the deformed wireframe. Carries over the input's data dicts
  (shallow copies, per the branch-ownership rule).
- `CrinkleClip(plane_or_function)` — re-homed `crinkle_clip` (decision/neighbor
  logic untouched, starting from the input's `visible`); `ClipPlane` kept.
  Chained clips now compose correctly for *data* (fresh cache per stage,
  each consistent with its own `visible`); the geometric boundary logic is the
  existing one.
- `Refine(n=1)` — re-homed `uniform_refinement`, with the fix that refined
  vertices' physical coords are computed by mapping the refined reference
  coords through the geometric interpolation (not by averaging physical
  coords) — `Refine` now genuinely improves curved geometry too.
- `FirstOrderRefinement()` — re-homed `for_discretization` (single
  subdofhandler + single field, asserted with a clear message; the misleading
  `field` argument from rev 1 is dropped). Builds the LOR grid/dh once; `u_new`
  lifted from `ds.u`.
- `Gradient(field=:default; copy_fields=Symbol[])` — builds `dh_gradient`
  once; `u_gradient` lifted from `ds.u`. Output **shares** the input's
  coords/buffers/mesh (same grid ⇒ identical tessellation), so an upstream warp
  survives: `ds |> WarpByVector(:u,2) |> Gradient(:u) |> VonMises()` works.
  Output's default field is `:gradient`.

Data-derivation filters — each inserts a named array lifted from an input
array; **association-preserving** (input resolved from point data first, then
cell data; output lands in the same dict):
- `Component(i; input=:default, output=Symbol("x",i))` — column i.
- `Magnitude(; input, output=:magnitude)` — row-wise 2-norm.
- `Norm1(; input, output=:norm1)` — row-wise 1-norm.
- `VonMises(; input, output=:vonMises)` — rows = dim² tensor components
  (Tensors linear order) → `vonmises(σ) = √(3/2 dev(σ) ⊡ dev(σ))`; the scalar
  helper `vonmises` is exported and reused by the docs examples (lifted from
  `plasticity.jl:50`).
- `Deviator(; input, output=:deviator)` — rows → `dev` components (dim² cols).
- `Threshold(; input, output=:threshold, min=-Inf, max=Inf)` — values outside
  [min,max] → NaN (renders as `nan_color`).
- `Derive(f; input=:default, output)` — generic row-wise map (ParaView
  "Calculator"): `f(row_as_reshaped_value) -> scalar/Vec/Tensor`. This is what
  the docs' stress example uses to apply the constitutive law
  (`Derive(∇u -> vonmises(σ(∇u)); output=:σvM)`), preserving the current
  example's mathematics — `VonMises` directly on a displacement gradient would
  *not* be a stress, as codex noted.

Filter-order semantics (documented): geometry-rebuilding filters
(`Refine`, `FirstOrderRefinement`, `CrinkleClip`) rebuild from base geometry —
put `WarpByVector` after them (`Gradient` is the designed exception).

## Layer 3 — representations.jl (thin recipes)

Shared helpers:
- `base_fe_attributes()` merged into each recipe's `Attributes`.
- `resolve_color(ds, colorattr)::Observable` — handles `Symbol` *and*
  `Observable{Symbol}` color attrs. Symbol resolution: `:default`/field
  name/named array → scalar data (1-column → `Vector`; informative error
  pointing at `Component`/`Magnitude` for multi-column); anything else passes
  through as a Makie color. Dynamic switching uses an explicit
  `_switching_observable(ds, name_obs)` flatten helper: an output Observable
  rewired via `on(name_obs)` — disconnect previous inner listener
  (`Observables.off`), connect the newly selected array, push its current
  value. No Observable-of-Observable, no leaked listeners.
- `cellset_data(grid)::Vector{Float64}` — deduped cellset coloring.

Recipes (all take `FEData`; `convert_arguments` keeps `recipe(dh, u)` sugar):
- `solutionplot(ds; color=:default, ...)` → `mesh!(ds.mesh, color=resolve_color(...))`.
  No `field`/`process`/`deformation_field` kwargs.
- `cellplot(ds, values::Vector; ...)` — registers the per-cell vector and
  colors by it; `cellplot(ds; color=:name)` works for named cell arrays.
- `meshplot(ds; ...)` / `meshplot(grid; ...)` — edges + nodes from
  `ds.gridnodes` (Observable → deformed when warped upstream), labels,
  cellsets. **[DEVIATION — resolved]** the `# FIXME const wireframe = meshplot`
  alias is deleted and `meshplot` *is* the properly-named replacement the brief
  asks for; a `wireframe` binding cannot be reintroduced without colliding with
  `Makie.wireframe` (imported via `using Makie`).
- `surfaceplot(ds; color=:default)` (2D only) — explicit reactive geometry:
  `positions = lift((c, s) -> Point3f.(first.(c), last.(c), s), ds.coords, scalar)`,
  drawn with `mesh!(positions, ds.vis_triangles)` so visibility and upstream
  warps (via `ds.coords`) and field updates (via the scalar array) all
  propagate; covered by the live-update test.
- `arrowplot(ds; field=:default, color=:default, ...)` — **renamed from
  `arrows`**: Makie 0.24 exports `Arrows` (alias of `Arrows2D/3D`) and
  deprecates the generic `arrows` entry point, so `@recipe(Arrows)` would
  collide with the imported binding. Internally dispatches to
  `Makie.arrows2d!`/`Makie.arrows3d!` by spatial dim; the `normalize` attribute
  is finally honored (`LinearAlgebra.normalize` applied to the vectors when
  set) — closes the `#TODO: broken`.
- `Elementinfo` — the two ~60-line methods collapse into one
  `_element_info_plot!(Ele, gip, node_ip, refshape, labelprefix)`; the
  `AbstractCell` method labels "N", the `Interpolation` method "D"; face labels
  from `reference_faces` + `facet_to_element_transformation`.
- `ferriteviewer(ds[, u_history])` — rewired: builds
  `apply(WarpByVector(field, scale_obs), ds)` once, toggle drives `scale_obs`
  (0 ↔ slider value); field/colormap/derivation menus drive the recipes' color
  attr through `_switching_observable`. Timeslider variant drives `update!`.
- Debug leftovers (`@show` at makieplotting.jl:247, `@info` at :253) deleted;
  `src/` imports only `Makie` (grep-guarded in tests).

## Public surface (breaking, v0.3.0)

Exported: `FEData`, `apply`, filters (`WarpByVector`, `Gradient`,
`CrinkleClip`, `Refine`, `FirstOrderRefinement`, `Component`, `Magnitude`,
`Norm1`, `VonMises`, `Deviator`, `Threshold`, `Derive`), `ClipPlane`,
`vonmises`, `ferriteviewer`, `set_point_data!`, `set_cell_data!`, **and the
recipe functions** `solutionplot`, `cellplot`, `meshplot`, `surfaceplot`,
`arrowplot`, `elementinfo` (+ `!` variants) — none of these names are exported
by Makie, and exporting `meshplot` satisfies the brief's "properly exported
name" for the old `wireframe` alias.
`update!` name: **unexported, called qualified as `FerriteViz.update!(ds, u)`**
— exactly today's convention (v0.2.3 doesn't export it either and the docs
already write `FerriteViz.update!`). Both Makie 0.24 and Ferrite export their
own distinct `update!`; owning an unexported `FerriteViz.update!` (defined as a
fresh function, extending neither) is the only collision-free choice, and user
scripts that also `using Ferrite` keep calling `Ferrite.update!`/`update!` for
constraints unambiguously. Docs/examples use the qualified form throughout.
Documented-but-unexported (docs use `FerriteViz.` prefix): the extension
interface `ReferenceTessellation`, `reference_tessellation`,
`facet_based_tessellation`, `first_order_subcells`, `linear_celltype`,
`AbstractFilter`, `point_data`, `cell_data`.
Gone: `MakiePlotter`, `for_discretization`, `crinkle_clip[!]`,
`uniform_refinement`, `wireframe`, `arrows`. CHANGELOG gets a migration table.

## Tests

Keep the analytic-field correctness structure (`f_ana`, Tri/Quad/Tet/Hex,
orders 2–3, scalar+vector, gradient checks); assertions go through the public
path (`point_data(ds, :u)`). New:
1. **Extensibility (mandatory, no fallback)**: a minimal test-only 2D refshape
   fixture (`RefDummy <: AbstractRefShape{2}` + cell + a tiny interpolation
   delegating its shape functions to `Lagrange{RefTriangle,1}` formulas)
   implementing the Ferrite-side contract; assert `FEData` errors informatively
   before registration, then that registering **one** `reference_tessellation`
   method makes `ntriangles`/`FEData`/`solutionplot` work with no other
   FerriteViz edits.
2. **Pipeline**: `FEData |> WarpByVector(:u,2.0) |> Gradient(:u) |> VonMises()`
   with an analytic field: `:vonMises` matches hand-computed values per visible
   vertex; warped coords = base + 2u; **a `solutionplot` figure is constructed
   from the full filtered dataset**; then `FerriteViz.update!(root, u2)` and assert the
   named array, the coords buffer, *and* the constructed plot's color observable
   all changed (live propagation).
3. **Recipe smoke tests**: solutionplot/meshplot/surfaceplot/arrowplot/
   cellplot/elementinfo figures constructed headless. Makie is declared in
   test extras/targets; if plot construction turns out to need a backend,
   CairoMakie is declared there instead (decided empirically, stated in the
   test file).
4. **Guards**: no `@show`/`@info` in `src/`; `src/` imports only Makie.

Cadence: numeric test file via the `jld` daemon per phase; full `Pkg.test` at
the end.

## Docs & examples

- `api.md`/`devdocs.md`: new surface; devdocs gains "adding your own cell
  type" (RefPyramid shown as the worked example of `facet_based_tessellation`,
  even though it now ships) + pipeline architecture notes + data-layout
  contract.
- `tutorial.md`/`atopics.md`: pipeline API throughout; stress example via
  `Gradient |> Derive(constitutive law)`; high-order via
  `FirstOrderRefinement()`/`Refine(n)`; clip via `CrinkleClip`.
- `plasticity-live.jl`/`plasticity-timeslider.jl` ported to Ferrite 1.0 and
  deduped: `plasticity.jl` keeps model + `solve(; liveplotting=false)`; the
  other two become thin includes. Local `vonMises` replaced by
  `FerriteViz.vonmises`.
- **`docs/make.jl` must build strictly** — hard acceptance. The blanket
  `warnonly=true` in `docs/make.jl:6` is removed (at most narrowed to
  documented-justified categories) so example-block errors fail the build. Run
  locally with the docs Manifest (WGLMakie/Bonito serialize scenes without a
  GPU); any failure is a bug to fix, not to report away.
- The two plasticity wrapper scripts (`plasticity-live.jl`,
  `plasticity-timeslider.jl`) are not part of the Documenter build — they get
  explicit headless smoke runs (CairoMakie) as part of verification.

## Project.toml / hygiene

- version → `0.3.0`; no new runtime deps and no removals of used ones — Makie
  of course **stays under `[deps]`** (`src/` imports it); drop `StaticArrays`
  only if verified unused. The `[extras]`/`[targets]` test list *additionally*
  declares Makie (and CairoMakie if the smoke tests need a backend) so test
  code may load them under `Pkg.test` isolation.
- CHANGELOG.md with migration guide. Branch `kc/pipeline-rewrite`, commits per
  phase.

## Phasing (each ends green on the numeric tests)

1. `tessellation.jl` + rebuilt construction inside the old shell; delete
   `decompose!`/`ntriangles`/`linear_face_cell` dispatch.
2. `FEData` (rename, `source_u`, named data, `cell_vertex_offsets`), `update!`.
3. `filters.jl` + `gradient.jl`; pipeline test.
4. `representations.jl`; smoke tests.
5. Tests/docs/examples/CHANGELOG/version; full `Pkg.test`; grep guards; docs
   build.

## Round-1 findings → resolutions

1. Tessellation invariant contradictory → **fixed**: indexed mesh, explicit
   invariant, `cell_vertex_offsets` added, transfer iterates vertex ranges.
2. Named-data model incomplete → **fixed**: `set_point_data!`/`set_cell_data!`,
   per-cell cell_data, association-preserving derivation filters, explicit
   survival/collision rules.
3. Topology filters break cached data → **fixed differently than suggested**:
   visibility made immutable per dataset + fresh caches per filter output keeps
   every cached array consistent by construction, *without* paying full
   evaluation of hidden interior cells on large 3D grids (the current skip is a
   deliberate optimization worth keeping; a clip stage recomputes only its own
   arrays).
4. Warp mis-specified → **fixed**: lifts from `point_data(ds, field)`,
   component-count validated, Float32 preserved, positional `(field, scale)`;
   nodal path via `evaluate_at_grid_nodes` (replacing `dof_to_node`).
5. Dynamic selection Observable-of-Observable → **fixed**:
   `_switching_observable` flatten helper with explicit listener rewiring/`off`.
6. `Arrows`/Makie 0.24 collision → **fixed**: `arrowplot` recipe using
   `arrows2d!`/`arrows3d!`. `wireframe` removal now argued, not silent: the
   brief asks the FIXME alias be replaced by a properly exported name;
   `meshplot` is that name, and `wireframe` cannot be defined without colliding
   with the imported `Makie.wireframe`.
7. `surfaceplot` under-specified → **fixed**: reactive positions lifted from
   `ds.coords` + scalar array, drawn against `vis_triangles`, in the live test.
8. RefLine missing / RefPyramid withheld → **fixed**: RefLine ships with an
   empty triangle list (documented meaning); RefPyramid ships; extensibility
   test uses a test-only fixture (with a stated contingency).
9. Midpoint/mapping via corners unsound → **fixed**: full geometric
   interpolation for both; midpoint maps the reference centroid.
10. Docs build weakened; pipeline test didn't construct a representation;
    test-dep declaration → **fixed**: docs build is hard acceptance; pipeline
    test constructs `solutionplot`; Makie/CairoMakie declared in test targets.
Other: `reference_coords` now Float64 (independent of dof type); 1-column →
`Vector` before Makie color; per-cell vs expanded cell data separated; tensor
layout documented as the Tensors.jl linear-order contract (a full typed
data-array abstraction is intentionally out of scope); `Derive` filter added so
docs keep computing real stresses; extension interface added to the documented
public surface; `first_order_subcells` keyed by interpolation.

## Round-2 findings → resolutions

- Named-data normalization + branch ownership → setters canonicalize at the
  boundary (adopted Observables through a canonicalizing lift); every `FEData`
  owns its dicts, carry-over is by shallow dict copy sharing entry Observables.
- `meshplot`/`update!` bindings → recipe functions are exported (`meshplot` is
  the brief's "properly exported name"); `update!` stays an unexported fresh
  `FerriteViz.update!` (extending neither Makie's nor Ferrite's), always called
  qualified.
- Warp baseline → lifted from the input's `coords` observable (no phantom
  `coords_base` field); warps compose.
- Extensibility fixture mandatory → contingency dropped. `order` removal from
  `reference_tessellation` now an explicit, argued [DEVIATION].
- Strict docs validation → `warnonly=true` removed/narrowed; plasticity wrapper
  scripts get explicit headless smoke runs.
