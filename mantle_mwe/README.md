# Adaptive tessellation of curved FEM cells — a minimal example

A self-contained reduction of what FerriteViz wants to do on the GPU, for
prototyping against [Mantle.jl](https://github.com/SimonDanisch/Mantle.jl).
No FerriteViz, no Ferrite, no Makie — two hard-coded curved quadratic cells, a
scalar field on them, three passes, and a software rasterizer producing the
image a GPU implementation has to reproduce.

```
julia --project=. isubd_mwe.jl      # writes out_*.ppm next to the script
julia --project=. convergence.jl    # what per-fragment evaluation is worth
```

## The problem, in one picture

`out_geomonly_vertex.ppm` and `out_geomonly_fragment.ppm` are the *same*
70-triangle mesh, shaded two ways. Per-vertex colours (what a pipeline without
per-cell coefficients can do) show the field as facets and streaks;
per-fragment evaluation shows it smooth, leaving only the silhouette polygonal
— on this mesh, a peak difference of 0.117 on a 0–1 colour scale.

That split is the whole design:

| | resolved by | needs |
|---|---|---|
| **geometry** — curved cells, deformation | subdividing triangles | a compute pass that rewrites the index/vertex buffer per frame |
| **solution** — the field's colour | evaluating the polynomial per pixel | a per-cell coefficient buffer bound to the fragment stage |

### What per-fragment evaluation does and does not buy

It is tempting to say "triangles are for the geometry, not for the solution".
That is half true, and `convergence.jl` measures which half.

A fragment gets its reference coordinate ξ by interpolating the triangle's
corner ξs *linearly*. On an **affine** cell that interpolation is exact, so the
polynomial is evaluated exactly where it should be and the colour is right at
any triangle count — the field's own curvature never asks for a triangle. On a
**curved** cell the drawn triangle is a chord of the true surface, so ξ is off
by the geometry's approximation error and the colour inherits it: the colour
still improves with refinement, but bounded by the *geometry* tolerance rather
than by the field's curvature.

```
CURVED cells                          AFFINE cells, same curved field
geo_tol  tris  per-frag  per-vertex   geo_tol  tris  per-frag  per-vertex
3e-02      52   0.02578     0.11581   3e-02       8   0.00000     0.12021
1e-02     120   0.01516     0.11581   1e-02       8   0.00000     0.12021
3e-03     472   0.00437     0.02904   3e-03       8   0.00000     0.12021
1e-03    1382   0.00203     0.01593   1e-03       8   0.00000     0.12021
```

(worst colour difference against a 27660-triangle reference, 0–1 scale;
refinement driven by geometry only, so nothing here refines *for* the field)

The consequence for the design is the useful part: **with per-fragment
evaluation the refinement criterion needs only its geometry term.** The
solution term — which is what FerriteViz's CPU path spends most of its time on
— exists only because per-vertex colours have to resolve the field themselves.
On the affine rows that is stark: per-vertex is stuck at 0.12 forever, because
geometry-driven refinement never adds a triangle, while per-fragment is exact
from the base mesh.

## The three passes

### A — `update_keys` (compute dispatch)

A triangle is one `UInt64`: the base triangle it descends from, and the path of
longest-edge bisections that produced it. Parent and children are bit shifts,
so the tree is implicit — no connectivity is stored, and a triangle's geometry
is reconstructed from its key alone.

Each key is classified independently (split into two children / keep / merge
into its parent) by an error estimator, then the survivors are compacted.
The MWE uses **classify → exclusive scan → scatter** rather than an atomic
append: the output keeps the input's order, which keeps the buffer sorted, and
a scan is a pass every backend already has. An atomic-append variant is fine
too if ordering does not matter to you.

*What Mantle needs:* a compute dispatch over the key buffer, a prefix-sum (or
an atomic counter), and a second dispatch writing the compacted output.
Double-buffered — the pass reads one key buffer and writes another.

### B — `decode_keys` (mesh shader, or vertex pull + indirect draw)

Each surviving key expands to three vertices. Every vertex carries

* its **position** — the cell's geometry polynomial evaluated at the corner,
* its **reference coordinate ξ** — the varying the fragment stage needs,
* its **cell id** — flat, selects the coefficient buffer.

*What Mantle needs:* either a mesh shader emitting three vertices per key, or
an indirect draw whose vertex shader pulls key `gl_VertexIndex ÷ 3` and
computes corner `gl_VertexIndex % 3`. The draw count comes from pass A's
counter, so an **indirect** draw (count read from a buffer) is required either
way.

### C — `shade` (fragment shader)

Per pixel, from the interpolated ξ and the cell's coefficients:

```glsl
in vec2 xi;  flat in int cell;
layout(std430) buffer Coeffs { float c[]; };   // 9 floats per cell here
const ivec2 E[9] = ...;                        // monomial exponents

void main() {
    float v = 0.0;
    for (int m = 0; m < 9; ++m)
        v += c[cell*9 + m] * pow(xi.x, float(E[m].x)) * pow(xi.y, float(E[m].y));
    fragColor = colormap(v);
}
```

No textures, no extra geometry: the exact field, per pixel. The same function
serves a ray-tracing hit shader unchanged — it is a pure
`(coefficients, ξ) → value`, which is why we care about it being orthogonal to
the rest of the pipeline.

*What Mantle needs:* a storage buffer readable from the fragment stage, one
interpolated `vec2` varying, and one flat integer varying.

## Where the coefficients come from

On a fixed cell an FE field *is* a polynomial in the reference coordinate, so
it is rewritten once as monomial coefficients

    u(ξ) = Σ_m c[m] ξ₁^e₁[m] ξ₂^e₂[m]

by solving `V c = u`, with `V` the monomials sampled at the element's own
nodes. That is a small dense solve per cell, done on the host whenever the
solution changes, and the result is the only per-cell data the GPU needs. For
the quadratic quadrilaterals here that is 9 coefficients for the field and 9
each for x and y.

Conditioning is worth watching at high order — a Bernstein basis is the usual
answer above roughly cubic.

## What the MWE deliberately leaves out

* **Conforming refinement.** Pass A here is the per-key rule, exactly as in
  [demo-isubd-terrain](https://github.com/jdupuy/opengl-framework/tree/master/demo-isubd-terrain).
  Refinement levels can then differ across an edge, which leaves T-vertices —
  visible as pinholes in `out_fine_fragment.ppm`. FerriteViz splits a triangle
  together with the leaf across its split edge (forcing that partner down
  first when it is coarser), which removes them structurally; on the GPU that
  neighbour query is the natural job of a concurrent binary tree. See
  `src/isubd.jl` in FerriteViz for the full rule.
* **3D.** The same construction works per surface facet; only the base
  triangles change.
* **Rate limiting.** The MWE refines to a fixed point; a real frame does one
  level per frame and lets the mesh chase the criterion.

## Numbers from the reference run

```
base triangles: 8 (2 cells × 4-triangle fan)
coarse  geo_tol=2e-02 fld_tol=5e-02 →   112 triangles, depths 3–4
fine    geo_tol=2e-03 fld_tol=5e-03 →  1106 triangles, depths 6–8
geometry-only mesh: 70 triangles; max per-pixel colour difference
  between per-fragment and per-vertex shading: 0.117 (0–1 scale)
```

Output is binary PPM (P6) so the script needs no image dependency; any viewer
or `magick out.ppm out.png` will convert it.

## Two things learned building the real thing

* **Sample every edge of a triangle in the estimator, not just the split
  edge.** A multilinear field is *exactly* linear along element edges and
  curves only across the fan diagonals; an estimator that looks only at the
  split edge stops refining trilinear hexahedra altogether.
* **Keep the estimator in double precision** even when positions render as
  `Float32`. It measures deviations far below `Float32` noise, and rounding
  them makes the criterion decide on garbage.
