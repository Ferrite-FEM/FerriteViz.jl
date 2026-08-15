# Does the *colour* still need triangles once it is evaluated per fragment?
#
# Partly — and the split is worth being precise about, because it decides what
# the refinement criterion has to contain.
#
# A fragment gets its reference coordinate ξ by interpolating the triangle's
# corner ξs *linearly*. On an affine cell that interpolation is exact, so the
# polynomial is evaluated at exactly the right place and the colour is exact at
# any triangle count: the field's own curvature never asks for a triangle. On a
# curved cell the drawn triangle is a chord of the true surface, so ξ is off by
# the geometry's approximation error, and the colour inherits that error —
# it improves with refinement, but bounded by the *geometry* tolerance rather
# than by the field's curvature.
#
# So: per-vertex shading needs triangles for the field itself; per-fragment
# shading needs them only for the geometry, and the colour follows along.
#
# Run:  julia --project=. convergence.jl

include("isubd_mwe.jl")

# An affine variant of the same scene: two parallelograms (mid-side nodes at
# the true midpoints), carrying the identical curved field.
const AFFINE_X = [[0.0, 1.0, 1.3, 0.3, 0.5, 1.15, 0.8, 0.15, 0.65],
                  [1.0, 2.0, 2.3, 1.3, 1.5, 2.15, 1.8, 1.15, 1.65]]
const AFFINE_Y = [[0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5],
                  [0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5]]

function coefficients_from(xs, ys, fs)
    Vinv = coefficient_matrix()
    cx = zeros(NMONO, NCELLS); cy = zeros(NMONO, NCELLS); cf = zeros(NMONO, NCELLS)
    for c in 1:NCELLS
        cx[:, c] = Vinv * xs[c]; cy[:, c] = Vinv * ys[c]; cf[:, c] = Vinv * fs[c]
    end
    return cx, cy, cf
end

"""Worst colour difference over the pixels both renders covered."""
function image_error(a, ca, b, cb)
    err = 0.0
    for i in eachindex(ca)
        (ca[i] && cb[i]) || continue
        err = max(err, maximum(abs.(a[i] .- b[i])))
    end
    return err
end

function study(name, cx, cy, cf, bounds)
    backend = CPU()
    # ground truth: refine hard, evaluate per fragment
    ref_keys = refine(backend, cx, cy, cf; geo_tol = 5.0e-5, fld_tol = 1.0e9, max_depth = 14)
    rp, rx, rc = decode(backend, ref_keys, cx, cy)
    ref, refcov = render(rp, rx, rc, cf; bounds, per_fragment = true)

    println("\n", name, "  (reference: ", length(ref_keys), " triangles)")
    @printf("  %-9s %-10s %-14s %-14s\n", "geo_tol", "triangles", "per-fragment", "per-vertex")
    for tol in (3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3)
        # refine on GEOMETRY ONLY, so the field never asks for a triangle
        keys = refine(backend, cx, cy, cf; geo_tol = tol, fld_tol = 1.0e9, max_depth = 14)
        pos, xi, cid = decode(backend, keys, cx, cy)
        f, fc = render(pos, xi, cid, cf; bounds, per_fragment = true)
        v, vc = render(pos, xi, cid, cf; bounds, per_fragment = false)
        @printf("  %-9.0e %-10d %-14.5f %-14.5f\n", tol, length(keys),
                image_error(f, fc, ref, refcov), image_error(v, vc, ref, refcov))
    end
end

function main_convergence()
    cx, cy, cf = build_coefficients()
    study("CURVED cells (the MWE scene)", cx, cy, cf, (-1.05, 1.05, -0.35, 1.25))
    ax, ay, af = coefficients_from(AFFINE_X, AFFINE_Y, CELL_F)
    study("AFFINE cells, same curved field", ax, ay, af, (-0.05, 2.35, -0.05, 1.05))
    println("""
    Read the columns: per-vertex shading needs triangles to get the colour
    right in both scenes. Per-fragment shading needs none at all on affine
    cells — the field's curvature is free — and on curved cells its error is
    the geometry's, so it falls with the geometry tolerance rather than with
    the field.""")
end

main_convergence()
