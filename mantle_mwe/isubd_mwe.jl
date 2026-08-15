# Adaptive tessellation of curved FEM cells, reduced to what a GPU backend has
# to provide. No FerriteViz, no Ferrite, no Makie: a hard-coded patch of two
# curved quadratic cells, the three passes, and a software rasterizer that
# renders the reference image the GPU version has to reproduce.
#
# Run:  julia --project=. isubd_mwe.jl
#
# The three passes, in the order a frame executes them:
#
#   A  update_keys   compute      per key: split / keep / merge, from an error
#                                 estimator; compacted with classify-scan-scatter
#   B  decode_keys   mesh/draw    per surviving key: three vertices, each with
#                                 its reference coordinate ξ and its cell id
#   C  shade         fragment     per pixel: colormap(Σ c[m] ξ^e[m]) from the
#                                 cell's coefficient buffer and the interpolated ξ
#
# A and B are written as KernelAbstractions kernels and run here on the CPU
# backend; they are the ones Mantle would dispatch. C is written twice: as the
# per-pixel function a fragment shader would run, and as per-vertex colours, so
# the two can be compared — that difference is the whole reason for wanting the
# fragment path.

using KernelAbstractions
using LinearAlgebra
using Printf

const KA = KernelAbstractions

# ---------------------------------------------------------------------------
# The scene: two quadratic quadrilaterals, curved, with a scalar field on them.
# ---------------------------------------------------------------------------
# Nine nodes per cell, in the reference order below. A real pipeline reads
# these from the mesh; here they are literals so nothing has to be installed.

# reference positions of the nine Q2 nodes, in [-1,1]²
const REF_NODES = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0),
                   (0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, 0.0)]

# cell 1 and cell 2, sharing the edge from (0,0) to (0,1): x, y and the field
const CELL_X = [[-1.0, 0.0, 0.0, -1.0, -0.55, 0.0, -0.5, -1.0, -0.5],
                [0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.45, 0.0, 0.5]]
const CELL_Y = [[0.0, 0.0, 1.0, 1.0, -0.25, 0.5, 1.15, 0.5, 0.45],
                [0.0, 0.0, 1.0, 1.0, 0.3, 0.5, 0.85, 0.5, 0.5]]
const CELL_F = [[0.05, 0.9, 0.15, 0.05, 0.35, 0.5, 0.1, 0.05, 0.85],
                [0.9, 0.1, 0.05, 0.15, 0.45, 0.05, 0.1, 0.5, 0.6]]
const NCELLS = 2

# ---------------------------------------------------------------------------
# Per-cell coefficient buffers
# ---------------------------------------------------------------------------
# On a fixed cell the geometry and the field are polynomials in the reference
# coordinate, so each is stored once as monomial coefficients:
#
#     u(ξ) = Σ_m c[m] ξ₁^e₁[m] ξ₂^e₂[m]
#
# The change of basis is V c = u with V the monomials sampled at the element's
# own nodes — a 9×9 solve per cell, done once here, and the only thing the GPU
# needs uploaded. This is what a fragment shader binds.

const EXPONENTS = [(i, j) for j in 0:2 for i in 0:2]          # 9 monomials
const NMONO = length(EXPONENTS)

monomials(ξ1, ξ2) = ntuple(m -> ξ1^EXPONENTS[m][1] * ξ2^EXPONENTS[m][2], NMONO)

function coefficient_matrix()
    V = zeros(NMONO, NMONO)
    for (k, (ξ1, ξ2)) in enumerate(REF_NODES)
        V[k, :] .= monomials(ξ1, ξ2)
    end
    return inv(V)
end

"""Coefficients of every cell, as `NMONO × NCELLS` matrices."""
function build_coefficients()
    Vinv = coefficient_matrix()
    cx = zeros(NMONO, NCELLS); cy = zeros(NMONO, NCELLS); cf = zeros(NMONO, NCELLS)
    for c in 1:NCELLS
        cx[:, c] = Vinv * CELL_X[c]
        cy[:, c] = Vinv * CELL_Y[c]
        cf[:, c] = Vinv * CELL_F[c]
    end
    return cx, cy, cf
end

@inline function eval_poly(coeffs, cell, ξ1, ξ2)
    mono = monomials(ξ1, ξ2)
    acc = 0.0
    @inbounds for m in 1:NMONO
        acc += coeffs[m, cell] * mono[m]
    end
    return acc
end

@inline eval_position(cx, cy, cell, ξ1, ξ2) =
    (eval_poly(cx, cell, ξ1, ξ2), eval_poly(cy, cell, ξ1, ξ2))

# ---------------------------------------------------------------------------
# Subdivision keys
# ---------------------------------------------------------------------------
# One UInt64 per triangle: the base triangle it descends from, and the path of
# longest-edge bisections that produced it, led by a sentinel bit. Parent and
# children are shifts, so the tree is implicit and a triangle's geometry is
# reconstructed from its key alone — nothing has to be stored per triangle,
# which is exactly why this fits a GPU.
#
#   triangle (v1, v2, v3) always splits across (v1, v3) at m = (v1+v3)/2 into
#       child 0 = (v1, m, v2)      child 1 = (v2, m, v3)

@inline key_base(k::UInt64) = Int(k >> 32)
@inline key_leb(k::UInt64) = k % UInt32
@inline make_key(b::Integer, leb::UInt32) = (UInt64(b) << 32) | UInt64(leb)
@inline root_key(b::Integer) = make_key(b, UInt32(1))
@inline key_depth(k::UInt64) = 31 - leading_zeros(key_leb(k))
@inline key_parent(k::UInt64) = make_key(key_base(k), key_leb(k) >> 1)
@inline key_child(k::UInt64, i) = make_key(key_base(k), key_leb(k) << 1 | UInt32(i))
@inline is_child0(k::UInt64) = iseven(key_leb(k))

# Base triangles: each cell is fanned from its centre over its four element
# edges, so a base triangle's split edge is always an element edge — the
# property a conforming implementation needs (see the README).
const BASE_CORNERS = let corners = NTuple{3,NTuple{2,Float64}}[]
    rim = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))
    for _ in 1:NCELLS, i in 1:4
        a, b = rim[i], rim[mod1(i + 1, 4)]
        push!(corners, (b, (0.0, 0.0), a))     # split edge = (b, a), apex = centre
    end
    corners
end
const BASE_CELL = [1 + (t - 1) ÷ 4 for t in 1:length(BASE_CORNERS)]
const NBASE = length(BASE_CORNERS)

# Corners of a key's triangle in the cell's reference coordinates: fold one
# 3×3 barycentric bisection per path bit.
@inline function key_corners(k::UInt64)
    b = key_base(k)
    c1, c2, c3 = BASE_CORNERS[b]
    w1 = (1.0, 0.0, 0.0); w2 = (0.0, 1.0, 0.0); w3 = (0.0, 0.0, 1.0)
    leb = key_leb(k)
    @inbounds for i in (key_depth(k) - 1):-1:0
        m = ((w1[1] + w3[1]) / 2, (w1[2] + w3[2]) / 2, (w1[3] + w3[3]) / 2)
        if ((leb >> i) & 1) == UInt32(0)
            w1, w2, w3 = w1, m, w2          # child 0 = (v1, m, v2)
        else
            w1, w2, w3 = w2, m, w3          # child 1 = (v2, m, v3)
        end
    end
    mix(w) = (w[1] * c1[1] + w[2] * c2[1] + w[3] * c3[1],
              w[1] * c1[2] + w[2] * c2[2] + w[3] * c3[2])
    return mix(w1), mix(w2), mix(w3)
end

# ---------------------------------------------------------------------------
# The error estimator
# ---------------------------------------------------------------------------
# How badly does the flat triangle approximate the truth? Sample the deviation
# of the linear interpolation over the triangle — at the edge midpoints and the
# centroid — for the geometry and for the field, and take whichever is worse
# relative to its tolerance. Note the *edge* midpoints all matter: a
# multilinear field is exactly linear along element edges and curves only
# across the fan diagonals, so sampling only the split edge silently stops the
# refinement.

const SAMPLES = ((0.5, 0.0, 0.5), (0.5, 0.5, 0.0), (0.0, 0.5, 0.5), (1 / 3, 1 / 3, 1 / 3))

@inline function excess_levels(k::UInt64, cx, cy, cf, geo_tol, fld_tol)
    cell = BASE_CELL[key_base(k)]
    p1, p2, p3 = key_corners(k)
    x1 = eval_position(cx, cy, cell, p1...); f1 = eval_poly(cf, cell, p1...)
    x2 = eval_position(cx, cy, cell, p2...); f2 = eval_poly(cf, cell, p2...)
    x3 = eval_position(cx, cy, cell, p3...); f3 = eval_poly(cf, cell, p3...)
    geo_err = 0.0; fld_err = 0.0
    @inbounds for (a, b, c) in SAMPLES
        ξ1 = a * p1[1] + b * p2[1] + c * p3[1]
        ξ2 = a * p1[2] + b * p2[2] + c * p3[2]
        ex, ey = eval_position(cx, cy, cell, ξ1, ξ2)
        lx = a * x1[1] + b * x2[1] + c * x3[1]
        ly = a * x1[2] + b * x2[2] + c * x3[2]
        geo_err = max(geo_err, sqrt((ex - lx)^2 + (ey - ly)^2))
        fld_err = max(fld_err, abs(eval_poly(cf, cell, ξ1, ξ2) - (a * f1 + b * f2 + c * f3)))
    end
    # a level buys a factor two: the deviation is O(h²) and an edge halves
    # every second bisection
    return max(log2(max(geo_err, 1e-16) / geo_tol), log2(max(fld_err, 1e-16) / fld_tol))
end

# ---------------------------------------------------------------------------
# PASS A — update the key buffer            (Mantle: a compute dispatch)
# ---------------------------------------------------------------------------
# Per key: split (emit two children), merge (child 0 emits the parent, child 1
# emits nothing), or keep. Written as classify → prefix sum → scatter rather
# than an atomic append: the output stays in the input's order, which keeps the
# buffer sorted, and a scan is a pass every backend already has.

@kernel function classify_kernel!(counts, @Const(keys), @Const(cx), @Const(cy), @Const(cf),
                                  geo_tol, fld_tol, max_depth)
    i = @index(Global)
    k = keys[i]
    d = key_depth(k)
    if d < max_depth && excess_levels(k, cx, cy, cf, geo_tol, fld_tol) > 0
        counts[i] = 2                                   # split
    elseif d > 0 && excess_levels(key_parent(k), cx, cy, cf, geo_tol, fld_tol) <= 0
        counts[i] = is_child0(k) ? 1 : 0                # merge: one of the pair emits
    else
        counts[i] = 1                                   # keep
    end
end

# A merging child 0 emits its *parent*, a kept key emits itself, and both
# carry a count of 1 — so the scatter re-tests the predicate rather than
# smuggling the decision through the count buffer. One extra evaluation per
# key; a real implementation would pack a two-bit code instead.
@kernel function scatter_kernel!(out, @Const(keys), @Const(counts), @Const(offsets),
                                            @Const(cx), @Const(cy), @Const(cf),
                                            geo_tol, fld_tol, max_depth)
    i = @index(Global)
    k = keys[i]
    n = counts[i]
    o = offsets[i]
    if n == 2
        out[o + 1] = key_child(k, 0)
        out[o + 2] = key_child(k, 1)
    elseif n == 1
        d = key_depth(k)
        merging = d > 0 && excess_levels(key_parent(k), cx, cy, cf, geo_tol, fld_tol) <= 0 &&
                  !(d < max_depth && excess_levels(k, cx, cy, cf, geo_tol, fld_tol) > 0)
        out[o + 1] = merging ? key_parent(k) : k
    end
end

exclusive_scan(counts) = (o = similar(counts); acc = 0;
                          for i in eachindex(counts); o[i] = acc; acc += counts[i]; end; (o, acc))

"""One update pass; returns the new key buffer."""
function update_keys(backend, keys, cx, cy, cf; geo_tol, fld_tol, max_depth)
    counts = similar(keys, Int)
    classify_kernel!(backend, 64)(counts, keys, cx, cy, cf, geo_tol, fld_tol, max_depth;
                                  ndrange = length(keys))
    KA.synchronize(backend)
    offsets, total = exclusive_scan(counts)            # Mantle: a scan pass
    out = similar(keys, total)
    scatter_kernel!(backend, 64)(out, keys, counts, offsets, cx, cy, cf,
                                            geo_tol, fld_tol, max_depth; ndrange = length(keys))
    KA.synchronize(backend)
    return out
end

"""Run passes until the key set stops changing."""
function refine(backend, cx, cy, cf; geo_tol, fld_tol, max_depth = 12)
    keys = UInt64[root_key(b) for b in 1:NBASE]
    for _ in 1:(max_depth + 1)
        new = update_keys(backend, keys, cx, cy, cf; geo_tol, fld_tol, max_depth)
        new == keys && break
        keys = new
    end
    return keys
end

# ---------------------------------------------------------------------------
# PASS B — decode keys into drawable triangles   (Mantle: mesh shader or draw)
# ---------------------------------------------------------------------------
# Three vertices per key. Each vertex carries its position *and* its reference
# coordinate ξ and cell id — ξ is the varying the fragment stage needs, and the
# cell id selects the coefficient buffer.

@kernel function decode_kernel!(pos, xi, cellid, @Const(keys), @Const(cx), @Const(cy))
    i = @index(Global)
    k = keys[i]
    cell = BASE_CELL[key_base(k)]
    c = key_corners(k)
    @inbounds for v in 1:3
        j = 3 * (i - 1) + v
        xi[j] = c[v]
        pos[j] = eval_position(cx, cy, cell, c[v]...)
        cellid[j] = cell
    end
end

function decode(backend, keys, cx, cy)
    n = length(keys)
    pos = Vector{NTuple{2,Float64}}(undef, 3n)
    xi = Vector{NTuple{2,Float64}}(undef, 3n)
    cellid = Vector{Int}(undef, 3n)
    decode_kernel!(backend, 64)(pos, xi, cellid, keys, cx, cy; ndrange = n)
    KA.synchronize(backend)
    return pos, xi, cellid
end

# ---------------------------------------------------------------------------
# PASS C — shading                                (Mantle: a fragment shader)
# ---------------------------------------------------------------------------
# This is the function a fragment shader runs, once per pixel, given the
# interpolated ξ and the cell's coefficient buffer. In GLSL it is:
#
#     in vec2 xi;  flat in int cell;
#     layout(std430) buffer Coeffs { float c[]; };
#     void main() {
#         float v = 0.0;
#         for (int m = 0; m < 9; ++m)
#             v += c[cell*9 + m] * pow(xi.x, E[m].x) * pow(xi.y, E[m].y);
#         fragColor = colormap(v);
#     }
#
# — no textures, no extra geometry: the exact field, per pixel.

@inline shade(cf, cell, ξ1, ξ2) = colormap(eval_poly(cf, cell, ξ1, ξ2))

# a small perceptual-ish blue→yellow ramp, so the reference image is comparable
@inline function colormap(t)
    s = clamp(t, 0.0, 1.0)
    return (0.15 + 0.75s, 0.15 + 0.65s, 0.55 - 0.45s)
end

# ---------------------------------------------------------------------------
# Software rasterizer: the reference image the GPU version must reproduce
# ---------------------------------------------------------------------------
# `per_fragment = true` evaluates the polynomial per pixel (pass C). `false`
# evaluates it at the three vertices and interpolates the colour, which is what
# a pipeline without the coefficient buffer can do — the difference between the
# two images is the argument for the fragment path.

function render(pos, xi, cellid, cf; width = 700, height = 700,
                bounds = (-1.05, 1.05, -0.35, 1.25), per_fragment = true)
    img = fill((1.0, 1.0, 1.0), height, width)
    covered = falses(height, width)
    x0, x1, y0, y1 = bounds
    px(x) = (x - x0) / (x1 - x0) * width
    py(y) = (1 - (y - y0) / (y1 - y0)) * height
    for t in 1:(length(pos) ÷ 3)
        a, b, c = pos[3t - 2], pos[3t - 1], pos[3t]
        ξa, ξb, ξc = xi[3t - 2], xi[3t - 1], xi[3t]
        cell = cellid[3t - 2]
        ax, ay = px(a[1]), py(a[2]); bx, by = px(b[1]), py(b[2]); cx_, cy_ = px(c[1]), py(c[2])
        area = (bx - ax) * (cy_ - ay) - (cx_ - ax) * (by - ay)
        abs(area) < 1e-12 && continue
        va = per_fragment ? nothing : shade(cf, cell, ξa...)
        vb = per_fragment ? nothing : shade(cf, cell, ξb...)
        vc = per_fragment ? nothing : shade(cf, cell, ξc...)
        for j in max(1, floor(Int, min(ax, bx, cx_))):min(width, ceil(Int, max(ax, bx, cx_)))
            for i in max(1, floor(Int, min(ay, by, cy_))):min(height, ceil(Int, max(ay, by, cy_)))
                x, y = j - 0.5, i - 0.5
                w1 = ((bx - x) * (cy_ - y) - (cx_ - x) * (by - y)) / area
                w2 = ((cx_ - x) * (ay - y) - (ax - x) * (cy_ - y)) / area
                w3 = 1 - w1 - w2
                (w1 < -1e-9 || w2 < -1e-9 || w3 < -1e-9) && continue
                covered[i, j] = true
                img[i, j] = if per_fragment
                    shade(cf, cell, w1 * ξa[1] + w2 * ξb[1] + w3 * ξc[1],
                                    w1 * ξa[2] + w2 * ξb[2] + w3 * ξc[2])
                else
                    (w1 * va[1] + w2 * vb[1] + w3 * vc[1],
                     w1 * va[2] + w2 * vb[2] + w3 * vc[2],
                     w1 * va[3] + w2 * vb[3] + w3 * vc[3])
                end
            end
        end
    end
    return img, covered
end

function write_ppm(path, img)
    h, w = size(img)
    open(path, "w") do io
        write(io, "P6\n$w $h\n255\n")
        for i in 1:h, j in 1:w
            for ch in img[i, j]
                write(io, UInt8(round(clamp(ch, 0, 1) * 255)))
            end
        end
    end
    return path
end

# ---------------------------------------------------------------------------
function main()
    backend = CPU()
    cx, cy, cf = build_coefficients()
    here = @__DIR__

    @printf("base triangles: %d (2 cells × 4-triangle fan)\n", NBASE)
    for (label, geo_tol, fld_tol) in (("coarse", 2.0e-2, 5.0e-2), ("fine", 2.0e-3, 5.0e-3))
        keys = refine(backend, cx, cy, cf; geo_tol, fld_tol)
        pos, xi, cellid = decode(backend, keys, cx, cy)
        depths = map(key_depth, keys)
        @printf("%-7s geo_tol=%.0e fld_tol=%.0e → %5d triangles, depths %d–%d\n",
                label, geo_tol, fld_tol, length(keys), minimum(depths), maximum(depths))
        write_ppm(joinpath(here, "out_$(label)_fragment.ppm"),
                  first(render(pos, xi, cellid, cf; per_fragment = true)))
        write_ppm(joinpath(here, "out_$(label)_vertex.ppm"),
                  first(render(pos, xi, cellid, cf; per_fragment = false)))
    end

    # what the fragment path is worth: the same mesh, shaded both ways
    keys = refine(backend, cx, cy, cf; geo_tol = 2.0e-2, fld_tol = 1.0e9)  # geometry only
    pos, xi, cellid = decode(backend, keys, cx, cy)
    a, _ = render(pos, xi, cellid, cf; per_fragment = true)
    b, _ = render(pos, xi, cellid, cf; per_fragment = false)
    diff = maximum(maximum(abs.(x .- y)) for (x, y) in zip(a, b))
    @printf("geometry-only mesh: %d triangles; max per-pixel colour difference\n", length(keys))
    @printf("  between per-fragment and per-vertex shading: %.3f (0–1 scale)\n", diff)
    write_ppm(joinpath(here, "out_geomonly_fragment.ppm"), a)
    write_ppm(joinpath(here, "out_geomonly_vertex.ppm"), b)
    println("wrote *.ppm next to this script")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
