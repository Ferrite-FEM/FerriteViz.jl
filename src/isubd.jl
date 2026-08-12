# CPU core of implicit adaptive subdivision via longest-edge bisection (LEB),
# in the spirit of jdupuy's demo-isubd-terrain (#161). Makie-free; every hot
# function is an element-wise pass over flat buffers, so a later GPU port
# (KernelAbstractions kernel / Mantle compute pass) is a lowering, not a
# rewrite: `update_keys!` is the per-key streaming kernel (the CPU `push!`
# maps to an atomic-counter append on the GPU), `decode_keys!` the
# mesh-emission pass. Nothing here is exported yet — the adaptive recipes
# (Phase 2) wire it into the plots' compute graphs.
#
# A subdivision key is one UInt64:
#
#     [ base triangle id : 32 | LEB path : 32 ]
#
# The LEB half carries a leading 1 sentinel; the bits below it are the
# bisection path (0 = first child, 1 = second). Children and parent are one
# shift away, so the whole binary tree is implicit and a triangle's geometry
# is reconstructed from its key alone — no connectivity is stored.
#
# Corner convention: a triangle (v₁, v₂, v₃) is always bisected across the
# edge (v₁, v₃) — its "hypotenuse" — at the midpoint m = (v₁ + v₃)/2, into
#
#     child 0: (v₁, m, v₂)        child 1: (v₂, m, v₃)
#
# Each child's split edge is one of the parent's legs, which is what keeps
# recursive bisection shape-stable (newest-vertex bisection). Base triangles
# must be ordered so their longest edge is (v₁, v₃) — `leb_order` does this —
# and neighbouring base triangles should share their split edges pairwise
# ("diamonds", like the demo's two-triangle quad): the split decision is a
# symmetric function of the physical split edge, so the two triangles sharing
# it always agree, and no T-junction opens along it. Across leg-neighbours,
# crack-freeness rests on the LoD varying by less than one level between
# adjacent triangles — the same smoothness argument the demo uses.

const LEB_MAX_DEPTH = 30

key_base(k::UInt64) = (k >> 32) % Int
key_leb(k::UInt64) = k % UInt32
make_key(base::Integer, leb::UInt32) = (UInt64(base) << 32) | UInt64(leb)
root_key(base::Integer) = make_key(base, UInt32(1))
key_depth(k::UInt64) = 31 - leading_zeros(key_leb(k))
key_parent(k::UInt64) = make_key(key_base(k), key_leb(k) >> 1)
key_children(k::UInt64) =
    (make_key(key_base(k), key_leb(k) << 1), make_key(key_base(k), key_leb(k) << 1 | UInt32(1)))
is_child0(k::UInt64) = iseven(key_leb(k))

# Child corners in parent barycentrics, one matrix per bisection bit; columns
# are the child's corners. Folding them over the key's path bits gives the
# sub-triangle of the base triangle in a handful of 3×3 multiplies.
const LEB_B0 = Tensor{2,3}((1.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 1.0, 0.0))
const LEB_B1 = Tensor{2,3}((0.0, 1.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 1.0))

function key_xform(k::UInt64)
    X = one(Tensor{2,3,Float64})
    l = key_leb(k)
    for i in (key_depth(k) - 1):-1:0
        X = X ⋅ (((l >> i) & 1) == UInt32(0) ? LEB_B0 : LEB_B1)
    end
    return X
end

"""
    IsubdBase(corners, mapping)

The base domain of an implicit LEB subdivision: one reference-space corner
triple per base triangle, ordered for bisection (see [`leb_order`](@ref)),
plus the geometry `mapping(base_id, ξ) -> physical point` (for FEM cells the
geometric map, optionally composed with a warp field). The mapping is applied
per decoded vertex, so curved cells subdivide into curved sub-triangles.
"""
struct IsubdBase{RV,F}
    corners::Vector{NTuple{3,RV}}
    mapping::F
end

root_keys(base::IsubdBase) = [root_key(i) for i in 1:length(base.corners)]

# Reference-space corners of the key's sub-triangle.
function key_corners(base::IsubdBase, k::UInt64)
    r1, r2, r3 = base.corners[key_base(k)]
    X = key_xform(k)
    return (X[1, 1] * r1 + X[2, 1] * r2 + X[3, 1] * r3,
            X[1, 2] * r1 + X[2, 2] * r2 + X[3, 2] * r3,
            X[1, 3] * r1 + X[2, 3] * r2 + X[3, 3] * r3)
end

"""
    leb_order(corners::NTuple{3})

Rotate a corner triple so its longest edge connects the first and third
corner — the edge [`update_keys!`](@ref) bisects. Cyclic, so orientation is
preserved.
"""
function leb_order(c::NTuple{3})
    l12 = LinearAlgebra.norm(c[2] - c[1])
    l23 = LinearAlgebra.norm(c[3] - c[2])
    l31 = LinearAlgebra.norm(c[1] - c[3])
    if l31 >= l12 && l31 >= l23
        return c
    elseif l12 >= l23
        return (c[2], c[3], c[1])
    else
        return (c[3], c[1], c[2])
    end
end

#######
# LoD #
#######

abstract type AbstractLoD end

"""
    excess_levels(lod, base, key) -> Float64

The extension point of a level-of-detail criterion: how many bisection levels
*below* `key` does the target resolution lie? Positive means the key is too
coarse (split), non-positive for the key's *parent* means the parent is
already fine enough (its children merge). For crack-free results the value
must be a symmetric function of the key's split edge, varying by less than
one level between adjacent triangles (see the file header).
"""
function excess_levels end

"""
    UniformLoD(target)

Refine everything to exactly `target` bisection levels. Deterministic;
useful for testing and as a static-refinement fallback.
"""
struct UniformLoD <: AbstractLoD
    target::Int
end

excess_levels(lod::UniformLoD, base::IsubdBase, k::UInt64) = Float64(lod.target - key_depth(k))

"""
    ScreenSpaceLoD(projectionview, eyeposition, resolution, px_target)

Split until a triangle's split edge occupies at most `px_target` pixels on
screen. The edge is measured *isotropically*: a view-facing segment of the
edge's physical length, placed at the edge midpoint, is projected — not the
edge itself. Projecting the actual edge would make the criterion
orientation-dependent under perspective (a foreshortened edge measures
shorter), and neighbouring triangles would disagree by several levels,
opening cracks; this is the same reason the isubd demo drives its LoD from
the distance to the edge midpoint. `projectionview` is a 4×4 camera matrix
(column-major linear indexing, e.g. a `Makie.Mat4` or plain `Matrix`),
`eyeposition` the camera position, `resolution` the viewport in pixels. Each
bisection roughly halves the split edge, so the excess level count is
`log2(pixel_length / px_target)`.
"""
struct ScreenSpaceLoD{M,E} <: AbstractLoD
    projectionview::M
    eyeposition::E
    resolution::NTuple{2,Float64}
    px_target::Float64
end

_xyz(p) = length(p) >= 3 ? Float64.((p[1], p[2], p[3])) : (Float64(p[1]), Float64(p[2]), 0.0)

function _project_px(lod::ScreenSpaceLoD, p::NTuple{3,Float64})
    x, y, z = p
    m = lod.projectionview
    cx = m[1] * x + m[5] * y + m[9] * z + m[13]
    cy = m[2] * x + m[6] * y + m[10] * z + m[14]
    cw = m[4] * x + m[8] * y + m[12] * z + m[16]
    # behind-camera clamp: degenerates to a huge pixel length, i.e. "split",
    # which max_depth caps — matches the demo's behaviour of over-refining
    # rather than dropping geometry near the eye
    w = max(cw, 1e-8)
    return ((cx / w + 1.0) / 2.0 * lod.resolution[1], (cy / w + 1.0) / 2.0 * lod.resolution[2])
end

function excess_levels(lod::ScreenSpaceLoD, base::IsubdBase, k::UInt64)
    corners = key_corners(base, k)
    a = _xyz(base.mapping(key_base(k), corners[1]))
    b = _xyz(base.mapping(key_base(k), corners[3]))
    m = (a .+ b) ./ 2
    len = sqrt(sum(abs2, a .- b))
    # view direction at the midpoint, and any unit vector perpendicular to it
    eye = _xyz(lod.eyeposition)
    v = m .- eye
    nv = sqrt(sum(abs2, v))
    v = nv > 1e-12 ? v ./ nv : (0.0, 0.0, 1.0)
    u = abs(v[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
    e = (v[2] * u[3] - v[3] * u[2], v[3] * u[1] - v[1] * u[3], v[1] * u[2] - v[2] * u[1])
    e = e ./ sqrt(sum(abs2, e))
    pa = _project_px(lod, m .- e .* (len / 2))
    pb = _project_px(lod, m .+ e .* (len / 2))
    px = sqrt((pa[1] - pb[1])^2 + (pa[2] - pb[2])^2)
    return log2(max(px, 1e-9) / lod.px_target)
end

##########
# Passes #
##########

"""
    update_keys!(out, keys, base, lod; max_depth=LEB_MAX_DEPTH, hysteresis=0.0)

One split/merge/keep streaming pass over the key buffer (the GPU compute
pass of the demo): every key either emits its two children (its
[`excess_levels`](@ref) is positive), its parent (the *parent's* excess is
at most `-hysteresis`; only child 0 emits it, so the pair collapses to one
key), or itself. Refinement moves at most one level per pass — drive it with
[`refine_keys!`](@ref) to reach the steady state.

Split (`excess(key) > 0`) and merge (`excess(parent) ≤ -hysteresis`) test
disjoint predicates even at `hysteresis = 0` — a parent that just split has
positive excess, so its children never immediately merge back; the
steady state is a true fixed point. A positive `hysteresis` additionally
keeps keys whose parent hovers around excess 0 from toggling under camera
jitter, at the price of the merged state lagging the split state by up to
that many levels.
"""
function update_keys!(out::Vector{UInt64}, keys::Vector{UInt64}, base::IsubdBase, lod::AbstractLoD;
                      max_depth::Int=LEB_MAX_DEPTH, hysteresis::Float64=0.0)
    max_depth <= LEB_MAX_DEPTH || throw(ArgumentError("max_depth must be ≤ $LEB_MAX_DEPTH"))
    empty!(out)
    for k in keys
        d = key_depth(k)
        if d < max_depth && excess_levels(lod, base, k) > 0
            c0, c1 = key_children(k)
            push!(out, c0, c1)
        elseif d > 0 && excess_levels(lod, base, key_parent(k)) <= -hysteresis
            # both children evaluate the same parent predicate, so exactly one
            # of them (child 0) re-emits the parent and the other vanishes
            is_child0(k) && push!(out, key_parent(k))
        else
            push!(out, k)
        end
    end
    return out
end

"""
    refine_keys!(keys, scratch, base, lod; max_depth, hysteresis, max_passes) -> keys

Run [`update_keys!`](@ref) passes until the key buffer is a fixed point of
the criterion (or `max_passes` is hit — one more than `max_depth` suffices
for any monotone criterion). `keys` is updated in place, `scratch` is the
ping-pong buffer.
"""
function refine_keys!(keys::Vector{UInt64}, scratch::Vector{UInt64}, base::IsubdBase, lod::AbstractLoD;
                      max_depth::Int=LEB_MAX_DEPTH, hysteresis::Float64=0.0,
                      max_passes::Int=max_depth + 1)
    sort!(keys)
    for _ in 1:max_passes
        update_keys!(scratch, keys, base, lod; max_depth, hysteresis)
        sort!(scratch)
        scratch == keys && break
        copy!(keys, scratch)
    end
    return keys
end

"""
    decode_keys!(mesh::IsubdMesh, keys, base) -> mesh

Emit the key buffer as a triangle soup: three fresh vertices per key (matching
`FEData`'s duplicated-vertex layout, so discontinuous data stays representable)
with their reference coordinates, mapped through `base.mapping`. Buffers are
resized in place, so a steady key count re-renders without allocating. Each
bisection flips the triangle's handedness, so odd depths emit the face
reversed to keep a consistent winding.
"""
function decode_keys!(mesh, keys::Vector{UInt64}, base::IsubdBase)
    n = length(keys)
    resize!(mesh.positions, 3n)
    resize!(mesh.refcoords, 3n)
    resize!(mesh.faces, n)
    for (i, k) in enumerate(keys)
        c1, c2, c3 = key_corners(base, k)
        b = key_base(k)
        j = 3 * (i - 1)
        mesh.refcoords[j + 1] = c1
        mesh.refcoords[j + 2] = c2
        mesh.refcoords[j + 3] = c3
        mesh.positions[j + 1] = base.mapping(b, c1)
        mesh.positions[j + 2] = base.mapping(b, c2)
        mesh.positions[j + 3] = base.mapping(b, c3)
        mesh.faces[i] = isodd(key_depth(k)) ? Int32.((j + 3, j + 2, j + 1)) :
                        Int32.((j + 1, j + 2, j + 3))
    end
    return mesh
end

"""
    IsubdMesh(base::IsubdBase)

The reusable output buffers of [`decode_keys!`](@ref): vertex positions (in
the mapping's image space), per-vertex reference coordinates, and one face
per key.
"""
struct IsubdMesh{P,RV}
    positions::Vector{P}
    refcoords::Vector{RV}
    faces::Vector{NTuple{3,Int32}}
end

function IsubdMesh(base::IsubdBase{RV}) where {RV}
    P = typeof(base.mapping(1, base.corners[1][1]))
    return IsubdMesh(P[], RV[], NTuple{3,Int32}[])
end
