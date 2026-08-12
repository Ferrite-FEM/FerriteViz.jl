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

# Edge types of a triangle (c1, c2, c3), in the corner convention above.
const EDGE_S = 1   # (c1, c3) — the split edge
const EDGE_L = 2   # (c1, c2)
const EDGE_R = 3   # (c2, c3)

# One entry of the base adjacency table: which base triangle sits across an
# edge, which of *its* edges that is, and whether the two traverse the shared
# segment in opposite directions. `base == 0` marks a boundary edge.
const BaseNeighbor = Tuple{Int,Int,Bool}
const NO_NEIGHBOR = (0, 0, false)

"""
    IsubdBase(corners, mapping[, adjacency])

The base domain of an implicit LEB subdivision: one reference-space corner
triple per base triangle, ordered for bisection (see [`leb_order`](@ref)),
plus the geometry `mapping(base_id, ξ) -> physical point` (for FEM cells the
geometric map, optionally composed with a warp field). The mapping is applied
per decoded vertex, so curved cells subdivide into curved sub-triangles.

`adjacency[b][e]` names the base triangle across edge `e` of base triangle `b`
(`EDGE_S`/`EDGE_L`/`EDGE_R`) as `(base, edge, reversed)`, or `NO_NEIGHBOR` for
a boundary. Supplying it enables conforming (watertight) refinement — see
[`refine_keys!`](@ref). It must be *compatible*: a split edge may only pair
with another split edge, and a leg only with legs. `_isubd_base` builds such a
table for `FEData`; pairings that violate compatibility are dropped to
boundaries, which costs conformity along those edges but nothing else.
"""
struct IsubdBase{RV,F}
    corners::Vector{NTuple{3,RV}}
    mapping::F
    adjacency::Vector{NTuple{3,BaseNeighbor}}
end

IsubdBase(corners, mapping) =
    IsubdBase(corners, mapping, [(NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR) for _ in corners])

is_conformable(base::IsubdBase) = any(t -> any(e -> e[1] != 0, t), base.adjacency)

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

# Where a triangle's linear interpolation is checked against the truth, in
# barycentric coordinates. The edge samples bound the deviation of the *drawn
# edges*, which is what a gap or colour seam between refinement levels would
# be; the interior samples bound what is actually rendered across the face,
# which for a smooth bump is where the deviation peaks. Measuring the interior
# is only sound because conformity is enforced structurally (see
# `force_split!`) — while split decisions had to agree between the two
# triangles sharing an edge, only shared-edge quantities could be used.
# Where the linear interpolation is checked: the midpoint of every edge — the
# vertices a full subdivision would introduce, and the worst points of a
# smooth function along them — plus the centroid for deviation that lives in
# the face rather than on its boundary.
#
# This criterion is *the* hot path (every leaf, every update, an FE evaluation
# per sample), so the temptation is to keep only the split edge and the
# centroid. Don't: a multilinear field is *exactly linear along element
# edges* and curves only across the fan diagonals, which are the legs. On a
# trilinear field over hexahedra that reduction stopped the refinement
# altogether — zero deviation on every edge it still looked at. The quarter
# points, on the other hand, add nothing measurable and are not sampled.
const DEVIATION_SAMPLES = (
    (0.5, 0.0, 0.5), (0.5, 0.5, 0.0), (0.0, 0.5, 0.5), (1 / 3, 1 / 3, 1 / 3),
)

"""
    DeviationLoD(f, tol)

Split until the triangle's linear interpolation approximates `f(base_id, ξ)`
to within `tol`. The deviation is sampled over the whole triangle — along
every edge and across the interior (see `DEVIATION_SAMPLES`) — against the
barycentric interpolation of the corner values, and `excess_levels` is
`log2(deviation / tol)`: the deviation of a smooth function under linear
interpolation is O(h²) and bisection halves an edge every *second* level, so
each level buys a factor 2. `f` may return points (geometry error: pass the
base's `mapping`) or scalars (solution error: pass the colour evaluation);
`tol` is absolute, in the units of `norm` of `f`'s values.

Sampling the interior is what makes the criterion bound what is actually
drawn — for curved geometry the deviation peaks in the middle of a face, and
an edge-only criterion happily leaves it there. Keeping the edge samples as
well bounds the width of any gap or colour seam at a refinement-level
boundary, which matters when conformity is disabled.

The criterion is not monotone in depth (a triangle's deviation can vanish
while a descendant's does not — e.g. a bilinear field is linear along a
quad's outer edges but curved along the fan diagonals), which makes the
refined state mildly path-dependent: a state merged down from finer keys may
stay finer than one refined up from the roots, because the passes never
discard detail whose deviation still exceeds the tolerance. The finer of the
two states is the more accurate one.
"""
struct DeviationLoD{F,T} <: AbstractLoD
    f::F
    tol::T
end

function excess_levels(lod::DeviationLoD, base::IsubdBase, k::UInt64)
    c1, c2, c3 = key_corners(base, k)
    b = key_base(k)
    f1, f2, f3 = lod.f(b, c1), lod.f(b, c2), lod.f(b, c3)
    err = 0.0
    for (a1, a2, a3) in DEVIATION_SAMPLES
        exact = lod.f(b, a1 * c1 + a2 * c2 + a3 * c3)
        linear = a1 * f1 + a2 * f2 + a3 * f3
        err = max(err, Float64(LinearAlgebra.norm(exact - linear)))
    end
    return log2(max(err, 1e-16) / max(lod.tol, 1e-16))
end

"""
    CachedLoD(inner)

Memoize a criterion per key. For a fixed solution the excess is a pure
function of the key, but the passes ask for the same keys again and again —
every refinement round re-tests the surviving leaves, and the conforming
closure additionally asks about parents and neighbours. With an FE evaluation
behind every query (a `PointValues` reinit per sample point) that repetition
dominates; [`refine_keys!`](@ref) therefore wraps its criterion in this for
the duration of the call.
"""
struct CachedLoD{L<:AbstractLoD} <: AbstractLoD
    inner::L
    cache::Dict{UInt64,Float64}
end
CachedLoD(inner::AbstractLoD) = CachedLoD(inner, Dict{UInt64,Float64}())

excess_levels(lod::CachedLoD, base::IsubdBase, k::UInt64) =
    get!(() -> excess_levels(lod.inner, base, k), lod.cache, k)

"""
    CombinedLoD(lods...)

Split when *any* member criterion wants to: the excess is the member maximum.
"""
struct CombinedLoD{T<:Tuple} <: AbstractLoD
    lods::T
end
CombinedLoD(lods::AbstractLoD...) = CombinedLoD(lods)

excess_levels(lod::CombinedLoD, base::IsubdBase, k::UInt64) =
    maximum(l -> excess_levels(l, base, k), lod.lods)

##########
# Passes #
##########

"""
    update_keys!(out, keys, base, lod; max_depth=LEB_MAX_DEPTH, hysteresis=0.0)

One split/merge/keep streaming pass over the *sorted* key buffer (the GPU
compute pass of the demo): every key either emits its two children (its
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

A pair merges only when *both* siblings are present as leaves, which the
sorted order makes an adjacent-element check. The demo omits this and relies
on its LoD never jumping levels between siblings; under a criterion with
sharp spatial variation (an error estimator on rough data), the unguarded
merge lets a parent overlap its sibling's still-deeper subtree, or drops a
child while the sibling subtree persists — converging to a state that
double-covers or holes the domain. (A GPU port checks sibling presence on the
concurrent-binary-tree bitfield instead of the sorted buffer.)
"""
function update_keys!(out::Vector{UInt64}, keys::Vector{UInt64}, base::IsubdBase, lod::AbstractLoD;
                      max_depth::Int=LEB_MAX_DEPTH, hysteresis::Float64=0.0)
    max_depth <= LEB_MAX_DEPTH || throw(ArgumentError("max_depth must be ≤ $LEB_MAX_DEPTH"))
    issorted(keys) || throw(ArgumentError("update_keys! requires a sorted key buffer"))
    empty!(out)
    for (i, k) in enumerate(keys)
        d = key_depth(k)
        if d < max_depth && excess_levels(lod, base, k) > 0
            c0, c1 = key_children(k)
            push!(out, c0, c1)
        elseif d > 0 && excess_levels(lod, base, key_parent(k)) <= -hysteresis
            # the pair collapses only when the sibling is a leaf too and does
            # not itself want splitting — evaluated symmetrically from both
            # sides, so child 0 emits the parent exactly when child 1 drops
            sib_idx = is_child0(k) ? i + 1 : i - 1
            sib = is_child0(k) ? k + 1 : k - 1
            mergeable = checkbounds(Bool, keys, sib_idx) && keys[sib_idx] == sib &&
                        !(d < max_depth && excess_levels(lod, base, sib) > 0)
            if !mergeable
                push!(out, k)
            elseif is_child0(k)
                push!(out, key_parent(k))
            end
        else
            push!(out, k)
        end
    end
    return out
end

###############################
# Neighbours and conformity   #
###############################

# Splitting a triangle alone subdivides only its split edge, so it is the only
# edge that can end up half-drawn: a triangle's legs become the *whole* split
# edges of its children, and stay full edges of whatever leaf lies across them.
# Conformity therefore reduces to one rule — never split a triangle without
# simultaneously splitting the leaf across its split edge ("diamond partner"),
# forcing that partner down first when it is coarser. This is classical
# newest-vertex bisection; the demo skips it because a smooth camera criterion
# keeps neighbouring levels within one of each other, which an error estimator
# on rough data does not.
#
# `key_neighbour` answers "which triangle at the same level lies across edge e"
# purely from the key and the base table, by recursing to the root. Its
# well-definedness rests on the compatibility invariant (split edges pair with
# split edges, legs with legs), which the recursion preserves at every level:
#
#   child 0 = (v1, m, v2)   split = parent's left leg
#   child 1 = (v2, m, v3)   split = parent's right leg
#   both     legs           = halves of the parent's split edge, or the new
#                             interior edge shared by the two children
#
# A GPU port replaces the base table lookup with the same recursion over a
# concurrent binary tree; nothing here needs the leaf set to be materialized.
function key_neighbour(base::IsubdBase, k::UInt64, e::Int)
    if key_depth(k) == 0
        b, e2, rev = base.adjacency[key_base(k)][e]
        return b == 0 ? nothing : (make_key(b, UInt32(1)), e2, rev)
    end
    p = key_parent(k)
    first = is_child0(k)
    if e == EDGE_S
        # the child's split edge is one of the parent's legs
        nb = key_neighbour(base, p, first ? EDGE_L : EDGE_R)
        nb === nothing && return nothing
        q, y, rev = nb
        y == EDGE_S && return nothing        # incompatible labelling: treat as boundary
        c0, c1 = key_children(q)
        return (y == EDGE_L ? c0 : c1, EDGE_S, rev)
    elseif e == EDGE_L
        first || return (key_children(p)[1], EDGE_R, true)   # interior edge, shared with the sibling
        nb = key_neighbour(base, p, EDGE_S)                  # half of the parent's split edge
        nb === nothing && return nothing
        r, y, rev = nb
        y == EDGE_S || return nothing
        c0, c1 = key_children(r)
        return rev ? (c1, EDGE_R, true) : (c0, EDGE_L, false)
    else # EDGE_R
        first && return (key_children(p)[2], EDGE_L, true)   # interior edge, shared with the sibling
        nb = key_neighbour(base, p, EDGE_S)
        nb === nothing && return nothing
        r, y, rev = nb
        y == EDGE_S || return nothing
        c0, c1 = key_children(r)
        return rev ? (c0, EDGE_L, true) : (c1, EDGE_R, false)
    end
end

"""
    diamond_partner(base, key) -> UInt64 | nothing

The triangle at `key`'s own level across its split edge — the only leaf that
has to split together with `key` to keep the mesh conforming.
"""
function diamond_partner(base::IsubdBase, k::UInt64)
    nb = key_neighbour(base, k, EDGE_S)
    return nb === nothing ? nothing : nb[1]
end

# The leaf covering `k`'s region when `k` itself is not one: its closest
# ancestor in the leaf set. `nothing` means the region is covered by leaves
# *finer* than `k`.
function leaf_ancestor(leaves, k::UInt64)
    while key_depth(k) > 0
        k = key_parent(k)
        k in leaves && return k
    end
    return nothing
end

function _split_leaf!(leaves::Set{UInt64}, k::UInt64)
    delete!(leaves, k)
    c0, c1 = key_children(k)
    push!(leaves, c0)
    push!(leaves, c1)
    return leaves
end

"""
    force_split!(leaves, base, key; max_depth)

Split `key` and everything that has to split with it to keep the mesh
conforming: the leaf across its split edge, recursively forced down when it is
coarser. Terminates because each forced neighbour is strictly coarser than the
triangle that asked for it.
"""
function force_split!(leaves::Set{UInt64}, base::IsubdBase, k::UInt64; max_depth::Int=LEB_MAX_DEPTH)
    for _ in 0:max_depth   # each retry resolves the neighbour one level deeper
        (k in leaves && key_depth(k) < max_depth) || return leaves
        n = diamond_partner(base, k)
        if n === nothing || n in leaves
            _split_leaf!(leaves, k)
            n === nothing || _split_leaf!(leaves, n)
            return leaves
        end
        a = leaf_ancestor(leaves, n)
        if a === nothing
            # the far side is already finer than us: our split edge is whole on
            # this side, so splitting alone keeps the mesh conforming
            _split_leaf!(leaves, k)
            return leaves
        end
        force_split!(leaves, base, a; max_depth)
    end
    return leaves
end

"""
    conforming_update!(leaves, base, lod; max_depth, hysteresis) -> changed::Bool

One conforming split/merge pass over the leaf set: every leaf whose
[`excess_levels`](@ref) is positive is split through [`force_split!`](@ref),
and every diamond of four leaves whose two parents both want to coarsen is
merged back. Both operations move whole diamonds, which is what keeps every
drawn edge a full edge of the leaf on the other side — no T-vertices, and
therefore no slivers or colour seams between refinement levels, on curved
geometry as much as on flat.
"""
function conforming_update!(leaves::Set{UInt64}, base::IsubdBase, lod::AbstractLoD;
                            max_depth::Int=LEB_MAX_DEPTH, hysteresis::Float64=0.0)
    max_depth <= LEB_MAX_DEPTH || throw(ArgumentError("max_depth must be ≤ $LEB_MAX_DEPTH"))
    before = length(leaves)
    wanted = UInt64[]
    for k in leaves
        key_depth(k) < max_depth && excess_levels(lod, base, k) > 0 && push!(wanted, k)
    end
    sort!(wanted)   # deterministic order: the Set's iteration order is not
    for k in wanted
        force_split!(leaves, base, k; max_depth)
    end
    split_count = length(leaves)

    # Merges: a diamond collapses only as a whole. Both parents are tested, so
    # the decision is symmetric and each parent is emitted once, by its child 0.
    mergeable = UInt64[]
    for k in leaves
        (is_child0(k) && key_depth(k) > 0) || continue
        (k + 1) in leaves || continue
        p = key_parent(k)
        excess_levels(lod, base, p) <= -hysteresis || continue
        r = diamond_partner(base, p)
        if r !== nothing
            c0, c1 = key_children(r)
            (c0 in leaves && c1 in leaves) || continue
            excess_levels(lod, base, r) <= -hysteresis || continue
        end
        push!(mergeable, p)
    end
    for p in mergeable
        c0, c1 = key_children(p)
        (c0 in leaves && c1 in leaves) || continue
        delete!(leaves, c0)
        delete!(leaves, c1)
        push!(leaves, p)
    end
    return length(leaves) != before || split_count != before
end

"""
    refine_keys!(keys, scratch, base, lod; max_depth, hysteresis, max_passes, conforming) -> keys

Run passes until the key buffer is a fixed point of the criterion (or
`max_passes` is hit — one more than `max_depth` suffices for any monotone
criterion). `keys` is updated in place, `scratch` is the ping-pong buffer.

With `conforming = true` (the default whenever `base` carries an adjacency
table) the passes are [`conforming_update!`](@ref) and the result is
watertight; otherwise they are the streaming [`update_keys!`](@ref), whose
refinement-level boundaries leave gaps bounded by the criterion's tolerance.
"""
function refine_keys!(keys::Vector{UInt64}, scratch::Vector{UInt64}, base::IsubdBase, lod::AbstractLoD;
                      max_depth::Int=LEB_MAX_DEPTH, hysteresis::Float64=0.0,
                      max_passes::Int=max_depth + 1, conforming::Bool=is_conformable(base))
    lod = lod isa CachedLoD ? lod : CachedLoD(lod)
    if conforming
        leaves = Set(keys)
        for _ in 1:max_passes
            conforming_update!(leaves, base, lod; max_depth, hysteresis) || break
        end
        resize!(keys, length(leaves))
        copyto!(keys, collect(leaves))
        sort!(keys)
        return keys
    end
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
    decode_keys!(mesh::IsubdMesh, keys, base; groups=nothing) -> mesh

Emit the key buffer as a drawable mesh: positions, their reference
coordinates, and one face per key, with each bisection's handedness flip
undone so the winding stays consistent.

Vertices are shared within a `groups` class and duplicated across classes.
`FEData` passes the owning cell, which shares everything inside a cell —
where the drawn field is continuous — while keeping element boundaries
duplicated, so discontinuous (L2/DG) fields keep their jumps exactly as in
the static tessellation. A conforming mesh has about half as many vertices as
triangles, against three per triangle unshared, and every vertex costs a
geometry evaluation here and a field evaluation downstream. With
`groups=nothing` each base triangle forms its own class, sharing within its
own subtree and duplicating along base edges.

Buffers (including the lookup table) are emptied rather than reallocated, so
re-decoding a steady mesh does not grow the heap.
"""
function decode_keys!(mesh, keys::Vector{UInt64}, base::IsubdBase; groups=nothing)
    decode_topology!(mesh, keys, base; groups)
    return decode_positions!(mesh, base)
end

"""
    decode_topology!(mesh, keys, base; groups=nothing) -> mesh

The connectivity half of [`decode_keys!`](@ref): reference coordinates, their
owning base triangle, and the faces. This is the part that depends only on the
key set, so a consumer whose mesh is unchanged can re-run just
[`decode_positions!`](@ref) and skip the vertex-sharing lookups entirely.
"""
function decode_topology!(mesh, keys::Vector{UInt64}, base::IsubdBase; groups=nothing)
    empty!(mesh.refcoords)
    empty!(mesh.vertex_base)
    empty!(mesh.lut)
    resize!(mesh.faces, length(keys))
    for (i, k) in enumerate(keys)
        c1, c2, c3 = key_corners(base, k)
        b = key_base(k)
        g = Int32(groups === nothing ? b : groups[b])
        i1 = _vertex!(mesh, b, g, c1)
        i2 = _vertex!(mesh, b, g, c2)
        i3 = _vertex!(mesh, b, g, c3)
        mesh.faces[i] = isodd(key_depth(k)) ? (i3, i2, i1) : (i1, i2, i3)
    end
    return mesh
end

"""
    decode_positions!(mesh, base) -> mesh

Map the vertices laid out by [`decode_topology!`](@ref) through
`base.mapping`. Cheap to repeat: it touches one point per vertex and no
lookup table, which is what an unchanged mesh under a changing solution
needs.
"""
function decode_positions!(mesh, base::IsubdBase)
    resize!(mesh.positions, length(mesh.refcoords))
    @inbounds for v in eachindex(mesh.refcoords)
        mesh.positions[v] = base.mapping(Int(mesh.vertex_base[v]), mesh.refcoords[v])
    end
    return mesh
end

@inline function _vertex!(mesh, b::Int, g::Int32, ξ)
    # reference coordinates are dyadic (corner averages), so they compare
    # exactly — no rounding, no tolerance
    idx = get(mesh.lut, (g, ξ), Int32(0))
    idx == 0 || return idx
    push!(mesh.refcoords, ξ)
    push!(mesh.vertex_base, Int32(b))
    new_idx = Int32(length(mesh.refcoords))
    mesh.lut[(g, ξ)] = new_idx
    return new_idx
end

"""
    IsubdMesh(base::IsubdBase)

The reusable output buffers of [`decode_keys!`](@ref): vertex positions (in
the mapping's image space), per-vertex reference coordinates and owning base
triangle, one face per key, and the lookup table backing vertex sharing.
"""
struct IsubdMesh{P,RV}
    positions::Vector{P}
    refcoords::Vector{RV}
    vertex_base::Vector{Int32}
    faces::Vector{NTuple{3,Int32}}
    lut::Dict{Tuple{Int32,RV},Int32}
end

function IsubdMesh(base::IsubdBase{RV}) where {RV}
    P = typeof(base.mapping(1, base.corners[1][1]))
    return IsubdMesh(P[], RV[], Int32[], NTuple{3,Int32}[], Dict{Tuple{Int32,RV},Int32}())
end
