# Tests for the CPU isubd core (src/isubd.jl): LEB key codec, bisection
# geometry, split/merge passes, LoD criteria, crack-freeness, buffer reuse.
using Tensors: Tensors, Vec, Tensor, norm
import FerriteViz as FV

# The demo's base domain: a unit square as two triangles sharing their split
# edge (the diagonal) — the "diamond" pairing the crack-freeness argument
# rests on. Corners ordered by leb_order convention: split edge = (v1, v3).
function diamond_base(mapping = (b, ξ) -> ξ)
    t1 = FV.leb_order((Vec((0.0, 0.0)), Vec((1.0, 0.0)), Vec((1.0, 1.0))))
    t2 = FV.leb_order((Vec((1.0, 1.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0))))
    return FV.IsubdBase([t1, t2], mapping)
end

_area(a, b, c) = 0.5 * abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2]))
_signed_area(a, b, c) = 0.5 * ((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2]))

# A drawn edge is cracked when some vertex of the mesh sits exactly at its
# midpoint: the neighbour across the edge is finer and its corner opens a
# T-junction on our straight edge. Exact only where the mapping is affine per
# cell (curved cells displace the mapped reference midpoint off the straight
# edge by O(h²)) — run it on linear geometry.
# `+ 0.0` normalizes negative zero: -0.0 and 0.0 are the same point but not
# the same dictionary key, and evaluation paths differ in which one they
# produce.
function count_tjunctions(mesh; digits = 9)
    key(p) = (round(p[1]; digits = digits) + 0.0, round(p[2]; digits = digits) + 0.0)
    verts = Set(key(p) for p in mesh.positions)
    cracks = 0
    for f in mesh.faces
        for (i, j) in ((1, 2), (2, 3), (3, 1))
            a, b = mesh.positions[f[i]], mesh.positions[f[j]]
            mid = key((a + b) / 2)
            if mid in verts && mid != key(a) && mid != key(b)
                cracks += 1
            end
        end
    end
    return cracks
end

# Edge bookkeeping on the unit-square diamond: how many interior edges are
# drawn only once (a hole), how many sit on the domain boundary, and how many
# are shared by more than two triangles (an overlap).
function edge_multiplicities(mesh; digits = 9)
    key(p) = (round(p[1]; digits = digits) + 0.0, round(p[2]; digits = digits) + 0.0)
    counts = Dict{Tuple{Any,Any},Int}()
    for f in mesh.faces, (i, j) in ((1, 2), (2, 3), (3, 1))
        a, b = key(mesh.positions[f[i]]), key(mesh.positions[f[j]])
        counts[a <= b ? (a, b) : (b, a)] = get(counts, a <= b ? (a, b) : (b, a), 0) + 1
    end
    border(p) = isapprox(p[1], 0; atol = 1e-9) || isapprox(p[1], 1; atol = 1e-9) ||
                isapprox(p[2], 0; atol = 1e-9) || isapprox(p[2], 1; atol = 1e-9)
    interior_once = boundary_once = over = 0
    for ((a, b), c) in counts
        if c == 1
            (border(a) && border(b)) ? (boundary_once += 1) : (interior_once += 1)
        elseif c > 2
            over += 1
        end
    end
    return interior_once, boundary_once, over
end

@testset "isubd key codec" begin
    k = FV.root_key(7)
    @test FV.key_base(k) == 7
    @test FV.key_depth(k) == 0
    c0, c1 = FV.key_children(k)
    @test FV.key_depth(c0) == FV.key_depth(c1) == 1
    @test FV.key_parent(c0) == k && FV.key_parent(c1) == k
    @test FV.is_child0(c0) && !FV.is_child0(c1)
    @test FV.key_base(c0) == FV.key_base(c1) == 7
    # deep path round-trips
    key = FV.root_key(42)
    path = [0, 1, 1, 0, 1, 0, 0, 1]
    for b in path
        key = FV.key_children(key)[b + 1]
    end
    @test FV.key_depth(key) == length(path)
    for b in reverse(path)
        @test FV.key_children(FV.key_parent(key))[b + 1] == key
        key = FV.key_parent(key)
    end
    @test key == FV.root_key(42)
    # root xform is the identity
    @test FV.key_xform(FV.root_key(1)) ≈ one(Tensors.Tensor{2,3,Float64})
end

@testset "isubd bisection geometry" begin
    base = diamond_base()
    root = FV.root_key(1)
    v = FV.key_corners(base, root)
    c0, c1 = FV.key_children(root)
    w0 = FV.key_corners(base, c0)
    w1 = FV.key_corners(base, c1)
    m = (v[1] + v[3]) / 2
    # children = (v1, m, v2) and (v2, m, v3): exact tiling of the parent
    @test w0[1] ≈ v[1] && w0[2] ≈ m && w0[3] ≈ v[2]
    @test w1[1] ≈ v[2] && w1[2] ≈ m && w1[3] ≈ v[3]
    @test _area(w0...) + _area(w1...) ≈ _area(v...)
    # both children share the parent's handedness (flipped relative to parent)
    @test sign(_signed_area(w0...)) == sign(_signed_area(w1...)) == -sign(_signed_area(v...))
    # the child's split edge (corner 1 ↔ 3) is its longest edge, at every depth
    key = root
    for _ in 1:6
        key = FV.key_children(key)[1]
        c = FV.key_corners(base, key)
        l13 = norm(c[1] - c[3])
        @test l13 >= norm(c[1] - c[2]) - 1e-12 && l13 >= norm(c[2] - c[3]) - 1e-12
    end
end

@testset "isubd leb_order" begin
    tri = (Vec((0.0, 0.0)), Vec((1.0, 0.0)), Vec((1.0, 1.0)))
    for rot in 0:2
        rotated = ntuple(i -> tri[mod1(i + rot, 3)], 3)
        ordered = FV.leb_order(rotated)
        @test norm(ordered[1] - ordered[3]) >= norm(ordered[1] - ordered[2]) - 1e-12
        @test norm(ordered[1] - ordered[3]) >= norm(ordered[2] - ordered[3]) - 1e-12
        # cyclic rotation only: orientation preserved
        @test sign(_signed_area(ordered...)) == sign(_signed_area(rotated...))
    end
end

@testset "isubd uniform refinement and merge" begin
    base = diamond_base()
    keys, scratch = FV.root_keys(base), UInt64[]
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(3))
    @test length(keys) == 2 * 2^3
    @test all(k -> FV.key_depth(k) == 3, keys)
    mesh = FV.IsubdMesh(base)
    FV.decode_keys!(mesh, keys, base)
    @test sum(_area(mesh.positions[f[1]], mesh.positions[f[2]], mesh.positions[f[3]])
              for f in mesh.faces) ≈ 1.0     # tiles the unit square exactly
    @test count_tjunctions(mesh) == 0
    # a steady state is a fixed point of the pass
    before = sort(copy(keys))
    FV.update_keys!(scratch, keys, base, FV.UniformLoD(3))
    @test sort(scratch) == before
    # relaxing the criterion merges all the way back to the roots
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(1))
    @test sort(keys) == sort([FV.key_children(FV.root_key(1))..., FV.key_children(FV.root_key(2))...])
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(0))
    @test sort(keys) == sort(FV.root_keys(base))
    # max_depth caps the split branch
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(9); max_depth = 4)
    @test all(k -> FV.key_depth(k) == 4, keys)
end

# A deliberately non-smooth criterion: deep refinement left of x = 0.7, none
# to the right — level jumps far beyond 1 between neighbouring triangles.
struct StepLoD <: FV.AbstractLoD
    target::Int
end
function FV.excess_levels(lod::StepLoD, base::FV.IsubdBase, k::UInt64)
    c = FV.key_corners(base, k)
    mid = (c[1] + c[3]) / 2
    return mid[1] < 0.7 ? Float64(lod.target - FV.key_depth(k)) : Float64(-FV.key_depth(k))
end

@testset "isubd deviation and combined criteria" begin
    # paraboloid: linear interpolation error is O(h²), nonzero everywhere
    bump = (b, ξ) -> Vec((ξ[1], ξ[2], ξ[1] * (1 - ξ[1]) + ξ[2] * (1 - ξ[2])))
    base = diamond_base(bump)
    keys, scratch = FV.root_keys(base), UInt64[]

    # a linear (flat) mapping never asks for refinement
    flat = diamond_base((b, ξ) -> Vec((ξ[1], ξ[2], 0.3 * ξ[1] + 0.7 * ξ[2])))
    kflat = FV.root_keys(flat)
    FV.refine_keys!(kflat, scratch, flat, FV.DeviationLoD(flat.mapping, 1e-6))
    @test sort(kflat) == sort(FV.root_keys(flat))

    # the deviation is O(h²) and the split edge halves every second level, so
    # a 16× tighter tolerance buys about four more levels
    FV.refine_keys!(keys, scratch, base, FV.DeviationLoD(bump, 1e-3))
    d1 = maximum(FV.key_depth, keys)
    n1 = length(keys)
    FV.refine_keys!(keys, scratch, base, FV.DeviationLoD(bump, 1e-3 / 16))
    d2 = maximum(FV.key_depth, keys)
    @test n1 > 2 && d2 - d1 in 3:5
    # at the steady state, every leaf's split-edge deviation is within tolerance
    for k in keys
        c = FV.key_corners(base, k)
        fa, fb = bump(0, c[1]), bump(0, c[3])
        for t in (0.25, 0.5, 0.75)
            dev = norm(bump(0, c[1] + t * (c[3] - c[1])) - (fa + t * (fb - fa)))
            @test dev <= 1e-3 / 16 + 1e-12
        end
    end
    # crack-free
    mesh = FV.IsubdMesh(base)
    FV.decode_keys!(mesh, keys, base)
    @test count_tjunctions(mesh) == 0

    # combined: the maximum excess wins
    kc = FV.root_keys(base)
    FV.refine_keys!(kc, scratch, base, FV.CombinedLoD(FV.UniformLoD(2), FV.UniformLoD(0)))
    @test all(k -> FV.key_depth(k) == 2, kc)

    # On curved data, wherever leaves of different depth meet, the visible gap
    # is the coarse edge's deviation — bounded by tol because every leaf edge
    # is measured. Force unequal depths with a step criterion and check every
    # reference-space T-vertex's physical gap.
    tol = 1e-3
    kg = FV.root_keys(base)
    # the step must out-refine the deviation criterion (uniform depth 9 here)
    # somewhere, or no unequal-depth boundaries exist
    FV.refine_keys!(kg, scratch, base, FV.CombinedLoD(FV.DeviationLoD(bump, tol), StepLoD(12)); max_depth=14)
    meshg = FV.IsubdMesh(base)
    FV.decode_keys!(meshg, kg, base)
    rk(p) = (round(p[1]; digits = 12), round(p[2]; digits = 12))
    refverts = Set(rk(p) for p in meshg.refcoords)
    ntv, maxgap = 0, 0.0
    for (fi, f) in enumerate(meshg.faces), (i, j) in ((1, 2), (2, 3), (3, 1))
        a, b = meshg.refcoords[f[i]], meshg.refcoords[f[j]]
        mid = (a + b) / 2
        if rk(mid) in refverts && rk(mid) != rk(a) && rk(mid) != rk(b)
            ntv += 1
            chord = (meshg.positions[f[i]] + meshg.positions[f[j]]) / 2
            maxgap = max(maxgap, norm(bump(FV.key_base(kg[fi]), mid) - chord))
        end
    end
    @test ntv > 0                 # the step forces unequal-depth boundaries
    @test maxgap <= tol + 1e-12
end

@testset "DeviationLoD: configurable sample points" begin
    base = diamond_base()
    f0 = (b, ξ) -> 0.0
    @test FV.DeviationLoD(f0, 1e-3).samples == FV.DEVIATION_SAMPLES
    @test_throws ArgumentError FV.DeviationLoD(f0, 1e-3; samples = ())
    @test_throws ArgumentError FV.DeviationLoD(f0, 1e-3; samples = ((0.5, 0.5, 0.5),))
    @test_throws ArgumentError FV.DeviationLoD(f0, 1e-3; samples = ((1.5, -0.5, 0.0),))

    # A sampled deviation is only a lower bound: this quintic vanishes at the
    # corners, the edge midpoints and the centroid of both root triangles, so
    # the default samples miss it entirely — a denser set catches it.
    p(t) = t * (t - 1 / 3) * (t - 1 / 2) * (t - 2 / 3) * (t - 1)
    blind = (b, ξ) -> p(ξ[1])
    tol = 1e-4
    lod_default = FV.DeviationLoD(blind, tol)
    lod_dense = FV.DeviationLoD(blind, tol;
                                samples = (FV.DEVIATION_SAMPLES..., (0.7, 0.2, 0.1)))
    for k in FV.root_keys(base)
        @test FV.deviation(lod_default, base, k) == 0.0
        @test FV.deviation(lod_dense, base, k) > tol
    end
    keys, scratch = FV.root_keys(base), UInt64[]
    FV.refine_keys!(keys, scratch, base, lod_default)
    @test maximum(FV.key_depth, keys) == 0
    FV.refine_keys!(keys, scratch, base, lod_dense; max_depth = 6)
    @test maximum(FV.key_depth, keys) > 0
end

@testset "isubd merge under sharp level jumps" begin
    base = diamond_base()
    mesh = FV.IsubdMesh(base)
    refarea(keys) = (FV.decode_keys!(mesh, keys, base);
                     sum(_area(mesh.refcoords[f[1]], mesh.refcoords[f[2]], mesh.refcoords[f[3]])
                         for f in mesh.faces))
    keys, scratch = FV.root_keys(base), UInt64[]
    FV.refine_keys!(keys, scratch, base, StepLoD(6))
    # exact tiling: no overlaps, no holes, even across the sharp jump
    @test refarea(keys) ≈ 1.0
    @test maximum(FV.key_depth, keys) == 6 && minimum(FV.key_depth, keys) <= 2
    # collapsing back through wildly unequal sibling depths reaches the roots
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(0))
    @test sort(keys) == sort(FV.root_keys(base))
    # and every intermediate pass conserves the tiling
    FV.refine_keys!(keys, scratch, base, StepLoD(6))
    lod0 = FV.UniformLoD(0)
    for _ in 1:8
        FV.update_keys!(scratch, keys, base, lod0)
        sort!(scratch)
        copy!(keys, scratch)
        @test refarea(keys) ≈ 1.0
    end
    @test sort(keys) == sort(FV.root_keys(base))
end

# the same diamond, with the adjacency that makes it conformable: the two
# triangles share the diagonal as both their split edges, traversed opposite
function diamond_base_adj(mapping = (b, ξ) -> ξ)
    b = diamond_base(mapping)
    adjacency = [((2, FV.EDGE_S, true), FV.NO_NEIGHBOR, FV.NO_NEIGHBOR),
                 ((1, FV.EDGE_S, true), FV.NO_NEIGHBOR, FV.NO_NEIGHBOR)]
    return FV.IsubdBase(b.corners, b.mapping, adjacency)
end

@testset "isubd neighbour algebra" begin
    base = diamond_base_adj()
    @test FV.is_conformable(base)
    @test !FV.is_conformable(diamond_base())
    @test FV.diamond_partner(base, FV.root_key(1)) == FV.root_key(2)
    @test FV.diamond_partner(base, FV.root_key(2)) == FV.root_key(1)

    # the relation is an involution wherever it is defined, at every depth,
    # and partners always sit at the same level
    keys = FV.root_keys(base)
    for _ in 1:6
        keys = vcat((collect(FV.key_children(k)) for k in keys)...)
        for k in keys
            n = FV.diamond_partner(base, k)
            n === nothing && continue
            @test FV.key_depth(n) == FV.key_depth(k)
            @test FV.diamond_partner(base, n) == k
        end
    end

    # a partner shares the split edge as a whole — check the endpoints coincide
    for k in keys[1:17:end]
        n = FV.diamond_partner(base, k)
        n === nothing && continue
        a = FV.key_corners(base, k)
        b = FV.key_corners(base, n)
        @test (a[1] ≈ b[1] && a[3] ≈ b[3]) || (a[1] ≈ b[3] && a[3] ≈ b[1])
    end
end

@testset "isubd conforming refinement is watertight" begin
    base = diamond_base_adj()
    mesh = FV.IsubdMesh(base)
    scratch = UInt64[]
    # the step criterion jumps 7 levels across x = 0.7 — the hard case for
    # conformity, and the one the unguarded scheme leaves full of T-vertices
    for (lod, name) in ((StepLoD(7), "step"), (FV.UniformLoD(4), "uniform"))
        keys = FV.root_keys(base)
        FV.refine_keys!(keys, scratch, base, lod; max_depth = 10, conforming = true)
        FV.decode_keys!(mesh, keys, base)
        @test count_tjunctions(mesh) == 0
        @test sum(_area(mesh.positions[f[1]], mesh.positions[f[2]], mesh.positions[f[3]])
                  for f in mesh.faces) ≈ 1.0        # exact tiling, no overlap or hole
        # watertight in the strict sense: every drawn edge is shared by exactly
        # two triangles unless it lies on the domain boundary
        interior, boundary, over = edge_multiplicities(mesh)
        @test interior == 0 && over == 0 && boundary > 0
    end
    # the non-conforming pass on the same criterion does leave T-vertices
    keys = FV.root_keys(base)
    FV.refine_keys!(keys, scratch, base, StepLoD(7); max_depth = 10, conforming = false)
    FV.decode_keys!(mesh, keys, base)
    @test count_tjunctions(mesh) > 0

    # conforming merge collapses whole diamonds and returns to the base
    keys = FV.root_keys(base)
    FV.refine_keys!(keys, scratch, base, StepLoD(7); max_depth = 10, conforming = true)
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(0); max_depth = 10, conforming = true)
    @test sort(keys) == sort(FV.root_keys(base))
    # ... and every intermediate state stays watertight
    keys = FV.root_keys(base)
    FV.refine_keys!(keys, scratch, base, StepLoD(7); max_depth = 10, conforming = true)
    leaves = Set(keys)
    for _ in 1:12
        FV.conforming_update!(leaves, base, FV.UniformLoD(0); max_depth = 10) || break
        FV.decode_keys!(mesh, sort!(collect(leaves)), base)
        @test count_tjunctions(mesh) == 0
        @test sum(_area(mesh.positions[f[1]], mesh.positions[f[2]], mesh.positions[f[3]])
                  for f in mesh.faces) ≈ 1.0
    end
    @test sort(collect(leaves)) == sort(FV.root_keys(base))
end

@testset "isubd decode buffers and mapping" begin
    # a curved (paraboloid) mapping: decode must evaluate it per vertex
    bump = (b, ξ) -> Vec((ξ[1], ξ[2], ξ[1] * (1 - ξ[1]) + ξ[2] * (1 - ξ[2])))
    base = diamond_base(bump)
    keys, scratch = FV.root_keys(base), UInt64[]
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(2))
    mesh = FV.IsubdMesh(base)
    FV.decode_keys!(mesh, keys, base)
    # vertices are shared within each base triangle's subtree
    @test length(mesh.positions) == length(mesh.refcoords) == length(mesh.vertex_base)
    @test length(mesh.positions) < 3 * length(keys)
    @test length(mesh.faces) == length(keys)
    @test all(mesh.positions[i] ≈ bump(0, mesh.refcoords[i]) for i in eachindex(mesh.positions))
    # faces have consistent winding in reference space despite the depth-parity flip
    refarea(f) = _signed_area(mesh.refcoords[f[1]], mesh.refcoords[f[2]], mesh.refcoords[f[3]])
    @test all(sign(refarea(f)) == sign(refarea(mesh.faces[1])) for f in mesh.faces)
    # steady re-decode reuses the buffers
    p, r, f = mesh.positions, mesh.refcoords, mesh.faces
    FV.decode_keys!(mesh, keys, base)
    @test mesh.positions === p && mesh.refcoords === r && mesh.faces === f
    # zero-allocation steady state, measured behind a function barrier: on
    # Julia 1.10 an `@allocated` directly in the testset body measures the
    # body's own dynamic-dispatch overhead (48 bytes) on top of the call
    measure_decode(mesh, keys, base) = @allocated FV.decode_keys!(mesh, keys, base)
    measure_decode(mesh, keys, base)
    @test measure_decode(mesh, keys, base) == 0
end
