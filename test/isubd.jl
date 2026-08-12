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
function count_tjunctions(mesh; digits = 9)
    key(p) = (round(p[1]; digits = digits), round(p[2]; digits = digits))
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

# a hand-rolled perspective camera looking at the unit square from `eye`
function perspective_projview(eye::NTuple{3,Float64}; target = (0.5, 0.5, 0.0))
    # look-at basis
    f = [target[1] - eye[1], target[2] - eye[2], target[3] - eye[3]]
    f ./= sqrt(sum(abs2, f))
    upv = abs(f[3]) > 0.99 ? [0.0, 1.0, 0.0] : [0.0, 0.0, 1.0]
    s = [f[2] * upv[3] - f[3] * upv[2], f[3] * upv[1] - f[1] * upv[3], f[1] * upv[2] - f[2] * upv[1]]
    s ./= sqrt(sum(abs2, s))
    u = [s[2] * f[3] - s[3] * f[2], s[3] * f[1] - s[1] * f[3], s[1] * f[2] - s[2] * f[1]]
    view = [s[1] s[2] s[3] -sum(s .* collect(eye));
            u[1] u[2] u[3] -sum(u .* collect(eye));
            -f[1] -f[2] -f[3] sum(f .* collect(eye));
            0.0 0.0 0.0 1.0]
    fovfac = 1.0 / tan(π / 8)
    near, far = 0.01, 100.0
    proj = [fovfac 0.0 0.0 0.0;
            0.0 fovfac 0.0 0.0;
            0.0 0.0 -(far + near) / (far - near) -2far * near / (far - near);
            0.0 0.0 -1.0 0.0]
    return proj * view
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

@testset "isubd screen-space LoD" begin
    base = diamond_base((b, ξ) -> Vec((ξ[1], ξ[2], 0.0)))   # square in the z=0 plane
    resolution = (800.0, 600.0)
    near = FV.ScreenSpaceLoD(perspective_projview((0.05, 0.05, 0.4)), (0.05, 0.05, 0.4), resolution, 30.0)
    far = FV.ScreenSpaceLoD(perspective_projview((0.5, 0.5, 20.0)), (0.5, 0.5, 20.0), resolution, 30.0)

    keys, scratch = FV.root_keys(base), UInt64[]
    FV.refine_keys!(keys, scratch, base, near; max_depth = 16)
    nnear = length(keys)
    @test nnear > 2
    # view-dependent: refinement concentrates near the eye
    mesh = FV.IsubdMesh(base)
    FV.decode_keys!(mesh, keys, base)
    depth_near_eye = maximum((FV.key_depth(k) for k in keys if
                              all(c -> norm(c - Vec((0.05, 0.05))) < 0.4, FV.key_corners(base, k)));
                             init=-1)
    depth_far_corner = maximum((FV.key_depth(k) for k in keys if
                                all(c -> norm(c - Vec((1.0, 1.0))) < 0.3, FV.key_corners(base, k)));
                               init=-1)
    @test depth_near_eye >= 0 && depth_far_corner >= 0
    @test depth_near_eye > depth_far_corner
    # crack-free at the converged state
    @test count_tjunctions(mesh) == 0
    # fixed point: another pass with the same camera changes nothing
    before = sort(copy(keys))
    FV.update_keys!(scratch, keys, base, near; max_depth = 16)
    @test sort(scratch) == before
    # retreating the camera coarsens the mesh again, still crack-free
    FV.refine_keys!(keys, scratch, base, far; max_depth = 16)
    @test length(keys) < nnear
    FV.decode_keys!(mesh, keys, base)
    @test count_tjunctions(mesh) == 0
end

@testset "isubd decode buffers and mapping" begin
    # a curved (paraboloid) mapping: decode must evaluate it per vertex
    bump = (b, ξ) -> Vec((ξ[1], ξ[2], ξ[1] * (1 - ξ[1]) + ξ[2] * (1 - ξ[2])))
    base = diamond_base(bump)
    keys, scratch = FV.root_keys(base), UInt64[]
    FV.refine_keys!(keys, scratch, base, FV.UniformLoD(2))
    mesh = FV.IsubdMesh(base)
    FV.decode_keys!(mesh, keys, base)
    @test length(mesh.positions) == 3 * length(keys) == length(mesh.refcoords)
    @test all(mesh.positions[i] ≈ bump(0, mesh.refcoords[i]) for i in eachindex(mesh.positions))
    # faces have consistent winding in reference space despite the depth-parity flip
    refarea(f) = _signed_area(mesh.refcoords[f[1]], mesh.refcoords[f[2]], mesh.refcoords[f[3]])
    @test all(sign(refarea(f)) == sign(refarea(mesh.faces[1])) for f in mesh.faces)
    # steady re-decode reuses the buffers
    p, r, f = mesh.positions, mesh.refcoords, mesh.faces
    FV.decode_keys!(mesh, keys, base)
    @test mesh.positions === p && mesh.refcoords === r && mesh.faces === f
    allocs = @allocated FV.decode_keys!(mesh, keys, base)
    @test allocs == 0
end
