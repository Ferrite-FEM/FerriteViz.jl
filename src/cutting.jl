# Layer 1c: cutting infrastructure shared by Clip and ExtractIsosurfaces.
#
# Everything here is expressed in terms of *affine combinations of parent
# vertices*: a derived dataset's vertex is Σ wₖ·parent(pₖ). Storing the sparse
# (parents, weights) table once gives every consumer the same uniform rule —
# coordinates stay `lift`ed from the parent's (so update!-driven deformation
# keeps propagating), registered point data is linearly interpolated with the
# identical weights, and reference coordinates are combined statically (so dof
# fields re-evaluate exactly at the derived positions). All arithmetic is done
# in Float64 even though coordinates are stored as Float32.

"""
Sparse affine combinations mapping parent vertices to derived vertices:
derived vertex `v` is `Σ weights[k]·parent[parents[k]]` for
`k in offsets[v]+1:offsets[v+1]`. Built by [`combo_identity!`](@ref) and
[`combo_pair!`](@ref); applied with
[`combine_points`](@ref)/[`combine_rows`](@ref).
"""
struct AffineCombinations
    offsets::Vector{Int}
    parents::Vector{Int}
    weights::Vector{Float64}
end
AffineCombinations() = AffineCombinations([0], Int[], Float64[])

ncombos(ac::AffineCombinations) = length(ac.offsets) - 1

function _push_combo!(ac::AffineCombinations, parents, weights)
    for (p, w) in zip(parents, weights)
        w == 0 && continue
        push!(ac.parents, p)
        push!(ac.weights, w)
    end
    push!(ac.offsets, length(ac.parents))
    return ncombos(ac)
end

"Append the identity combination (the derived vertex is parent `p`)."
combo_identity!(ac::AffineCombinations, p::Int) = _push_combo!(ac, (p,), (1.0,))

"Append the point at parameter `t` on the parent edge `a`→`b`."
combo_pair!(ac::AffineCombinations, a::Int, b::Int, t::Float64) = _push_combo!(ac, (a, b), (1.0 - t, t))

"""
    combine_points(ac, parent_coords) -> Vector{Point}

Apply the combinations to a vector of `GeometryBasics.Point`s (accumulating in
Float64). This is the function a derived dataset's coordinate `lift` applies to
its parent's coordinates.
"""
function combine_points(ac::AffineCombinations, parent_coords::AbstractVector{GeometryBasics.Point{dim,T}}) where {dim,T}
    out = Vector{GeometryBasics.Point{dim,T}}(undef, ncombos(ac))
    for v in 1:ncombos(ac)
        # anchored form  x₁ + Σₖ wₖ(xₖ - x₁)  (valid since Σ wₖ = 1): exact when
        # all parents coincide, which keeps piecewise-constant data (Voronoi
        # regions) exactly flat through a cut
        first_k = ac.offsets[v] + 1
        anchor = GeometryBasics.Point{dim,Float64}(parent_coords[ac.parents[first_k]])
        acc = anchor
        for k in (first_k+1):(ac.offsets[v+1])
            acc = acc .+ ac.weights[k] .* (GeometryBasics.Point{dim,Float64}(parent_coords[ac.parents[k]]) .- anchor)
        end
        out[v] = GeometryBasics.Point{dim,T}(acc)
    end
    return out
end

"""
    combine_rows(ac, A::AbstractMatrix) -> Matrix{Float64}

Apply the combinations rowwise to a per-vertex data matrix (registered point
data, reference coordinates).
"""
function combine_rows(ac::AffineCombinations, A::AbstractMatrix)
    out = zeros(Float64, ncombos(ac), size(A, 2))
    for v in 1:ncombos(ac)
        # anchored form, see combine_points
        first_k = ac.offsets[v] + 1
        a = ac.parents[first_k]
        for d in 1:size(A, 2)
            out[v, d] = A[a, d]
        end
        for k in (first_k+1):(ac.offsets[v+1])
            w, p = ac.weights[k], ac.parents[k]
            for d in 1:size(A, 2)
                out[v, d] += w * (A[p, d] - A[a, d])
            end
        end
    end
    return out
end

###################
# Cut vertex pool #
###################

# Builds the output vertex set of a cut: kept parent vertices enter through
# `out_vertex!` (memoized identity combinations), vertices on cut edges through
# `cut_vertex!` (memoized on the unordered parent pair plus a tag, so every
# primitive of a cell sharing an edge shares the cut vertex — the tag separates
# cuts of the same edge at different isosurface levels). Because edges never
# span cells, processing cell by cell keeps the output vertices of a cell
# contiguous. Positions are tracked in Float64 for the degeneracy checks;
# parent coordinates are converted on access rather than materialized (a large
# dataset's coordinate array is sizable, and a cut touches only part of it).
struct CutVertexPool{dim,PC<:AbstractVector}
    parent_points::PC
    combos::AffineCombinations
    pos::Vector{Tensors.Vec{dim,Float64}}
    remap::Vector{Int}               # parent vertex -> output vertex; 0 = not emitted
    cuts::Dict{NTuple{3,Int},Int}
    # scratch buffers of the per-primitive kernels below (the pool is used
    # strictly sequentially), so the hot loops allocate nothing
    poly_buf::Vector{Int}
    ins_buf::Vector{Int}
    outs_buf::Vector{Int}
end

function CutVertexPool(parent_coords::PC) where {dim,T,PC<:AbstractVector{GeometryBasics.Point{dim,T}}}
    return CutVertexPool{dim,PC}(parent_coords, AffineCombinations(), Tensors.Vec{dim,Float64}[],
                                 zeros(Int, length(parent_coords)), Dict{NTuple{3,Int},Int}(),
                                 Int[], Int[], Int[])
end

nvertices(pool::CutVertexPool) = ncombos(pool.combos)

@inline _parent_pos(pool::CutVertexPool{dim}, i::Int) where {dim} =
    Tensors.Vec{dim,Float64}(NTuple{dim,Float64}(pool.parent_points[i]))

# Reserve capacity for `n` output vertices (a clip's output is close to its
# input size; growing the columnar arrays incrementally dominates otherwise).
function Base.sizehint!(pool::CutVertexPool, n::Integer)
    sizehint!(pool.pos, n)
    sizehint!(pool.combos.offsets, n + 1)
    sizehint!(pool.combos.parents, n)
    sizehint!(pool.combos.weights, n)
    return pool
end

function out_vertex!(pool::CutVertexPool, i::Int)
    r = pool.remap[i]
    r != 0 && return r
    push!(pool.pos, _parent_pos(pool, i))
    r = combo_identity!(pool.combos, i)
    pool.remap[i] = r
    return r
end

function cut_vertex!(pool::CutVertexPool, i::Int, j::Int, gi::Float64, gj::Float64, tag::Int)
    a, b = i < j ? (i, j) : (j, i)
    return get!(pool.cuts, (a, b, tag)) do
        ga, gb = a == i ? (gi, gj) : (gj, gi)
        t = ga / (ga - gb)   # deterministic: always parametrized from the lower index
        push!(pool.pos, (1.0 - t) * _parent_pos(pool, a) + t * _parent_pos(pool, b))
        combo_pair!(pool.combos, a, b, t)
    end
end

######################
# Geometric measures #
######################

_signed_tet_volume(a, b, c, d) = (((b - a) × (c - a)) ⋅ (d - a)) / 6
_tri_area(a::Tensors.Vec{3}, b, c) = LinearAlgebra.norm((b - a) × (c - a)) / 2
_tri_area(a::Tensors.Vec{2}, b, c) = abs((b - a)[1] * (c - a)[2] - (b - a)[2] * (c - a)[1]) / 2

_tri_area_at(pool, t) = _tri_area(pool.pos[t[1]], pool.pos[t[2]], pool.pos[t[3]])

# Snap scalar values within `tol` of the cut to exactly 0. Classification below
# is a pure function of the snapped values ("zero counts as inside"), so it is
# identical for coincident duplicated vertices across neighboring cells.
function snap!(g::Vector{Float64}, tol::Float64)
    for i in eachindex(g)
        abs(g[i]) <= tol && (g[i] = 0.0)
    end
    return g
end

_inside(g) = g <= 0

####################
# Half-space clips #
####################

# Sutherland–Hodgman clip of one triangle against g ≤ 0; appends the kept part
# (0–2 triangles, fan-triangulated) with degenerate slivers dropped.
function clip_triangle!(pool::CutVertexPool, out_tris::Vector{NTuple{3,Int}}, out_cells::Vector{Int},
                        tri::NTuple{3,Int}, cell::Int, g::Vector{Float64}, area_tol::Float64)
    poly = empty!(pool.poly_buf)
    for k in 1:3
        i, j = tri[k], tri[mod1(k + 1, 3)]
        gi, gj = g[i], g[j]
        _inside(gi) && push!(poly, out_vertex!(pool, i))
        if _inside(gi) ⊻ _inside(gj)
            push!(poly, cut_vertex!(pool, i, j, gi, gj, 0))
        end
    end
    for t in 2:(length(poly) - 1)
        a, b, c = poly[1], poly[t], poly[t+1]
        (a == b || b == c || a == c) && continue
        _tri_area(pool.pos[a], pool.pos[b], pool.pos[c]) <= area_tol && continue
        push!(out_tris, (a, b, c))
        push!(out_cells, cell)
    end
    return nothing
end

# Clip one wireframe edge segment against g ≤ 0 (0–1 segments appended).
function clip_edge!(pool::CutVertexPool, out_edges::Vector{NTuple{2,Int}}, out_cells::Vector{Int},
                    edge::NTuple{2,Int}, cell::Int, g::Vector{Float64})
    i, j = edge
    gi, gj = g[i], g[j]
    if _inside(gi) && _inside(gj)
        p, q = out_vertex!(pool, i), out_vertex!(pool, j)
    elseif _inside(gi)
        p, q = out_vertex!(pool, i), cut_vertex!(pool, i, j, gi, gj, 0)
    elseif _inside(gj)
        p, q = cut_vertex!(pool, i, j, gi, gj, 0), out_vertex!(pool, j)
    else
        return nothing
    end
    # an endpoint exactly on the plane yields a coincident cut vertex (the
    # duplication contract) — such a zero-length remainder is dropped
    (p == q || pool.pos[p] == pool.pos[q]) && return nothing
    push!(out_edges, (p, q))
    push!(out_cells, cell)
    return nothing
end

function _push_tet!(pool, out_tets, out_cells, tet, cell, vol_tol)
    (tet[1] == tet[2] || tet[1] == tet[3] || tet[1] == tet[4] ||
     tet[2] == tet[3] || tet[2] == tet[4] || tet[3] == tet[4]) && return nothing
    vol = _signed_tet_volume(pool.pos[tet[1]], pool.pos[tet[2]], pool.pos[tet[3]], pool.pos[tet[4]])
    abs(vol) <= vol_tol && return nothing
    # normalize to positive orientation (the simplices invariant)
    push!(out_tets, vol > 0 ? tet : (tet[2], tet[1], tet[3], tet[4]))
    push!(out_cells, cell)
    return nothing
end

function _push_tri!(pool, out_tris, out_cells, tri, cell, area_tol)
    (tri[1] == tri[2] || tri[2] == tri[3] || tri[1] == tri[3]) && return nothing
    _tri_area_at(pool, tri) <= area_tol && return nothing
    push!(out_tris, tri)
    push!(out_cells, cell)
    return nothing
end

# Split a convex planar quad (cyclic order) into two triangles along the
# canonical diagonal anchored at the smallest output index, so a face shared by
# two primitives splits identically.
function _push_quad!(pool, out_tris, out_cells, quad::NTuple{4,Int}, cell, area_tol)
    r = argmin(quad) - 1
    q = ntuple(k -> quad[mod1(k + r, 4)], 4)
    _push_tri!(pool, out_tris, out_cells, (q[1], q[2], q[3]), cell, area_tol)
    _push_tri!(pool, out_tris, out_cells, (q[1], q[3], q[4]), cell, area_tol)
    return nothing
end

"""
Clip one tetrahedron against the half-space `g ≤ 0`: the kept region is
re-tetrahedralized into `out_tets` and the cut face goes to `out_caps` (both
tagged with `cell`). The case table is exact for planar cuts of straight-edged
tets, which is what "exact clipping" means for linear (and linearized) cells.
"""
function clip_tet!(pool::CutVertexPool, out_tets::Vector{NTuple{4,Int}}, out_tet_cells::Vector{Int},
                   out_caps::Vector{NTuple{3,Int}}, out_cap_cells::Vector{Int},
                   tet::NTuple{4,Int}, cell::Int, g::Vector{Float64}, vol_tol::Float64, area_tol::Float64)
    nin = count(v -> _inside(g[v]), tet)
    nin == 0 && return nothing
    if nin == 4
        t = (out_vertex!(pool, tet[1]), out_vertex!(pool, tet[2]), out_vertex!(pool, tet[3]), out_vertex!(pool, tet[4]))
        push!(out_tets, t)
        push!(out_tet_cells, cell)
        return nothing
    end
    ins = empty!(pool.ins_buf)
    outs = empty!(pool.outs_buf)
    for v in tet
        _inside(g[v]) ? push!(ins, v) : push!(outs, v)
    end
    cut(i, j) = cut_vertex!(pool, i, j, g[i], g[j], 0)
    if nin == 1
        a = ins[1]
        pab, pac, pad = cut(a, outs[1]), cut(a, outs[2]), cut(a, outs[3])
        _push_tet!(pool, out_tets, out_tet_cells, (out_vertex!(pool, a), pab, pac, pad), cell, vol_tol)
        _push_tri!(pool, out_caps, out_cap_cells, (pab, pac, pad), cell, area_tol)
    elseif nin == 3
        d = outs[1]
        A, B, C = out_vertex!(pool, ins[1]), out_vertex!(pool, ins[2]), out_vertex!(pool, ins[3])
        pa, pb, pc = cut(ins[1], d), cut(ins[2], d), cut(ins[3], d)
        _push_tet!(pool, out_tets, out_tet_cells, (A, B, C, pa), cell, vol_tol)
        _push_tet!(pool, out_tets, out_tet_cells, (B, C, pa, pb), cell, vol_tol)
        _push_tet!(pool, out_tets, out_tet_cells, (C, pa, pb, pc), cell, vol_tol)
        _push_tri!(pool, out_caps, out_cap_cells, (pa, pb, pc), cell, area_tol)
    else # nin == 2
        a, b = ins[1], ins[2]
        c, d = outs[1], outs[2]
        A, B = out_vertex!(pool, a), out_vertex!(pool, b)
        pac, pad, pbc, pbd = cut(a, c), cut(a, d), cut(b, c), cut(b, d)
        # prism (A, pac, pad | B, pbc, pbd): all quad faces are planar (two lie
        # on original tet faces, one is the cap), so the 3-tet pattern is exact.
        # The pattern's cap faces use one diagonal of the cap quad — pick the
        # prism ordering whose diagonal matches the canonical `_push_quad!`
        # split, so the rendered cap is exactly the kept tets' boundary.
        quad = (pac, pad, pbd, pbc)
        m = quad[argmin(quad)]
        if m == pac || m == pbd  # canonical diagonal pac–pbd
            _push_tet!(pool, out_tets, out_tet_cells, (A, pad, pac, B), cell, vol_tol)
            _push_tet!(pool, out_tets, out_tet_cells, (pad, pac, B, pbd), cell, vol_tol)
            _push_tet!(pool, out_tets, out_tet_cells, (pac, B, pbd, pbc), cell, vol_tol)
        else                     # canonical diagonal pad–pbc
            _push_tet!(pool, out_tets, out_tet_cells, (A, pac, pad, B), cell, vol_tol)
            _push_tet!(pool, out_tets, out_tet_cells, (pac, pad, B, pbc), cell, vol_tol)
            _push_tet!(pool, out_tets, out_tet_cells, (pad, B, pbc, pbd), cell, vol_tol)
        end
        _push_quad!(pool, out_caps, out_cap_cells, quad, cell, area_tol)
    end
    return nothing
end

######################
# Marching simplices #
######################

"""
March one tetrahedron for the level set `g = 0` (0–2 triangles into `out_tris`).
The side rule "zero counts as inside" makes a face lying exactly on the level
set be emitted by exactly one of its two adjacent tets (the one whose remaining
vertex is strictly outside).
"""
function march_tet!(pool::CutVertexPool, out_tris::Vector{NTuple{3,Int}}, out_cells::Vector{Int},
                    tet::NTuple{4,Int}, cell::Int, g::Vector{Float64}, area_tol::Float64, tag::Int)
    nin = count(v -> _inside(g[v]), tet)
    (nin == 0 || nin == 4) && return nothing
    ins = empty!(pool.ins_buf)
    outs = empty!(pool.outs_buf)
    for v in tet
        _inside(g[v]) ? push!(ins, v) : push!(outs, v)
    end
    cut(i, j) = cut_vertex!(pool, i, j, g[i], g[j], tag)
    if nin == 1
        a = ins[1]
        _push_tri!(pool, out_tris, out_cells, (cut(a, outs[1]), cut(a, outs[2]), cut(a, outs[3])), cell, area_tol)
    elseif nin == 3
        d = outs[1]
        _push_tri!(pool, out_tris, out_cells, (cut(ins[1], d), cut(ins[2], d), cut(ins[3], d)), cell, area_tol)
    else
        a, b = ins[1], ins[2]
        c, d = outs[1], outs[2]
        _push_quad!(pool, out_tris, out_cells, (cut(a, c), cut(a, d), cut(b, d), cut(b, c)), cell, area_tol)
    end
    return nothing
end

"""
March one triangle for the level set `g = 0` (0–1 segments into `out_segs`).
"""
function march_triangle!(pool::CutVertexPool, out_segs::Vector{NTuple{2,Int}}, out_cells::Vector{Int},
                         tri::NTuple{3,Int}, cell::Int, g::Vector{Float64}, len_tol::Float64, tag::Int)
    nin = count(v -> _inside(g[v]), tri)
    (nin == 0 || nin == 3) && return nothing
    ins = empty!(pool.ins_buf)
    outs = empty!(pool.outs_buf)
    for v in tri
        _inside(g[v]) ? push!(ins, v) : push!(outs, v)
    end
    cut(i, j) = cut_vertex!(pool, i, j, g[i], g[j], tag)
    p, q = nin == 1 ? (cut(ins[1], outs[1]), cut(ins[1], outs[2])) :
                      (cut(ins[1], outs[1]), cut(ins[2], outs[1]))
    p == q && return nothing
    LinearAlgebra.norm(pool.pos[q] - pool.pos[p]) <= len_tol && return nothing
    push!(out_segs, (p, q))
    push!(out_cells, cell)
    return nothing
end
