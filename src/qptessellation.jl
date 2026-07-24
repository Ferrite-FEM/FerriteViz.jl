# Layer 1b: the quadrature-point Voronoi tessellation.
#
# Internal variables are L2 functions known only at the quadrature points, so the
# honest way to draw them is to partition each cell into the Voronoi regions of
# its quadrature points and paint every region with a single (flat) value.
#
# The partition is built once per reference shape, in reference space, and then
# instantiated per cell through the geometric map — exactly like
# `reference_tessellation`. Because `Ferrite.reference_faces` is uniform across
# dimensions (a 2D shape reports itself as its single face, a 3D shape its
# boundary faces) one algorithm covers every shape: clip each face polygon by the
# bisector half-spaces of the quadrature points.

"""
    QPTessellation{refdim}

Voronoi partition of a reference shape induced by a quadrature rule. `coords` are
vertices in reference space, `triangles` index into them, and `vertex_qp[v]` is
the quadrature point whose region vertex `v` belongs to.

Vertices are *not* shared between regions: every region carries its own copy, so
assigning each vertex its quadrature point's value renders the region flat.
"""
struct QPTessellation{refdim}
    coords::Vector{Ferrite.Vec{refdim,Float64}}
    triangles::Vector{NTuple{3,Int}}
    vertex_qp::Vector{Int}
end

nvertices(tess::QPTessellation) = length(tess.coords)
ntriangles(tess::QPTessellation) = length(tess.triangles)

# Sutherland–Hodgman clip of a convex polygon against the half-space `n ⋅ x ≤ c`.
# The polygon is planar but may be embedded in 2D or 3D — the bisector of two
# points is a plane, so clipping a planar polygon by it stays planar.
function _clip_halfspace(poly::Vector{V}, n::V, c::Float64, tol::Float64) where {V<:Ferrite.Vec}
    m = length(poly)
    m < 3 && return V[]
    out = V[]
    for i in 1:m
        a = poly[i]
        b = poly[mod1(i + 1, m)]
        da = (n ⋅ a) - c
        db = (n ⋅ b) - c
        da <= tol && push!(out, a)
        # only add an intersection for a genuine sign change (points on the
        # plane are already kept above)
        if (da > tol && db < -tol) || (da < -tol && db > tol)
            push!(out, a + (da / (da - db)) * (b - a))
        end
    end
    return out
end

# Drop consecutive (and wrap-around) duplicates so clipping cannot emit
# zero-area slivers.
function _dedup_polygon(poly::Vector{V}, tol::Float64) where {V<:Ferrite.Vec}
    isempty(poly) && return poly
    out = V[]
    for p in poly
        (isempty(out) || LinearAlgebra.norm(p - out[end]) > tol) && push!(out, p)
    end
    while length(out) > 1 && LinearAlgebra.norm(out[end] - out[1]) <= tol
        pop!(out)
    end
    return out
end

"""
    qp_voronoi_tessellation(::Type{<:Ferrite.AbstractRefShape}, qr::Ferrite.QuadratureRule) -> QPTessellation

Partition a reference shape into the Voronoi regions of `qr`'s quadrature points
and triangulate them. For 3D shapes the partition is intersected with the
boundary faces, which is what the surface renderer draws.
"""
function qp_voronoi_tessellation(::Type{RS}, qr::Ferrite.QuadratureRule) where {RS<:Ferrite.AbstractRefShape}
    corners = Ferrite.reference_coordinates(Ferrite.Lagrange{RS,1}())
    V = eltype(corners)
    refdim = length(first(corners))
    ξs = Ferrite.getpoints(qr)
    length(first(ξs)) == refdim ||
        error("quadrature rule has reference dimension $(length(first(ξs))), expected $refdim for $RS")

    coords = V[]
    triangles = NTuple{3,Int}[]
    vertex_qp = Int[]
    # a line has no surface to triangulate (mirrors reference_tessellation(RefLine))
    refdim < 2 && return QPTessellation{refdim}(coords, triangles, vertex_qp)

    nqp = length(ξs)
    tol = 1e-12
    for face in Ferrite.reference_faces(RS)
        base = V[corners[k] for k in face]
        for i in 1:nqp
            poly = base
            for j in 1:nqp
                j == i && continue
                # keep the side closer to ξᵢ: (ξⱼ-ξᵢ)⋅x ≤ (|ξⱼ|²-|ξᵢ|²)/2
                poly = _clip_halfspace(poly, ξs[j] - ξs[i],
                                       (sum(abs2, ξs[j]) - sum(abs2, ξs[i])) / 2, tol)
                length(poly) < 3 && break
            end
            poly = _dedup_polygon(poly, tol)
            length(poly) < 3 && continue   # this region does not meet this face
            offset = length(coords)
            append!(coords, poly)
            append!(vertex_qp, fill(i, length(poly)))
            for t in 2:(length(poly)-1)    # fan-triangulate the convex region
                push!(triangles, (offset + 1, offset + t, offset + t + 1))
            end
        end
    end
    return QPTessellation{refdim}(coords, triangles, vertex_qp)
end
