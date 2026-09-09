# Layer 1b: the quadrature-point Voronoi tessellation.
#
# Internal variables carry a value only at the quadrature points — there is no
# interpolation attached to them that would say what the value is anywhere else
# (in FEM terms: they live in an L2 space, with no continuity to exploit). The
# honest way to draw them is therefore to partition each cell into the Voronoi
# regions of its quadrature points and paint every region with a single (flat)
# value, rather than to invent a nodal field by averaging or smoothing.
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
vertices in reference space, `triangles` (the rendered boundary surface) and
`simplices` (the refdim-dimensional volume decomposition of the regions) index
into them, and `vertex_qp[v]` is the quadrature point whose region vertex `v`
belongs to.

Vertices are *not* shared between regions: every region carries its own copy, so
assigning each vertex its quadrature point's value renders the region flat — and
no simplex ever spans two regions, which is what keeps an exact [`Clip`](@ref)
of quadrature-point data region-faithful.
"""
struct QPTessellation{refdim,S}
    coords::Vector{Ferrite.Vec{refdim,Float64}}
    triangles::Vector{NTuple{3,Int}}
    simplices::Vector{S}
    vertex_qp::Vector{Int}
end

nvertices(tess::QPTessellation) = length(tess.coords)
ntriangles(tess::QPTessellation) = length(tess.triangles)
nsimplices(tess::QPTessellation) = length(tess.simplices)

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

# The Voronoi regions of a rule on the faces of a reference shape, as
# `(face_index, qp_index, polygon)` triples — the polygon's vertices are in
# ring order. The shared core of `qp_voronoi_tessellation` (which fills the
# regions with flat triangles for the static renderer) and the adaptive base
# construction (which fans each region so its boundaries become split edges).
function _qp_face_regions(::Type{RS}, qr::Ferrite.QuadratureRule) where {RS<:Ferrite.AbstractRefShape}
    corners = Ferrite.reference_coordinates(Ferrite.Lagrange{RS,1}())
    V = eltype(corners)
    refdim = length(first(corners))
    ξs = Ferrite.getpoints(qr)
    length(first(ξs)) == refdim ||
        error("quadrature rule has reference dimension $(length(first(ξs))), expected $refdim for $RS")
    regions = Tuple{Int,Int,Vector{V}}[]
    refdim < 2 && return regions
    tol = 1e-12
    for (fi, face) in enumerate(Ferrite.reference_faces(RS))
        base = V[corners[k] for k in face]
        for i in 1:length(ξs)
            poly = base
            for j in 1:length(ξs)
                j == i && continue
                # keep the side closer to ξᵢ: (ξⱼ-ξᵢ)⋅x ≤ (|ξⱼ|²-|ξᵢ|²)/2
                poly = _clip_halfspace(poly, ξs[j] - ξs[i],
                                       (sum(abs2, ξs[j]) - sum(abs2, ξs[i])) / 2, tol)
                length(poly) < 3 && break
            end
            poly = _dedup_polygon(poly, tol)
            length(poly) < 3 && continue   # this region does not meet this face
            push!(regions, (fi, i, poly))
        end
    end
    return regions
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
    coords = V[]
    triangles = NTuple{3,Int}[]
    vertex_qp = Int[]
    for (_, i, poly) in _qp_face_regions(RS, qr)
        offset = length(coords)
        append!(coords, poly)
        append!(vertex_qp, fill(i, length(poly)))
        for t in 2:(length(poly)-1)    # fan-triangulate the convex region
            push!(triangles, (offset + 1, offset + t, offset + t + 1))
        end
    end
    simplices = NTuple{refdim + 1,Int}[]
    # in 2D the face fans already tile the region areas — they are the volume
    # decomposition; in 3D the regions are tetrahedralized separately below
    refdim == 2 && append!(simplices, triangles)
    refdim == 3 && _qp_voronoi_volume!(coords, simplices, vertex_qp, RS, corners,
                                       Ferrite.getpoints(qr), 1e-12)
    return QPTessellation{refdim,NTuple{refdim + 1,Int}}(coords, triangles, simplices, vertex_qp)
end

# Simple normal of a planar convex polygon (the reference faces are planar).
function _poly_normal(poly::Vector{V}) where {V<:Ferrite.Vec{3}}
    return Tensors.cross(poly[2] - poly[1], poly[3] - poly[1])
end

# A rectangle lying on the plane n ⋅ x = c, large enough to cover the cell;
# clipping it by the cell faces and the other bisectors yields the Voronoi wall.
function _plane_rectangle(n::V, c::Float64, corners::Vector{V}) where {V<:Ferrite.Vec{3}}
    n̂ = n / LinearAlgebra.norm(n)
    centroid = sum(corners) / length(corners)
    p0 = centroid + (c - n ⋅ centroid) / (n ⋅ n) * n
    e = abs(n̂[1]) < 0.9 ? Ferrite.Vec(1.0, 0.0, 0.0) : Ferrite.Vec(0.0, 1.0, 0.0)
    u = Tensors.cross(n̂, e)
    u /= LinearAlgebra.norm(u)
    v = Tensors.cross(n̂, u)
    R = 4.0 * maximum(LinearAlgebra.norm(x - centroid) for x in corners)
    return V[p0 + R * u + R * v, p0 - R * u + R * v, p0 - R * u - R * v, p0 + R * u - R * v]
end

# Volumetric 3D Voronoi partition: every region is decomposed into tets by
# fanning its boundary polygons — the cell-face portions (which are also the
# rendered triangles) and the bisector wall polygons — from its quadrature
# point, which lies in (the closure of) its own convex region, so the cone fan
# tiles it. Every region carries its own vertex copies, so no tet ever spans a
# region wall.
function _qp_voronoi_volume!(coords::Vector{V}, simplices, vertex_qp, ::Type{RS}, corners::Vector{V}, ξs, tol) where {V,RS}
    nqp = length(ξs)
    # Degenerate rules render fine on the surface partition but have no valid
    # volume decomposition; skip it (with a hint) instead of failing the whole
    # apply — Clip/ExtractIsosurfaces then treat these cells as volume-less.
    for i in 1:nqp, j in (i+1):nqp
        LinearAlgebra.norm(ξs[i] - ξs[j]) > 1e-9 && continue
        @warn "quadrature points $i and $j coincide; the volumetric Voronoi partition needs pairwise " *
              "distinct points, so cells with this rule get no volume decomposition (Clip cuts their " *
              "surface without caps, ExtractIsosurfaces skips them)" maxlog = 1
        return nothing
    end
    base_faces = [V[corners[k] for k in face] for face in Ferrite.reference_faces(RS)]
    centroid = sum(corners) / length(corners)
    face_hs = map(base_faces) do poly
        n = _poly_normal(poly)
        n /= LinearAlgebra.norm(n)
        c = n ⋅ poly[1]
        n ⋅ centroid > c ? (-n, -c) : (n, c)   # oriented outward: inside is n ⋅ x ≤ c
    end
    # The point must lie in (the closure of) the cell so its convex Voronoi
    # region can be fanned from it. Boundary points are fine — standard rules
    # have them (e.g. the order-2 prism rule) — the cone over a face incident
    # to the point is just degenerate and pruned below.
    for (q, ξ) in enumerate(ξs), (n, c) in face_hs
        ξ ⋅ n <= c + 1e-9 && continue
        @warn "quadrature point $q at $ξ lies outside the reference cell; the volumetric Voronoi " *
              "partition fans each region from its point, so cells with this rule get no volume " *
              "decomposition (Clip cuts their surface without caps, ExtractIsosurfaces skips them)" maxlog = 1
        return nothing
    end
    # keep the side closer to ξᵢ: (ξⱼ-ξᵢ)⋅x ≤ (|ξⱼ|²-|ξᵢ|²)/2
    bisector(i, j) = (ξs[j] - ξs[i], (sum(abs2, ξs[j]) - sum(abs2, ξs[i])) / 2)
    for i in 1:nqp
        polys = Vector{V}[]
        for base in base_faces
            poly = base
            for j in 1:nqp
                j == i && continue
                poly = _clip_halfspace(poly, bisector(i, j)..., tol)
                length(poly) < 3 && break
            end
            poly = _dedup_polygon(poly, tol)
            length(poly) >= 3 && push!(polys, poly)
        end
        for j in 1:nqp
            j == i && continue
            wall = _plane_rectangle(bisector(i, j)..., corners)
            for (nf, cf) in face_hs
                wall = _clip_halfspace(wall, nf, cf, tol)
                length(wall) < 3 && break
            end
            for k in 1:nqp
                (k == i || k == j || length(wall) < 3) && continue
                wall = _clip_halfspace(wall, bisector(i, k)..., tol)
            end
            wall = _dedup_polygon(wall, tol)
            length(wall) >= 3 && push!(polys, wall)
        end
        center = length(coords) + 1
        push!(coords, ξs[i])
        push!(vertex_qp, i)
        for poly in polys
            offset = length(coords)
            append!(coords, poly)
            append!(vertex_qp, fill(i, length(poly)))
            for t in 2:(length(poly)-1)
                # cones over faces incident to a boundary quadrature point are
                # flat — skip them; the rest are stored positively oriented
                vol = ((poly[t] - poly[1]) × (poly[t+1] - poly[1])) ⋅ (ξs[i] - poly[1]) / 6
                abs(vol) <= 1e-12 && continue
                push!(simplices, vol > 0 ? (offset + 1, offset + t, offset + t + 1, center) :
                                           (offset + t, offset + 1, offset + t + 1, center))
            end
        end
    end
    return nothing
end
