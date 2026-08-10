# Layer 1: the tessellation interface.
#
# A cell's rendering geometry is described entirely by the `ReferenceTessellation`
# of its reference shape. Everything else (triangle counts, physical coordinates,
# midpoints, ...) is derived from it, so supporting a new reference shape means
# implementing a single `reference_tessellation` method.

"""
    ReferenceTessellation{refdim,T}

Surface triangulation of a reference shape: `coords` are vertices in reference
space, `triangles` index into `coords`. Coordinates may be shared between
triangles; the per-cell vertex duplication that makes discontinuous (L2) fields
render correctly is applied by [`FEData`](@ref), not here.

`edges` are the wireframe segments drawn by [`meshplot`](@ref): polylines along
the *finite element cell's* edges (from `Ferrite.reference_edges`), also
indexing into `coords`. They are deliberately separate from the triangles —
the triangulation contains interior diagonals (e.g. the quadrilateral's center
fan) that are not cell edges and must not show up in the wireframe. A
tessellation without edges renders no wireframe.
"""
struct ReferenceTessellation{refdim,T}
    coords::Vector{Ferrite.Vec{refdim,T}}
    triangles::Vector{NTuple{3,Int}}
    edges::Vector{NTuple{2,Int}}
end

# backwards compatible: a tessellation without edges draws no wireframe
ReferenceTessellation(coords, triangles) = ReferenceTessellation(coords, triangles, NTuple{2,Int}[])

nvertices(tess::ReferenceTessellation) = length(tess.coords)
ntriangles(tess::ReferenceTessellation) = length(tess.triangles)
nedges(tess::ReferenceTessellation) = length(tess.edges)

"""
    reference_tessellation(::Type{<:Ferrite.AbstractRefShape}) -> ReferenceTessellation

The extension point for custom cell types: return the surface triangulation of
the reference shape. Implementing this one method makes `FEData` and all
representations work for cells with that reference shape. For 3D shapes,
[`facet_based_tessellation`](@ref) builds a valid tessellation from
`Ferrite.reference_faces`.
"""
function reference_tessellation(::Type{RS}) where {RS<:Ferrite.AbstractRefShape}
    error("No tessellation known for reference shape $RS. \
           Define `FerriteViz.reference_tessellation(::Type{$RS})` to visualize cells with this shape \
           (for 3D shapes `FerriteViz.facet_based_tessellation($RS)` is usually all you need).")
end

# Lines have no surface to triangulate; an empty triangle list still lets
# FEData and meshplot handle grids containing line cells.
function reference_tessellation(::Type{Ferrite.RefLine})
    coords = Ferrite.reference_coordinates(Ferrite.Lagrange{Ferrite.RefLine,1}())
    return ReferenceTessellation(coords, NTuple{3,Int}[], [(1, 2)])
end

function reference_tessellation(::Type{Ferrite.RefTriangle})
    coords = Ferrite.reference_coordinates(Ferrite.Lagrange{Ferrite.RefTriangle,1}())
    return ReferenceTessellation(coords, [(1, 2, 3)], collect(Ferrite.reference_edges(Ferrite.RefTriangle)))
end

# The quadrilateral deliberately decomposes into 4 triangles through the center
# vertex: the cheaper 2-triangle split misses a (bilinear) solution mode.
function reference_tessellation(::Type{Ferrite.RefQuadrilateral})
    coords = Ferrite.reference_coordinates(Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}())
    push!(coords, zero(eltype(coords)))
    return ReferenceTessellation(coords, [(1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 1, 5)],
                                 collect(Ferrite.reference_edges(Ferrite.RefQuadrilateral)))
end

"""
    facet_based_tessellation(::Type{<:Ferrite.AbstractRefShape{3}}) -> ReferenceTessellation

Build the surface tessellation of a 3D reference shape from its faces: each
face's 2D tessellation (triangle or quadrilateral, chosen by vertex count) is
mapped into the element via `Ferrite.facet_to_element_transformation`. The
wireframe edges come from `Ferrite.reference_edges`, with their own (duplicated)
endpoint vertices appended after the face vertices. This is the default
building block for `reference_tessellation` of volumetric shapes.
"""
function facet_based_tessellation(::Type{RS}) where {RS<:Ferrite.AbstractRefShape{3}}
    coords = Ferrite.Vec{3,Float64}[]
    triangles = NTuple{3,Int}[]
    for (face_idx, face) in enumerate(Ferrite.reference_faces(RS))
        face_shape = length(face) == 3 ? Ferrite.RefTriangle :
                     length(face) == 4 ? Ferrite.RefQuadrilateral :
                     error("faces with $(length(face)) vertices are not supported")
        face_tess = reference_tessellation(face_shape)
        offset = length(coords)
        for ξ in face_tess.coords
            push!(coords, Ferrite.facet_to_element_transformation(ξ, RS, face_idx))
        end
        for tri in face_tess.triangles
            push!(triangles, tri .+ offset)
        end
    end
    edges = NTuple{2,Int}[]
    corners = Ferrite.reference_coordinates(Ferrite.Lagrange{RS,1}())
    for (v1, v2) in Ferrite.reference_edges(RS)
        push!(coords, corners[v1], corners[v2])
        push!(edges, (length(coords) - 1, length(coords)))
    end
    return ReferenceTessellation(coords, triangles, edges)
end

reference_tessellation(::Type{Ferrite.RefTetrahedron}) = facet_based_tessellation(Ferrite.RefTetrahedron)
reference_tessellation(::Type{Ferrite.RefHexahedron}) = facet_based_tessellation(Ferrite.RefHexahedron)
reference_tessellation(::Type{Ferrite.RefPrism}) = facet_based_tessellation(Ferrite.RefPrism)
reference_tessellation(::Type{Ferrite.RefPyramid}) = facet_based_tessellation(Ferrite.RefPyramid)

"""
Number of triangles a cell tessellates into.
"""
ntriangles(cell::Ferrite.AbstractCell) = ntriangles(reference_tessellation(Ferrite.getrefshape(cell)))

"""
    subdivide(tess::ReferenceTessellation, n::Int) -> ReferenceTessellation

Subdivide a tessellation in reference space, `n` times: every triangle into 4
(orientation preserving) and every edge segment into 2. Midpoints are
deduplicated by coordinate, so subdivided triangles share vertices with their
neighbours and with edge segments running along the same reference line. New
vertices are appended, so indices into the input tessellation stay valid.

This is what resolves curved geometry: the subdivided reference vertices are
mapped through the cell's geometric interpolation when the tessellation is
instantiated by [`FEData`](@ref), so surfaces and wireframe edges of high-order
(or nonlinearly deformed) cells bend instead of being drawn as flat facets and
straight chords.
"""
function subdivide(tess::ReferenceTessellation{refdim,T}, n::Int) where {refdim,T}
    n <= 0 && return tess
    coords = copy(tess.coords)
    index = Dict{Ferrite.Vec{refdim,T},Int}()
    for (i, ξ) in enumerate(coords)
        get!(index, ξ, i)
    end
    function midpoint(a, b)
        ξ = (coords[a] + coords[b]) / 2
        return get!(index, ξ) do
            push!(coords, ξ)
            length(coords)
        end
    end
    triangles = NTuple{3,Int}[]
    for (a, b, c) in tess.triangles
        mab, mbc, mca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
        push!(triangles, (a, mab, mca), (b, mbc, mab), (c, mca, mbc), (mab, mbc, mca))
    end
    edges = NTuple{2,Int}[]
    for (a, b) in tess.edges
        m = midpoint(a, b)
        push!(edges, (a, m), (m, b))
    end
    return subdivide(ReferenceTessellation(coords, triangles, edges), n - 1)
end

# One edge-only subdivision round: edges are 1D and therefore much cheaper to
# refine than the surface, so their resolution can exceed the triangles'.
function _subdivide_edges(tess::ReferenceTessellation)
    coords = copy(tess.coords)
    edges = NTuple{2,Int}[]
    for (a, b) in tess.edges
        push!(coords, (coords[a] + coords[b]) / 2)
        push!(edges, (a, length(coords)), (length(coords), b))
    end
    return ReferenceTessellation(coords, tess.triangles, edges)
end

# The tessellation a cell is actually instantiated with: `surface_rounds`
# subdivision rounds for the triangles, and at least `edge_rounds` rounds for
# the wireframe edges. Each surface round already halves the edge segments, so
# only the missing edge-only rounds are applied on top and every edge ends up
# split into 2^max(surface_rounds, edge_rounds) segments.
function _subdivided_tessellation(base::ReferenceTessellation, surface_rounds::Int, edge_rounds::Int)
    tess = subdivide(base, surface_rounds)
    for _ in (surface_rounds + 1):edge_rounds
        tess = _subdivide_edges(tess)
    end
    return tess
end

# Self-contained wireframe geometry of a tessellation: only the coordinates
# referenced by edge segments, reindexed from 1. Used by layouts that rebuild
# their vertices from scratch (AddQuadraturePointData) and still want the FE
# cell edges drawn.
function edge_geometry(tess::ReferenceTessellation{refdim,T}) where {refdim,T}
    remap = Dict{Int,Int}()
    coords = Ferrite.Vec{refdim,T}[]
    edges = NTuple{2,Int}[]
    function vertexof(i)
        return get!(remap, i) do
            push!(coords, tess.coords[i])
            length(coords)
        end
    end
    for (a, b) in tess.edges
        push!(edges, (vertexof(a), vertexof(b)))
    end
    return coords, edges
end

# Map a reference coordinate through the cell's geometric interpolation.
function geometric_map(ip_geo::Ferrite.ScalarInterpolation, node_coords, ξ)
    x = Ferrite.reference_shape_value(ip_geo, ξ, 1) * node_coords[1]
    for i in 2:Ferrite.getnbasefunctions(ip_geo)
        x += Ferrite.reference_shape_value(ip_geo, ξ, i) * node_coords[i]
    end
    return x
end

# Cell centroid mapped through the geometric interpolation, for label placement.
function midpoint(cell::Ferrite.AbstractCell, points::AbstractVector{PT}) where {dim,T,PT<:GeometryBasics.Point{dim,T}}
    ip_geo = Ferrite.geometric_interpolation(typeof(cell))
    refcoords = Ferrite.reference_coordinates(ip_geo)
    ξc = sum(refcoords) / length(refcoords)
    x = zero(PT)
    for (i, node) in enumerate(cell.nodes)
        x = x .+ T(Ferrite.reference_shape_value(ip_geo, ξc, i)) .* points[node]
    end
    return x
end

