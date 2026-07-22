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
"""
struct ReferenceTessellation{refdim,T}
    coords::Vector{Ferrite.Vec{refdim,T}}
    triangles::Vector{NTuple{3,Int}}
end

nvertices(tess::ReferenceTessellation) = length(tess.coords)
ntriangles(tess::ReferenceTessellation) = length(tess.triangles)

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
    return ReferenceTessellation(coords, NTuple{3,Int}[])
end

function reference_tessellation(::Type{Ferrite.RefTriangle})
    coords = Ferrite.reference_coordinates(Ferrite.Lagrange{Ferrite.RefTriangle,1}())
    return ReferenceTessellation(coords, [(1, 2, 3)])
end

# The quadrilateral deliberately decomposes into 4 triangles through the center
# vertex: the cheaper 2-triangle split misses a (bilinear) solution mode.
function reference_tessellation(::Type{Ferrite.RefQuadrilateral})
    coords = Ferrite.reference_coordinates(Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}())
    push!(coords, zero(eltype(coords)))
    return ReferenceTessellation(coords, [(1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 1, 5)])
end

"""
    facet_based_tessellation(::Type{<:Ferrite.AbstractRefShape{3}}) -> ReferenceTessellation

Build the surface tessellation of a 3D reference shape from its faces: each
face's 2D tessellation (triangle or quadrilateral, chosen by vertex count) is
mapped into the element via `Ferrite.facet_to_element_transformation`. This is
the default building block for `reference_tessellation` of volumetric shapes.
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
    return ReferenceTessellation(coords, triangles)
end

reference_tessellation(::Type{Ferrite.RefTetrahedron}) = facet_based_tessellation(Ferrite.RefTetrahedron)
reference_tessellation(::Type{Ferrite.RefHexahedron}) = facet_based_tessellation(Ferrite.RefHexahedron)
reference_tessellation(::Type{Ferrite.RefPrism}) = facet_based_tessellation(Ferrite.RefPrism)
reference_tessellation(::Type{Ferrite.RefPyramid}) = facet_based_tessellation(Ferrite.RefPyramid)

"""
Number of triangles a cell tessellates into.
"""
ntriangles(cell::Ferrite.AbstractCell) = ntriangles(reference_tessellation(Ferrite.getrefshape(cell)))
nvertices_tess(cell::Ferrite.AbstractCell) = nvertices(reference_tessellation(Ferrite.getrefshape(cell)))

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

#####################################################################
# First-order refinement registry (consumed by FirstOrderRefinement)
#####################################################################

"""
    first_order_subcells(ip::Ferrite.Interpolation)

Connectivity of the first-order sub-cells spanned by the nodes of `ip`, used to
build a low-order discretization of a high-order field. Keyed by interpolation
(not just refshape/order) since the node layout is the interpolation's.
This is a second, optional extension point: rendering a custom cell only needs
[`reference_tessellation`](@ref); first-order refinement additionally needs this.
"""
function first_order_subcells end

"""
    linear_celltype(::Type{<:Ferrite.AbstractRefShape})

The linear `Ferrite.AbstractCell` type for a reference shape, used when
constructing first-order refined grids.
"""
linear_celltype(::Type{Ferrite.RefTriangle}) = Ferrite.Triangle
linear_celltype(::Type{Ferrite.RefQuadrilateral}) = Ferrite.Quadrilateral
linear_celltype(::Type{Ferrite.RefTetrahedron}) = Ferrite.Tetrahedron
linear_celltype(::Type{Ferrite.RefHexahedron}) = Ferrite.Hexahedron

# Triangle
first_order_subcells(::Union{Lagrange{RefTriangle,1},DiscontinuousLagrange{RefTriangle,1}}) = (
    (3,1,2),
)
# Quadratic Triangle
first_order_subcells(::Union{Lagrange{RefTriangle,2},DiscontinuousLagrange{RefTriangle,2}}) = (
    (6,1,4),
    (5,6,4),
    (3,6,5),
    (5,4,2),
)
# Cubic Triangle
first_order_subcells(::Union{Lagrange{RefTriangle,3},DiscontinuousLagrange{RefTriangle,3}}) = (
    (3,8,7),
    (7,8,10),
    (8,9,10),
    (10,9,4),
    (9,1,4),
    (7,10,6),
    (6,10,5),
    (6,5,2),
    (10,4,5),
)
# Biquadratic Triangle
first_order_subcells(::Union{Lagrange{RefTriangle,4},DiscontinuousLagrange{RefTriangle,4}}) = (
    (3,10,9),
    (13,9,10),
    (10,11,13),
    (14,13,11),
    (11,12,14),
    (4,14,12),
    (12,1,4),
    (9,13,8),
    (15,8,13),
    (13,14,15),
    (5,15,14),
    (14,4,5),
    (8,15,7),
    (6,7,15),
    (15,5,6),
    (7,6,2),
)
# Quintic Triangle
first_order_subcells(::Union{Lagrange{RefTriangle,5},DiscontinuousLagrange{RefTriangle,5}}) = (
    (3,12,11),
    (16,11,12),
    (12,13,16),
    (17,16,13),
    (13,14,17),
    (18,17,14),
    (14,15,18),
    (4,18,15),
    (15,1,4),
    (11,16,10),
    (19,10,16),
    (16,17,19),
    (20,19,17),
    (17,18,20),
    (5,20,18),
    (18,4,5),
    (10,19,9),
    (21,9,19),
    (19,20,21),
    (6,21,20),
    (20,5,6),
    (9,21,8),
    (7,8,21),
    (21,6,7),
    (8,7,2),
)
# Tetrahedron
first_order_subcells(::Union{Lagrange{RefTetrahedron,1},DiscontinuousLagrange{RefTetrahedron,1}}) = (
    (1,2,3,4),
)
# Quadratic Tetrahedron
first_order_subcells(::Union{Lagrange{RefTetrahedron,2},DiscontinuousLagrange{RefTetrahedron,2}}) = (
    (5,2,6,9),
    (7,6,3,10),
    (8,9,10,4),
    (8,5,6,9),
    (8,6,7,10),
    (5,8,1,6),
    (7,6,1,8),
    (9,10,8,6),
)
# Quadrilateral
first_order_subcells(::Union{Lagrange{RefQuadrilateral,1},DiscontinuousLagrange{RefQuadrilateral,1}}) = (
    (1,2,3,4),
)
# Quadratic Quadrilateral
first_order_subcells(::Union{Lagrange{RefQuadrilateral,2},DiscontinuousLagrange{RefQuadrilateral,2}}) = (
    (1,5,9,8),
    (5,2,6,9),
    (9,6,3,7),
    (8,9,7,4),
)
# Hexahedron
first_order_subcells(::Union{Lagrange{RefHexahedron,1},DiscontinuousLagrange{RefHexahedron,1}}) = (
    (1,2,3,4,5,6,7,8),
)
# Quadratic Hexahedron
first_order_subcells(::Union{Lagrange{RefHexahedron,2},DiscontinuousLagrange{RefHexahedron,2}}) = (
    (1,9,21,12,17,22,27,25),
    (17,22,27,25,5,13,26,16),
    (9,2,10,21,22,18,23,27),
    (22,18,23,27,13,6,14,26),
    (12,21,11,4,25,27,24,20),
    (25,27,24,20,16,26,15,8),
    (21,10,3,11,27,23,19,24),
    (27,23,19,24,26,14,7,15),
)
