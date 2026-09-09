# Custom cells: cohesive/interface elements

FerriteViz can render *any* cell type, not just the ones Ferrite ships with. A
cell's rendering geometry is described entirely by the
[`FerriteViz.ReferenceTessellation`](@ref) of its reference shape, so teaching
FerriteViz about a new cell means implementing a single method,
[`FerriteViz.reference_tessellation`](@ref) (see the
[developer documentation](devdocs.md) for the full extension-point contract).

This page walks through a complete example: **cohesive zone (interface)
elements**, as used for modelling delamination and cohesive fracture (e.g. the
[FerriteCohesiveZones.jl](https://github.com/kimauth/FerriteCohesiveZones.jl)
package). A cohesive element has two facets that *coincide* in the undeformed
configuration and *separate* ("open") under load. That makes its node
numbering unusual — the two facets are listed one after another — so it needs
its own reference shape, parametrized so the cell renders as a
non-self-intersecting quad.

The whole example only depends on `Ferrite` and `FerriteViz`.

```@example cohesive
import WGLMakie, Bonito # hide
Bonito.Page() # hide
WGLMakie.activate!() # hide
WGLMakie.Makie.inline!(true) # hide
nothing # hide
```

## Defining the cell

A four-node cohesive quad has its facets on nodes `(1,2)` and `(3,4)`, so its
node loop is `1 → 2 → 4 → 3` rather than the quadrilateral's `1 → 2 → 3 → 4`.
We give it a dedicated reference shape, and a geometric interpolation that
places each node at the reference corner it belongs to:

```@example cohesive
using Ferrite
import FerriteViz
using FerriteViz: FEData
import WGLMakie

# 1. a reference shape for the interface element
struct RefCohesiveQuad <: Ferrite.AbstractRefShape{2} end

# 2. the cell: four nodes, two facets (1,2) and (3,4) that coincide when closed
struct CohesiveQuadrilateral <: Ferrite.AbstractCell{RefCohesiveQuad}
    nodes::NTuple{4,Int}
end

# 3. the geometric interpolation: facet (1,2) is the bottom edge η = -1 of the
# reference square and facet (3,4) the top edge η = +1, so nodes 3 and 4 sit
# at the corners (-1, 1) and (1, 1) — the bilinear map then never folds
struct CohesiveLagrange <: Ferrite.ScalarInterpolation{RefCohesiveQuad,1} end
Ferrite.getnbasefunctions(::CohesiveLagrange) = 4
function Ferrite.reference_shape_value(::CohesiveLagrange, ξ::Vec{2}, i::Int)
    x, y = ξ
    i == 1 && return (1 - x) * (1 - y) / 4     # node 1 at (-1, -1)
    i == 2 && return (1 + x) * (1 - y) / 4     # node 2 at ( 1, -1)
    i == 3 && return (1 - x) * (1 + y) / 4     # node 3 at (-1,  1)
    i == 4 && return (1 + x) * (1 + y) / 4     # node 4 at ( 1,  1)
    throw(ArgumentError("no shape function $i"))
end
Ferrite.reference_coordinates(::CohesiveLagrange) =
    [Vec((-1.0, -1.0)), Vec((1.0, -1.0)), Vec((-1.0, 1.0)), Vec((1.0, 1.0))]

Ferrite.geometric_interpolation(::Type{CohesiveQuadrilateral}) = CohesiveLagrange()
nothing # hide
```

The interpolation is the important design decision. Plugging the cohesive node
order into a plain `Lagrange{RefQuadrilateral,1}` would put node 3 at the
corner `(1, 1)` and node 4 at `(-1, 1)`, and the bilinear map would fold into
a bow-tie: the corners land in the right places, but the interior of the map
degenerates. The static tessellation would never notice (it only ever samples
the corners and the centre), but everything that evaluates the geometry
*between* the corners — the error-adaptive refinement in particular — would
faithfully render the fold. Placing the nodes onto their proper reference
corners keeps the map clean everywhere.

## The one method FerriteViz needs

A [`FerriteViz.ReferenceTessellation`](@ref) has three parts, all in
*reference* space: the vertices, the triangles indexing into them, and the
wireframe edges drawn by [`meshplot`](@ref). Vertex indices refer to
positions in the reference square, not to cell nodes — the interpolation
above is what ties the two together. The four corners plus a centre vertex
give a fan of four triangles; the four sides are the cell edges:

```@example cohesive
function FerriteViz.reference_tessellation(::Type{RefCohesiveQuad})
    coords = [Vec((-1.0, -1.0)), Vec((1.0, -1.0)), Vec((1.0, 1.0)), Vec((-1.0, 1.0)),  # corners
              Vec((0.0, 0.0))]                                                          # centre, vertex 5
    triangles = [(1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 1, 5)]   # fan around the centre
    edges = [(1, 2), (2, 3), (3, 4), (4, 1)]                   # the wireframe
    return FerriteViz.ReferenceTessellation(coords, triangles, edges)
end
nothing # hide
```

Through the interpolation, reference edge `(2, 3)` runs from node 2 to node 4
and edge `(4, 1)` from node 3 to node 1 — the cohesive `1 → 2 → 4 → 3` loop,
without any renumbering in the tessellation itself.

The tessellation is written out here to show what the extension point consists
of. Whenever a custom shape's reference geometry coincides with one FerriteViz
already knows — as it does here, since the cohesive cell is parametrized as a
plain quadrilateral — the method can just as well forward to the existing
tessellation:

```julia
FerriteViz.reference_tessellation(::Type{RefCohesiveQuad}) =
    FerriteViz.reference_tessellation(Ferrite.RefQuadrilateral)
```

Both definitions are equivalent; the explicit one is the template for shapes
that have no built-in counterpart.

That is the entire extension. `FEData` and every representation now work for
grids containing `CohesiveQuadrilateral`s — including the wireframe and the
error-adaptive refinement, whose conforming base is fanned over exactly the
element edges the tessellation lists. (A tessellation may also omit the edge
list; such cells then draw no wireframe and their datasets keep the static
tessellation.)

## A tiny grid with an interface

Two square blocks share an interface bridged by a single cohesive element. The
parameter `Δ` shifts the right block, opening the interface — the cohesive cell
goes from zero thickness (closed) to a real quad of width `Δ` (opened). We tag
each cell with its opening as per-cell data:

```@example cohesive
function cohesive_demo(Δ)
    nodes = [
        Node((0.0, 0.0)), Node((1.0, 0.0)), Node((1.0, 1.0)), Node((0.0, 1.0)),          # left block
        Node((1.0 + Δ, 0.0)), Node((1.0 + Δ, 1.0)), Node((2.0 + Δ, 1.0)), Node((2.0 + Δ, 0.0)), # right block
    ]
    cells = Ferrite.AbstractCell[
        Quadrilateral((1, 2, 3, 4)),
        Quadrilateral((5, 6, 7, 8)),
        CohesiveQuadrilateral((2, 3, 5, 6)),   # facet (2,3) on the left block, (5,6) on the right
    ]
    grid = Grid(cells, nodes)
    ds = FEData(DofHandler(grid), Float64[])   # no dof field needed, we colour by cell data
    FerriteViz.set_cell_data!(ds, :opening, [0.0, 0.0, Δ])
    return ds
end
nothing # hide
```

## Plotting it

```@example cohesive
f = WGLMakie.Figure(size = (1000, 380))
ax1 = WGLMakie.Axis(f[1, 1], aspect = WGLMakie.DataAspect(), title = "closed")
ax2 = WGLMakie.Axis(f[1, 2], aspect = WGLMakie.DataAspect(), title = "opened (Δ = 0.5)")

# closed configuration: the interface is a line between the two blocks
FerriteViz.meshplot!(ax1, cohesive_demo(0.0))

# opened configuration: the cohesive cell is now a real quad, coloured by its opening
opened = cohesive_demo(0.5)
p = FerriteViz.cellplot!(ax2, opened; color = :opening, colormap = :inferno)
FerriteViz.meshplot!(ax2, opened; plotnodes = false)
WGLMakie.Colorbar(f[1, 3], p, label = "interface opening")

f
```

The cohesive cell tessellates, colours and refines just like a native Ferrite
cell — FerriteViz never needed to know what a cohesive element *is*, only how
its reference shape is parametrized and triangulated.
