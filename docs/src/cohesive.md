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

A four-node cohesive quad has its facets on nodes `(1,2)` and `(3,4)`. Drawing
it as a filled quad therefore means walking the node loop `1 → 2 → 4 → 3`
(not the standard quadrilateral loop `1 → 2 → 3 → 4`, which would produce a
bow-tie). We give it a dedicated reference shape so its tessellation is
independent of the ordinary quadrilateral:

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

# 3. the geometric interpolation: node i sits at reference corner PERM[i], so
# the facets (1,2) and (3,4) are the bottom and top edges of the reference
# quadrilateral and the bilinear map never folds
struct CohesiveLagrange <: Ferrite.ScalarInterpolation{RefCohesiveQuad,1} end
const PERM = (1, 2, 4, 3)
Ferrite.getnbasefunctions(::CohesiveLagrange) = 4
Ferrite.reference_shape_value(::CohesiveLagrange, ξ, i::Int) =
    Ferrite.reference_shape_value(Lagrange{RefQuadrilateral,1}(), ξ, PERM[i])
Ferrite.reference_coordinates(::CohesiveLagrange) =
    Ferrite.reference_coordinates(Lagrange{RefQuadrilateral,1}())[collect(PERM)]

Ferrite.geometric_interpolation(::Type{CohesiveQuadrilateral}) = CohesiveLagrange()
nothing # hide
```

The interpolation is the important design decision. Plugging the cohesive node
order into a plain `Lagrange{RefQuadrilateral,1}` would trace the loop
`1 → 2 → 3 → 4` and fold the bilinear map into a bow-tie: the corners land in
the right places, but the interior of the map self-intersects. The static
tessellation would never notice (it only ever samples the corners), but
everything that evaluates the geometry *between* the corners — the
error-adaptive refinement in particular — would faithfully render the fold.
Permuting the nodes onto the reference corners keeps the map clean everywhere.

## The one method FerriteViz needs

With a clean parametrization the reference geometry *is* just a quadrilateral,
so its tessellation — surface triangles, and the element edges that give
[`meshplot`](@ref) its wireframe — can simply be reused:

```@example cohesive
FerriteViz.reference_tessellation(::Type{RefCohesiveQuad}) =
    FerriteViz.reference_tessellation(Ferrite.RefQuadrilateral)
nothing # hide
```

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
