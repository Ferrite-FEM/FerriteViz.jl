# Tutorial

FerriteViz follows the ParaView model: a **source** holds your solution, **filters**
transform it, and **representations** draw the result.

```
FEData(dh, u)  |>  filter  |>  filter  |>  solutionplot(...)
   source              transformations          representation
```

This tutorial introduces the representations first, then shows how filters are chained
to get from a raw solution vector to derived quantities like stresses, pressures and
internal variables.

```@example 1
import WGLMakie, Bonito # hide
Bonito.Page() # hide
WGLMakie.activate!() # hide
WGLMakie.Makie.inline!(true) # hide
nothing # hide
```

!!! tip "Plotting functions"
    [`FerriteViz.solutionplot`](@ref), [`FerriteViz.meshplot`](@ref),
    [`FerriteViz.surfaceplot`](@ref), [`FerriteViz.arrowplot`](@ref),
    [`FerriteViz.cellplot`](@ref) and their mutating analogues with `!` are defined for
    `FEData`. The docs need `WGLMakie`; locally you can replace every `WGLMakie` call by
    `GLMakie` (or `CairoMakie` for 2D).

## Plotting recipes

### The mesh

Solve a boundary value problem as you usually would with Ferrite, and keep the
`DofHandler` and the solution vector — those are what `FEData` needs. Before that, a grid
alone is already plottable:

```@example 1
import FerriteViz
using FerriteViz: FEData, WarpByVector, Gradient, Derive, AddQuadraturePointData,
                  CrinkleClip, ClipPlane, vonmises
using Ferrite
import WGLMakie #activating the backend, switch to GLMakie or CairoMakie (for 2D) locally
WGLMakie.set_theme!(size=(800, 400)) # hide

grid = generate_grid(Hexahedron,(3,3,3))
FerriteViz.meshplot(grid,markersize=10,linewidth=2)
```

Node and cell labels as well as cellsets can be shown too, which is handy while debugging
a mesh:

```@example 1
grid = generate_grid(Quadrilateral,(3,3))
addcellset!(grid,"s1",Set((1,4,7)))
addcellset!(grid,"s2",Set((2,5,8)))
addcellset!(grid,"s3",Set((3,6,9)))
FerriteViz.meshplot(grid,markersize=10,linewidth=1,nodelabels=true,celllabels=true,cellsets=true)
```

Curved (higher-order geometry) cells are rendered curved: the cell surfaces and the
wireframe edges are subdivided in reference space and mapped through the geometric
interpolation, so the element edges bend through their midside nodes instead of being
drawn as straight chords. A quarter annulus of quadratic quadrilaterals:

```@example 1
grid = generate_grid(QuadraticQuadrilateral,(6,3))
Ferrite.transform_coordinates!(grid, x -> begin
    r = 1.5 + 0.5x[2]        # x[2] ∈ [-1,1]  →  r ∈ [1,2]
    θ = π/4*(x[1] + 1)       # x[1] ∈ [-1,1]  →  θ ∈ [0,π/2]
    Vec(r*cos(θ), r*sin(θ))
end)
FerriteViz.meshplot(grid, markersize=8, linewidth=2, axis=(aspect=WGLMakie.DataAspect(),))
```

The subdivision kicks in automatically whenever the geometry *or* any field of the dof
handler is nonlinear — a quadratic displacement field warping a linear mesh bends the
wireframe just the same. How fine the surfaces and edges are resolved is controlled by
the `resolution` and `edge_resolution` keywords of [`FEData`](@ref)
(`FEData(dh, u; resolution=0, edge_resolution=0)` restores the flat tessellation).

### The solution field

[`FEData`](@ref) wraps a `DofHandler` together with a solution vector; every plotting
function and every filter operates on it. Here we use a mixed displacement/pressure
formulation of incompressible elasticity, which gives us two fields to play with.

```@example 1
include("ferrite-examples/incompressible-elasticity.jl") #defines dh_quadratic, u_quadratic and mp

ds = FEData(dh_quadratic,u_quadratic)
FerriteViz.arrowplot(ds)
```

By default every plot grabs the first field of the `DofHandler`, reduced to its magnitude
if it is vector valued. Naming a field colors by it instead — here the pressure:

```@example 1
FerriteViz.solutionplot(ds,color=:p)
```

For 2D problems the solution can also be lifted into the third dimension with
`surfaceplot`. The mutating variants let you combine representations in one scene:

```@example 1
FerriteViz.surfaceplot(ds)
FerriteViz.solutionplot!(ds,colormap=:magma)
WGLMakie.current_figure()
```

### Per-cell data

Quantities that live per cell rather than per node are registered with
[`FerriteViz.set_cell_data!`](@ref) and drawn with `cellplot`. We use the plastified
cantilever from the plasticity example, which returns the cell averaged von Mises stress:

```@example 1
include("ferrite-examples/plasticity.jl") #only defines the solving function
u, dh, u_history, mises, κ, states, qr = solve(; celltype = Hexahedron)
ds_p = FEData(dh,u)

FerriteViz.cellplot(ds_p,mises,colormap=:thermal)
WGLMakie.current_figure()
```

## Chaining filters

Filters map an `FEData` to a new `FEData` and are applied by piping. Because the result is
again an `FEData`, they compose freely, and the whole chain stays reactive: a
[`FerriteViz.update!`](@ref) anywhere in the pipeline propagates into every open plot.

### Deforming the geometry

[`WarpByVector`](@ref) displaces the geometry by a vector field, so everything drawn
downstream of it appears in the deformed configuration:

```@example 1
warped = ds_p |> WarpByVector(:u, 2.0)
FerriteViz.cellplot(warped,mises,colormap=:thermal)
FerriteViz.meshplot!(warped,markersize=10,linewidth=1)
WGLMakie.current_figure()
```

### Looking inside a 3D domain

[`CrinkleClip`](@ref) hides cells on one side of a decision function, revealing the
interior. Filters chain, so we clip and then warp:

```@example 1
clipped = ds_p |> CrinkleClip(ClipPlane(Vec((0.0,0.5,0.5)), 0.7)) |> WarpByVector(:u, 2.0)
FerriteViz.solutionplot(clipped,colormap=:thermal)
WGLMakie.current_figure()
```

The plane can be replaced by any function of the grid and a cell index returning whether
the cell stays visible.

### Derived fields: `Gradient` into `Derive`

This is where chaining pays off. [`Gradient`](@ref) turns a field into its piecewise
discontinuous gradient — for a displacement field that is ∇u, from which strains and
stresses follow. [`Derive`](@ref) then maps that array through an arbitrary function,
which is where the constitutive law goes.

Back to the mixed formulation: `copy_fields` carries the pressure along, so the pressure
of the mixed formulation and the stress derived from the displacement gradient live in
*one* pipeline:

```@example 1
ε(∇u) = symmetric(∇u)
stress(∇u) = 2*mp.G*dev(ε(∇u)) + mp.K*tr(ε(∇u))*one(ε(∇u))

mixed = FEData(dh_quadratic, u_quadratic) |>
        Gradient(:u; copy_fields = [:p]) |>
        Derive(∇u -> vonmises(stress(∇u)); input = :gradient, output = :σvM)

f = WGLMakie.Figure(size = (900, 330))
ax1 = WGLMakie.Axis(f[1,1], title = "von Mises stress (derived)")
ax2 = WGLMakie.Axis(f[1,3], title = "pressure (solved for)")
p1 = FerriteViz.solutionplot!(ax1, mixed; color = :σvM, colormap = :jet)
p2 = FerriteViz.solutionplot!(ax2, mixed; color = :p,   colormap = :jet)
WGLMakie.Colorbar(f[1,2], p1)
WGLMakie.Colorbar(f[1,4], p2)
f
```

After `Gradient` the pipeline's primary field is `:gradient`, which is why `Derive` is
told `input = :gradient` explicitly (that is also its default via `:default`). `Derive`
can be chained repeatedly to build up several named arrays from the same gradient.

### Quadrature point data: `AddQuadraturePointData` into `Derive`

Internal variables such as plastic strain or stress are only known at the quadrature
points. [`AddQuadraturePointData`](@ref) puts them on the mesh without averaging, by
partitioning each cell into the Voronoi regions of its quadrature points.

Because the filter output is ordinary point data, it feeds straight into `Derive` — and
`Derive` accepts *several* inputs, one argument per name. That lets us combine two
quadrature point quantities, here the stress and the plastic strain, into the plastic work
density ``\sigma : \varepsilon^\mathrm{p}``:

```@example 1
dissipation = ds_p |>
    AddQuadraturePointData(qr, states; extract = s -> s.σ,  output = :σ) |>
    AddQuadraturePointData(qr, states; extract = s -> s.ϵᵖ, output = :εᵖ) |>
    Derive((σ, εᵖ) -> σ ⊡ εᵖ; input = [:σ, :εᵖ], output = :wᵖ) |>
    WarpByVector(:u, 2.0)

FerriteViz.solutionplot(dissipation; color = :wᵖ, colormap = :inferno)
FerriteViz.meshplot!(ds_p |> WarpByVector(:u, 2.0); plotnodes = false, linewidth = 1)
WGLMakie.current_figure()
```

The element outlines come from the plain warped dataset and make the resolution visible:
each element is filled by several flat patches, one per quadrature point, rather than a
single averaged colour.

`extract` pulls the quantity out of the material state struct, so `states` can be
handed over as it comes out of your solver. The two `AddQuadraturePointData` filters share the
same quadrature rule and therefore the same vertex layout, which is what allows their
arrays to be combined afterwards — and you can chain as many of them as you have
quantities, then take them all into one `Derive`.

!!! note
    We use `celltype = Hexahedron` above on purpose: linear tetrahedra are constant strain
    elements, so all of their quadrature points carry the same state and the plot would be
    indistinguishable from a cell averaged one. The trilinear hexahedron has bilinear modes
    whose strains genuinely differ between quadrature points.

## What's next?

The [Recommended Practices](atopics.md) page goes deeper: what the discontinuous gradient buys
you over an L2 projection, how the quadrature point partition is built, high order fields,
the composable viewer and live plotting during a simulation.

This package also provides an interactive viewer, `ferriteviewer(ds)` and
`ferriteviewer(ds,u_history)` for time dependent views. Its layout and controls are built
from `Makie.SpecApi` and can be fully customized, see the
[composable viewer](atopics.md#Composable-viewer) section.
