# Recommended Practices

```@example 1
import WGLMakie, Bonito # hide
Bonito.Page() # hide
WGLMakie.activate!() # hide
WGLMakie.Makie.inline!(true) # hide
nothing # hide
```

## Why the gradient field is discontinuous

The [tutorial](tutorial.md) shows how `Gradient |> Derive` turns a displacement field into
stresses. What it does not explain is *why* [`Gradient`](@ref) interpolates into a
piecewise **discontinuous** field instead of doing what most codes do, an L2 projection
onto an H1 field.

The projection smooths the inter-element jumps away — and those jumps are one of the best
indicators of an under-resolved discretization. Keeping them makes the difference between
a linear and a quadratic discretization of the same problem immediately visible: the
linear one shows pronounced facets between elements, the quadratic one is nearly smooth.

```@example 1
using Ferrite
import FerriteViz
using FerriteViz: FEData, WarpByVector, Gradient, Derive, Refine, FirstOrderRefinement, CrinkleClip, ClipPlane
ε(∇u) = (∇u+transpose(∇u))/2
import WGLMakie #activating the backend, switch to GLMakie or CairoMakie (for 2D) locally

include("ferrite-examples/incompressible-elasticity.jl") #defines dh_linear, dh_quadratic, u_linear, u_quadratic and mp

σ(∇u) = 2*mp.G*dev(ε(∇u)) + mp.K*tr(ε(∇u))*one(ε(∇u)) #helper function to map gradient to stress
cmap = :jet

# pipelines: gradient field + named derived arrays
pipeline_linear    = FEData(dh_linear, u_linear)       |> Gradient(:u) |> Derive(∇u->norm(ε(∇u)), output=:εnorm) |> Derive(∇u->norm(σ(∇u)), output=:σnorm)
pipeline_quadratic = FEData(dh_quadratic, u_quadratic) |> Gradient(:u) |> Derive(∇u->norm(ε(∇u)), output=:εnorm) |> Derive(∇u->norm(σ(∇u)), output=:σnorm)
deformed_linear    = FEData(dh_linear, u_linear)       |> WarpByVector(:u)
deformed_quadratic = FEData(dh_quadratic, u_quadratic) |> WarpByVector(:u)

f = WGLMakie.Figure()
axs = [WGLMakie.Axis(f[1, 1], title="Strain norm (linear)"),WGLMakie.Axis(f[1, 2], title="Stress norm (linear)"),WGLMakie.Axis(f[1, 3], title="Pressure (deformed, linear)"),
       WGLMakie.Axis(f[3, 1], title="Strain norm (quadratic)"),WGLMakie.Axis(f[3, 2], title="Stress norm (quadratic)"),WGLMakie.Axis(f[3, 3], title="Pressure (deformed, quadratic)")]
p1 = FerriteViz.solutionplot!(axs[1], pipeline_linear, color=:εnorm, colormap=cmap)
p2 = FerriteViz.solutionplot!(axs[2], pipeline_linear, color=:σnorm, colormap=cmap)
p3 = FerriteViz.solutionplot!(axs[3], deformed_linear, color=:p, colormap=cmap)
f[2,1] = WGLMakie.Colorbar(f[1,1], p1, vertical=false)
f[2,2] = WGLMakie.Colorbar(f[1,2], p2, vertical=false)
f[2,3] = WGLMakie.Colorbar(f[1,3], p3, vertical=false)

p4 = FerriteViz.solutionplot!(axs[4], pipeline_quadratic, color=:εnorm, colormap=cmap)
p5 = FerriteViz.solutionplot!(axs[5], pipeline_quadratic, color=:σnorm, colormap=cmap)
p6 = FerriteViz.solutionplot!(axs[6], deformed_quadratic, color=:p, colormap=cmap)
f[4,1] = WGLMakie.Colorbar(f[3,1], p1, vertical=false)
f[4,2] = WGLMakie.Colorbar(f[3,2], p2, vertical=false)
f[4,3] = WGLMakie.Colorbar(f[3,3], p3, vertical=false)

f
```

For genuinely stress-valued arrays (e.g. registered per-cell data from quadrature-point states, see
[`FerriteViz.set_cell_data!`](@ref)), the [`VonMises`](@ref) and [`Deviator`](@ref) filters apply
those reductions directly.

An alternative to this approach is to compute gradient quantities at sample points and plot these via `arrowplot`.

## What the quadrature point partition buys you

The [tutorial](tutorial.md) uses [`AddQuadraturePointData`](@ref) to plot internal variables.
The reason it exists is that the two usual alternatives both destroy information: averaging
per cell throws away the sub-element variation, and projecting onto a nodal field invents
smoothness that smears exactly the localization one wants to see.

The filter instead partitions every cell into the **Voronoi regions of its quadrature
points** and fills each region with that point's value. A localisation band makes the
difference obvious — the cell average clips the peak, the partition resolves it:

```@example 1
using FerriteViz: AddQuadraturePointData
import WGLMakie

band(x) = exp(-((x[1] - x[2]) / 0.18)^2)   # a localisation band

grid_iv = generate_grid(Quadrilateral, (12, 12))
dh_iv = DofHandler(grid_iv); add!(dh_iv, :s, Lagrange{RefQuadrilateral,1}()); close!(dh_iv)
qr_iv = QuadratureRule{RefQuadrilateral}(2)
cv_iv = CellValues(qr_iv, Lagrange{RefQuadrilateral,1}(), Lagrange{RefQuadrilateral,1}())

# one value per (cell, quadrature point) — the layout of Ferrite's material states
qpvals = [zeros(getnquadpoints(qr_iv)) for _ in 1:getncells(grid_iv)]
for cell in CellIterator(dh_iv)
    reinit!(cv_iv, cell)
    coords = getcoordinates(cell)
    for q in 1:getnquadpoints(qr_iv)
        qpvals[cellid(cell)][q] = band(spatial_coordinate(cv_iv, q, coords))
    end
end

ds_iv = FEData(dh_iv, zeros(ndofs(dh_iv)))
resolved = ds_iv |> AddQuadraturePointData(qr_iv, qpvals; output = :iv)

averaged = FEData(dh_iv, zeros(ndofs(dh_iv)))
FerriteViz.set_cell_data!(averaged, :avg, [sum(v) / length(v) for v in qpvals])

f = WGLMakie.Figure(size = (900, 400))
ax1 = WGLMakie.Axis(f[1, 1], aspect = WGLMakie.DataAspect(), title = "cell averaged")
ax2 = WGLMakie.Axis(f[1, 2], aspect = WGLMakie.DataAspect(), title = "quadrature point Voronoi")
FerriteViz.cellplot!(ax1, averaged; color = :avg, colormap = :inferno, colorrange = (0, 1))
p = FerriteViz.solutionplot!(ax2, resolved; color = :iv, colormap = :inferno, colorrange = (0, 1))
WGLMakie.Colorbar(f[1, 3], p)
f
```

Beyond the `Vector` of per-cell vectors used above, values may be given as a `Matrix`
(`values[cell, qp]`) or as an `Observable` of either for live updating. Grids with mixed
cell types take one rule per reference shape, e.g.
`AddQuadraturePointData(Dict(RefTriangle => qr_tri, RefQuadrilateral => qr_quad), values)`.
Since the filter rebuilds the geometry, apply [`WarpByVector`](@ref) *after* it.

!!! note
    The region boundaries are a visualization choice, not physical discontinuities — they
    only mark where the closest quadrature point changes. Values are likewise extended
    from the (interior) quadrature points out to the element boundary, and in 3D you see
    the partition intersected with the surface of the domain.

## High-order fields

The investigation of high-order fields is currently only supported via a first-order refinement of the problem.
Here, the high-order approximation is replaced by a first order approximation of the field, which is
spanned by the nodes of the high-order approximation — the [`FirstOrderRefinement`](@ref) filter. For example, the first order refinement of a
heat problem on a square domain for Lagrange polynomials of order 4 looks like this:
```@example 1
include("ferrite-examples/heat-equation.jl"); #defines manufactured_heat_problem

f = WGLMakie.Figure()
axs = [WGLMakie.Axis3(f[1, 1], title="Coarse"), WGLMakie.Axis3(f[1, 2], title="Fine")]

dh,u = manufactured_heat_problem(Triangle, Lagrange{RefTriangle,4}(), 1)
FerriteViz.surfaceplot!(axs[1], FEData(dh, u) |> FirstOrderRefinement())

dh,u = manufactured_heat_problem(Triangle, Lagrange{RefTriangle,4}(), 3)
FerriteViz.surfaceplot!(axs[2], FEData(dh, u) |> FirstOrderRefinement())

f
```
Note that this method produces small artifacts due to the flattening of the nonlinearities of the high order ansatz.
However, it is still sufficient to investigate important features of the solution.
If users want to have higher resolution than the crude estimate given by the first order refinement (as well as enough RAM), then we also provide a uniform tessellation algorithm, the [`Refine`](@ref) filter:
```@example 1
f = WGLMakie.Figure()
axs = [WGLMakie.LScene(f[1, 1]), WGLMakie.LScene(f[1, 2])]

dh, u = manufactured_heat_problem(Hexahedron, Lagrange{RefHexahedron,2}(), 2);
clipped = FEData(dh,u) |> CrinkleClip(ClipPlane(Ferrite.Vec((0.0,0.5,0.5)), 0.1));

FerriteViz.solutionplot!(axs[1], clipped)
FerriteViz.solutionplot!(axs[2], clipped |> Refine(4))

f
```

In future we will also provide an adaptive tessellation algorithm to resolve the high-order fields with full detail.

## Pipeline semantics

Chaining itself is covered in the [tutorial](tutorial.md); two rules matter once pipelines
get longer.

**Ordering.** Geometry-rebuilding filters ([`Refine`](@ref),
[`FirstOrderRefinement`](@ref), [`AddQuadraturePointData`](@ref)) rebuild from the base
geometry, so apply [`WarpByVector`](@ref) *after* them. They also drop the point data
registered upstream, since it refers to vertices that no longer exist — the one exception
is a rebuild that reproduces the very same vertex layout, as when two
[`AddQuadraturePointData`](@ref) share a quadrature rule, which is what lets their arrays be
combined in a later [`Derive`](@ref). [`CrinkleClip`](@ref) and [`Gradient`](@ref) share
the — possibly warped — geometry of their input, so a warp survives those.

**Reactivity.** A [`FerriteViz.update!`](@ref) on *any* dataset of a pipeline updates the
root solution and propagates through every filter into all open plots, which is what makes
live plotting work. Filters can also be applied explicitly with
[`FerriteViz.apply`](@ref) instead of `|>`.

## Composable viewer

`ferriteviewer(ds)` gives a single panel with the usual menus and toggles, and
`ferriteviewer(ds, u_history)` adds a slider stepping through a solution history. Neither
is hard-wired though — the whole viewer is described declaratively with `Makie.SpecApi`,
so you can replace the layout and the controls.

### How SpecApi fits in

`Makie.SpecApi` describes a figure as **data** rather than by mutating a scene: a
`PlotSpec` is "this plot type, with these arguments and attributes", a `BlockSpec` an axis
or colorbar, and a `GridLayoutSpec` says how they are arranged. Building a spec draws
nothing. When Makie is handed a *new* spec it diffs it against the previous one and
updates only what actually changed, instead of rebuilding the scene. See
[Makie's documentation](https://docs.makie.org) for the full `SpecApi` reference.

FerriteViz plugs into that with three pieces:

* `layout(ds, state) -> GridLayoutSpec` — called again whenever the view state changes,
* [`FerriteViz.Control`](@ref)s — the widgets that feed `state`,
* spec helpers — [`solutionplotspec`](@ref) and friends wrap a representation into a
  `PlotSpec`, [`panelspec`](@ref) puts plots into an axis (or `LScene` for 3D) with an
  optional linked colorbar.

### A two-panel viewer

The mixed formulation solves for a displacement *and* a pressure, so a natural viewer
shows both at once. `state.colormap` comes from the colormap menu and drives both panels;
the deformation toggle feeds the pipeline, so it warps both as well:

```@example 1
using FerriteViz: FEData, ferriteviewer, panelspec, solutionplotspec, S,
                  ColormapMenu, DeformationToggle
import WGLMakie

function twopanels(ds, state)
    disp = solutionplotspec(ds; color = :default, colormap = state.colormap)
    pres = solutionplotspec(ds; color = :p,       colormap = state.colormap)
    return S.GridLayout([
        panelspec(disp; colorbar = disp, dim = 2, axis = (; title = "displacement magnitude"))
        panelspec(pres; colorbar = pres, dim = 2, axis = (; title = "pressure"))
    ])
end

ferriteviewer(FEData(dh_quadratic, u_quadratic);
              layout = twopanels, controls = [ColormapMenu(), DeformationToggle()])
```

Only the controls a layout actually consumes need to be listed — the default set is just
one such list. Writing your own control means returning the widgets together with the
state they contribute:

```julia
using FerriteViz: Control, ControlResult
scale = Control() do fig, ds
    slider = WGLMakie.Slider(fig, range = 0:0.5:5)
    ControlResult(Any[slider]; structural = [:scale => slider.value])
end
```

Structural state (colormap, which panels, anything the layout reads) re-diffs the spec, so
Makie updates only the changed attributes of the plots it reuses. Data streaming —
[`FerriteViz.update!`](@ref) and the deformation scale — bypasses the spec entirely and
mutates the shared GPU buffers instead, which is why live plotting stays cheap.

## Live plotting

Plotting while a computationally heavy simulation is performed can be easily achieved with FerriteViz.jl.
Every [`FEData`](@ref) holds the solution vector as an `Observable`.
If an `Observable` changes, all its dependencies are triggered to change as well — through the whole filter pipeline into the open plots.
For this purpose the function [`FerriteViz.update!`](@ref) is provided. It takes a `ds::FEData`
and a new solution vector `u_new` and updates the solution observable, thereby all open plots called with datasets of the pipeline are updated.

A summary of the needed steps for live plotting:
1. Create a `FEData` before your time stepping begins
2. Call a plot or the `ferriteviewer` and save the return in a variable, e.g. `fig`
3. `display(fig)` in order to force the plot/viewer to pop up, even if its called inside a function body
4. `FerriteViz.update!(ds,u_new)` where `u_new` corresponds to your new solution of the time step

As an illustrative example, let's consider a slightly modified [plasticity example of Ferrite.jl](https://github.com/Ferrite-FEM/FerriteViz.jl/blob/master/docs/src/ferrite-examples/plasticity-live.jl).
For the full source code, please refer to the link. In the following code we only highlight the necessary changes.

```julia
function solve(liveplotting=false)
    # set up your problem
    # lots of code
    dh = create_dofhandler(grid, interpolation) #helper function from script file
    n_dofs = ndofs(dh)  # total number of dofs
    u  = zeros(n_dofs)

    if liveplotting
        ####### Here we take care of the conceptual steps 1, 2 and 3 #######
        ds = FEData(dh,u)
        fig = ferriteviewer(ds)
        display(fig)
        ####################################################################
    end

    Δu = zeros(n_dofs)  # displacement correction
    r = zeros(n_dofs)   # residual
    K = allocate_matrix(dh); # tangent stiffness matrix

    nqp = getnquadpoints(cellvalues)
    states = [[MaterialState() for _ in 1:nqp] for _ in 1:getncells(grid)]

    # Newton-Raphson loop
    NEWTON_TOL = 1 # 1 N

    for timestep in 1:n_timesteps
        while true; newton_itr += 1

            if newton_itr > 8
                error("Reached maximum Newton iterations, aborting")
                break
            end
            K, r = doassemble(cellvalues, facevalues, K, grid, dh, material, u,
                             states, traction);
            norm_r = norm(r[Ferrite.free_dofs(dbcs)])

            if norm_r < NEWTON_TOL
                break
            end

            apply_zero!(K, r, dbcs)
            Δu = Symmetric(K) \ r
            u -= Δu
        end

        if liveplotting
            ####### Step 4 updating the current solution vector in ds #######
            FerriteViz.update!(ds,u)
            #################################################################
            sleep(0.1)
        end

        # Update all the material states after we have reached equilibrium
        for cell_states in states
            foreach(update_state!, cell_states)
        end
        u_max[timestep] = max(abs.(u)...) # maximum displacement in current timestep
    end

    # postprocessing
    # lots of code
    return u, dh, traction_magnitude
end

u, dh, traction_magnitude = solve();
```

Note that we create the `ds::FEData` object before the time stepping begins, as well as calling `ferriteviewer` on it.
The next function call is crucial to get the live plotting working. `display(fig)` forces the viewer to pop up, even if it's inside a function body.
Now, the only missing piece is the `FerriteViz.update!` of the dataset, which happens directly after the Newton iteration. The result for this code looks like this:

![liveplot](https://github.com/Ferrite-FEM/FerriteViz.jl/blob/master/docs/src/assets/liveplotting.gif?raw=true)

Since the computational load of one time step is in this example too low, the plotter would just update all the time and likely never display something, so we artificially increase the load of one time step by
`sleep`ing for 0.1s.

If you don't need the full viewer as a live plot, you can of course call instead `solutionplot` (or any other plot/plot combination) with appropriate keyword arguments to only have a specific live plot.
This can be beneficial performancewise.
