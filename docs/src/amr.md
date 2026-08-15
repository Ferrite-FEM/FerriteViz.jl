# Adaptive grids (AMR)

FerriteViz renders adaptively refined grids — Ferrite's `ForestBWG`
forest-of-octrees and the `NonConformingGrid`s it materializes — and follows
them **live** through refinement: one dataset, one set of plots, updated
across the whole AMR loop.

Two pieces make that work:

  * **Hanging nodes draw watertight.** The adaptive tessellation splits every
    cell rim at the grid's hanging nodes, so both sides of a 2:1 interface
    tile identically — no cracks, no z-fighting, and inter-element jumps of
    the (DG-rendered) field stay exactly on the element boundaries.
  * **[`FerriteViz.regrid!`](@ref)** swaps a dataset's grid world — the new
    dof handler and solution after a `refine!`/`balanceforest!`/`creategrid`
    step — while the dataset and its plots keep their identity. Adaptive
    plots rebuild their tessellation against the new grid on their next
    update; between adaptations, `FerriteViz.update!` costs exactly what it
    costs on a fixed grid.

This page runs the heat problem from Ferrite's
[adaptive heat tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/heat_adaptivity/)
— a manufactured Gaussian-ring solution driving refinement — with the ring
centered on a corner of the domain so the feature is visible from outside.

## The 3D example, as a film

Each frame is one AMR sweep: solve, estimate, mark, `refine!`,
`balanceforest!`, `creategrid`, rebuild the dof handler, `regrid!` — and the
same `solutionplot` and adaptive `meshplot` follow from 64 to ~17k elements
with ~10k hanging nodes, watertight throughout.

```@raw html
<video autoplay loop muted playsinline controls src="../assets/amr_heat.mp4" style="max-width: 100%; border-radius: 4px;"></video>
```

The generating script is
[`docs/src/ferrite-examples/heat-amr-video.jl`](https://github.com/Ferrite-FEM/FerriteViz.jl/blob/master/docs/src/ferrite-examples/heat-amr-video.jl)
— run it locally with GLMakie to reproduce the video.

## The loop, live in 2D

The same structure, small enough to run during the documentation build: the
2D Gaussian ring on a quadtree forest. First the solver ingredients — note
the [`Ferrite.ConformityConstraint`](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/amr/),
which constrains the hanging dofs; the rendered field is whatever `u` says,
so continuity across 2:1 interfaces is the constraint handler's job, exactly
as in the solve itself.

```@example amr
import WGLMakie # hide
WGLMakie.Makie.inline!(true) # hide
using Ferrite, FerriteViz, LinearAlgebra, SparseArrays

const C0 = Ferrite.Vec(1.0, 1.0)
u_exact(x) = exp(-((norm(x - C0) - 0.7) / 0.12)^2)
f_rhs(x) = -laplace(u_exact, x)     # Tensors.laplace, reexported by Ferrite

function solve_on(grid)
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefQuadrilateral, 1}())
    close!(dh)
    ch = ConstraintHandler(dh)
    add!(ch, Ferrite.ConformityConstraint(:u))
    add!(ch, Dirichlet(:u, union((getfacetset(grid, f) for f in
        ("left", "right", "top", "bottom"))...), u_exact))
    close!(ch)
    qr = QuadratureRule{RefQuadrilateral}(2)
    cv = CellValues(qr, Lagrange{RefQuadrilateral, 1}())
    K = allocate_matrix(dh, ch)
    fvec = zeros(ndofs(dh))
    assembler = start_assemble(K, fvec)
    ke = zeros(4, 4); fe = zeros(4)
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        fill!(ke, 0.0); fill!(fe, 0.0)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, getcoordinates(cell))
            dΩ = getdetJdV(cv, qp)
            for i in 1:4
                fe[i] += f_rhs(x) * shape_value(cv, qp, i) * dΩ
                for j in 1:4
                    ke[i, j] += shape_gradient(cv, qp, j) ⋅ shape_gradient(cv, qp, i) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    apply!(K, fvec, ch)
    u = K \ fvec
    apply!(u, ch)
    return dh, u, cv
end

# a simple gradient indicator with top-fraction marking; see Ferrite's
# adaptive heat tutorial for a proper ZZ estimator with Dörfler marking
function mark_cells(dh, u, cv; frac = 0.25)
    η = zeros(Ferrite.getncells(Ferrite.get_grid(dh)))
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        ue = u[celldofs(cell)]
        h = norm(getcoordinates(cell)[3] - getcoordinates(cell)[1])
        s = 0.0
        for qp in 1:getnquadpoints(cv)
            s += norm(function_gradient(cv, qp, ue))^2 * getdetJdV(cv, qp)
        end
        η[cellid(cell)] = h * sqrt(s)
    end
    n = max(1, round(Int, frac * length(η)))
    return partialsortperm(η, 1:n; rev = true)
end
nothing # hide
```

Now the AMR loop. The dataset and the plots are created **once**, before the
loop; every sweep ends in a `regrid!`, and the figure at the end shows the
final state of the same plot objects that started on the 16-cell base mesh.

```@example amr
forest = Ferrite.ForestBWG(generate_grid(Quadrilateral, (4, 4)), 6)
grid = Ferrite.creategrid(forest)
dh, u, cv = solve_on(grid)
ds = FEData(dh, u)

fig = WGLMakie.Figure(size = (800, 420))
ax = WGLMakie.Axis(fig[1, 1]; aspect = WGLMakie.DataAspect(), title = "solution")
solutionplot!(ax, ds; color = :u, adaptive = true, solution_tol = 2.0e-3,
              max_depth = 6, colormap = :inferno, colorrange = (0.0, 1.0))
ax2 = WGLMakie.Axis(fig[1, 2]; aspect = WGLMakie.DataAspect(), title = "adapted mesh")
solutionplot!(ax2, ds; color = :u, adaptive = true, solution_tol = 2.0e-3,
              max_depth = 6, colormap = :inferno, colorrange = (0.0, 1.0))
meshplot!(ax2, ds; adaptive = true, geometry_tol = 1.0e-3, color = (:white, 0.6),
          plotnodes = false)

for sweep in 1:5
    marked = mark_cells(dh, u, cv)
    Ferrite.refine!(forest, marked)
    Ferrite.balanceforest!(forest)
    global grid = Ferrite.creategrid(forest)
    global dh, u, cv = solve_on(grid)
    FerriteViz.regrid!(ds, dh, u)          # the live plots follow
end
WGLMakie.Label(fig[0, :], "after 5 AMR sweeps: $(Ferrite.getncells(grid)) elements, " *
               "$(length(grid.conformity_info)) hanging nodes")
fig
```

## What survives a `regrid!`, and what does not

  * Adaptive `solutionplot`s and `meshplot`s of the **root** dataset follow,
    as above.
  * **Static** plots made before the regrid keep showing the old grid (they
    hold the old coordinate observables — stale, but alive); delete and
    recreate them.
  * **Derived** datasets (warps, gradients, quadrature-point data) capture
    the dof handler they were built on: regrid the root and re-apply the
    filter chain.

Between adaptations the dataset behaves like any other: `FerriteViz.update!`
with a matching dof vector, at the cost of a fixed-grid update. One full
adaptation step (forest mutation + `creategrid` + dof handler rebuild +
`regrid!` + first re-draw) costs a handful of plain updates and is paid once
per `refine!`, not per frame.
