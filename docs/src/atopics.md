# Advanced Topics

```@example 1
import JSServe # hide
JSServe.Page(exportable=true, offline=true) # hide
```

## Gradient field visualization

FerriteViz also makes it easy to visualize gradient fields, like for example strain or stress fields.
A common approach to visualize stresses and strains is to compute the L2 projection onto a H1 field and plot this.
However, a big downside is that we loose the ability to investigate the jumps between elements, as they get smoothed out, hiding possible issues in the solution.
Therefore, we provide the ability to interpolate the gradient into a piecewise discontinuous field via `FerriteViz.interpolate_gradient_field`.
This function may be moved to Ferrite in the future.

In this quick example we show how to visualize strains and stresses side-by-side
```@example 1
using Ferrite
import FerriteViz
using FerriteViz: ε
import WGLMakie #activating the backend, switch to GLMakie or CairoMakie (for 2D) locally

include("ferrite-examples/incompressible-elasticity.jl") #defines dh_linear, dh_quadratic, u_linear, u_quadratic and mp

(dh_linear_grad, u_linear_grad) = FerriteViz.interpolate_gradient_field(dh_linear, u_linear, :u)
(dh_quadratic_grad, u_quadratic_grad) = FerriteViz.interpolate_gradient_field(dh_quadratic, u_quadratic, :u)
plotter_linear = FerriteViz.MakiePlotter(dh_linear_grad, u_linear_grad)
plotter_quadratic = FerriteViz.MakiePlotter(dh_quadratic_grad, u_quadratic_grad)
σ(∇u) = 2*mp.G*dev(ε(∇u)) + mp.K*tr(ε(∇u))*ones(ε(∇u)) #helper function to map gradient to stress
cmap = :jet

f = WGLMakie.Figure()
axs = [WGLMakie.Axis(f[1, 1], title="Strain norm (linear)"),WGLMakie.Axis(f[1, 2], title="Stress norm (linear)"),WGLMakie.Axis(f[1, 3], title="Pressure (deformed, linear)"),
       WGLMakie.Axis(f[3, 1], title="Strain norm (quadratic)"),WGLMakie.Axis(f[3, 2], title="Stress norm (quadratic)"),WGLMakie.Axis(f[3, 3], title="Pressure (deformed, quadratic)")]
p1 = FerriteViz.solutionplot!(axs[1], plotter_linear, process=∇u->norm(ε(∇u)), colormap=cmap, field=:gradient)
p2 = FerriteViz.solutionplot!(axs[2], plotter_linear, process=∇u->norm(σ(∇u)), colormap=cmap, field=:gradient)
p3 = FerriteViz.solutionplot!(axs[3], dh_linear, u_linear, field=:p, deformation_field=:u, colormap=cmap)
f[2,1] = WGLMakie.Colorbar(f[1,1], p1, vertical=false)
f[2,2] = WGLMakie.Colorbar(f[1,2], p2, vertical=false)
f[2,3] = WGLMakie.Colorbar(f[1,3], p3, vertical=false)

p4 = FerriteViz.solutionplot!(axs[4], plotter_quadratic, process=∇u->norm(ε(∇u)), colormap=cmap, field=:gradient)
p5 = FerriteViz.solutionplot!(axs[5], plotter_quadratic, process=∇u->norm(σ(∇u)), colormap=cmap, field=:gradient)
p6 = FerriteViz.solutionplot!(axs[6], dh_quadratic, u_quadratic, field=:p, deformation_field=:u, colormap=cmap)
f[4,1] = WGLMakie.Colorbar(f[3,1], p1, vertical=false)
f[4,2] = WGLMakie.Colorbar(f[3,2], p2, vertical=false)
f[4,3] = WGLMakie.Colorbar(f[3,3], p3, vertical=false)

f
```

An alternative to this approach is to compute gradient quantities at samples points and plot these via `arrows`.

## High-order fields

The investigation of high-order fields is currently only supported via a first-order refinment of the problem.
Here, the high-order approximation is replaced by a first order approximation of the field, which is
spanned by the nodes of the high-order approximation. For example, the first order refinement of a
heat problem on a square domain for Lagrange polynomials of order 5 looks like this:
```@example 1
include("ferrite-examples/heat-equation.jl"); #defines manufactured_heat_problem

f = WGLMakie.Figure()
axs = [WGLMakie.Axis3(f[1, 1], title="Coarse"), WGLMakie.Axis3(f[1, 2], title="Fine")]

dh,u = manufactured_heat_problem(Triangle, Lagrange{2,RefTetrahedron,5}(), 1)
dh_for,u_for = FerriteViz.for_discretization(dh, u)
plotter_for = FerriteViz.MakiePlotter(dh_for, u_for)
FerriteViz.surface!(axs[1], plotter_for)

dh,u = manufactured_heat_problem(Triangle, Lagrange{2,RefTetrahedron,5}(), 3)
dh_for,u_for = FerriteViz.for_discretization(dh, u)
plotter_for = FerriteViz.MakiePlotter(dh_for, u_for)
FerriteViz.surface!(axs[2], plotter_for)

f
```
Note that this method produces small artifacts due to the flattening of the nonlinearities of the high order ansatz.
However, it is still sufficient to investigate important features of the solution. 
In future we will also provide an adaptive tessellation algorithm to resolve the high-order fields with full detail.

## Live plotting

Plotting while a computational heavy simulation is performed can be easily achieved with FerriteViz.jl.
Every plotter object of type `MakiePlotter` holds a property called `u` which is a so called `Observable`.
If an `Observable` changes, all its dependencies are triggered to change as well. So, all we need to do is to update
the observable `plotter.u`.
For this purpose the function [`FerriteViz.update!`](@ref) is provided. It takes a `plotter:MakiePlotter`
and a new solutiuon vector `u_new` and updates `plotter.u`, thereby all open plots called with `plotter` are updated.

A summary of the needed steps for live plotting:
1. Create a plotter before your time stepping begins
2. Call a plot or the `ferriteviewer` and save the return in a variable, e.g. `fig`
3. `display(fig)` in order to force the plot/viewer to pop up, even if its called inside a function body
4. `FerriteViz.update!(plotter,u_new)` where `u_new` corresponds to your new solution of the time step

As an illustrative example, let's consider a slightly modified [plasticity example of Ferrite.jl](https://github.com/koehlerson/FerriteViz.jl/blob/master/docs/src/ferrite-examples/plasticity.jl).
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
        plotter = MakiePlotter(dh,u)
        fig = ferriteviewer(plotter)
        display(fig)
        ####################################################################
    end
    
    Δu = zeros(n_dofs)  # displacement correction
    r = zeros(n_dofs)   # residual
    K = create_sparsity_pattern(dh); # tangent stiffness matrix

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
            ####### Step 4 updating the current solution vector in plotter ####### 
            FerriteViz.update!(plotter,u)
            ###################################################################### 
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

Note that we create `plotter::MakiePlotter` object before the time stepping begins, as well as calling `ferriteviewer` on the `plotter`.
The next function call is crucial to get the live plotting working. `display(fig)` forces the viewer to pop up, even if it's inside a function body.
Now, the only missing piece is the `FerriteViz.update!` of the plotter, which happens directly after the Newton iteration. The result for this code looks like this:

![liveplot](https://github.com/Ferrite-FEM/FerriteViz.jl/blob/master/docs/src/assets/liveplotting.gif?raw=true)

Since the computational load of one time step is in this example too low, the plotter would just update all the time and likely never display something, so we artificially increase the load of one time step by
`sleep`ing for 0.1s.

If you don't need the full viewer as a live plot, you can of course call instead `solutionplot` (or any other plot/plot combination) with appropriate keyword arguments to only have a specific live plot.
This can be beneficial performancewise.
