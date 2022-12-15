# Advanced Topics

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

## Gradient field visualization

FerriteViz also makes it easy to visualize gradient fields, like for example strain or stress fields.
A common approach to visualize stresses and strains is to compute the L2 projection onto a H1 field and plot this.
However, a big downside is that we loose the ability to investigate the jumps between elements, as they get smoothed out, hiding possible issues in the solution.
Therefore, we provide the ability to interpolate the gradient into a piecewise discontinuous field via `FerriteViz.interpolate_gradient_field`.
This function may be moved to Ferrite in the future.


```@example 1
using FerriteViz: ε
include("ferrite-examples/incompressible-elasticity.jl") #only defines solving function
plotter = FerriteViz.MakiePlotter(dh,u)
(dh_grad, u_grad) = FerriteViz.interpolate_gradient_field(dh, u, :u)
FerriteViz.solutionplot(dh_grad, u_grad, process=x->norm(ε(x)))
WGLMakie.current_figure()
```

An alternative to this approach is to compute gradient quantities at samples points and plot these via `arrows`.
