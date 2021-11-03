# Advanced Topics

## Live plotting

Plotting while a computational heavy simulation is performed can be easily achieved with FerriteVis.jl.
Every plotter object of type `MakiePlotter` holds a property called `u` which is a so called `Observable`.
If an `Observable` changes, all its dependencies are triggered to change as well. So, all we need to do is to update
the observable `plotter.u`.
For this purpose the function [`FerriteVis.update!`](@ref) is provided. It takes a `plotter:MakiePlotter`
and a new solutiuon vector `u_new` and updates `plotter.u`, thereby all open plots called with `plotter` are updated.

A summary of the needed steps for live plotting:
1. Create a plotter before your time stepping begins
2. Call a plot or the `ferriteviewer` and save the return in a variable, e.g. `fig`
3. `display(fig)` in order to force the plot/viewer to pop up, even if its called inside a function body
4. `FerriteVis.update!(plotter,u_new)` where `u_new` corresponds to your new solution of the time step

As an illustrative example, let's consider the [plasticity example of Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/dev/examples/plasticity/).
In order to have live plotting we only need a few changes in the `solve` function.

```julia
using Ferrite, SparseArrays, LinearAlgebra, FerriteVis

struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}}
    G::T  # Shear modulus
    K::T  # Bulk modulus
    σ₀::T # Initial yield limit
    H::T  # Hardening modulus
    Dᵉ::S # Elastic stiffness tensor
end;

function J2Plasticity(E, ν, σ₀, H)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
    return J2Plasticity(G, K, σ₀, H, Dᵉ)
end;

mutable struct MaterialState{T, S <: SecondOrderTensor{3, T}}
    # Store "converged" values
    ϵᵖ::S # plastic strain
    σ::S # stress
    k::T # hardening variable

    # Store temporary values used during equilibrium iterations
    temp_ϵᵖ::S
    temp_σ::S
    temp_k::T
end

function MaterialState()
    return MaterialState(
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0,
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0)
end

function update_state!(state::MaterialState)
    state.ϵᵖ = state.temp_ϵᵖ
    state.σ = state.temp_σ
    state.k = state.temp_k
end;

function vonMises(σ)
    s = dev(σ)
    return sqrt(3.0/2.0 * s ⊡ s)
end;

function compute_stress_tangent(ϵ::SymmetricTensor{2, 3}, material::J2Plasticity, state::MaterialState)
    # unpack some material parameters
    G = material.G
    K = material.K
    H = material.H

    # We use (•)ᵗ to denote *trial*-values
    σᵗ = material.Dᵉ ⊡ (ϵ - state.ϵᵖ) # trial-stress
    sᵗ = dev(σᵗ)         # deviatoric part of trial-stress
    J₂ = 0.5 * sᵗ ⊡ sᵗ   # second invariant of sᵗ
    σᵗₑ = sqrt(3.0*J₂)   # effetive trial-stress (von Mises stress)
    σʸ = material.σ₀ + H * state.k # Previous yield limit

    φᵗ  = σᵗₑ - σʸ # Trial-value of the yield surface

    if φᵗ < 0.0 # elastic loading
        state.temp_σ = σᵗ
        return state.temp_σ, material.Dᵉ
    else # plastic loading
        h = H + 3G
        μ =  φᵗ / h   # plastic multiplier

        c1 = 1 - 3G * μ / σᵗₑ
        s = c1 * sᵗ           # updated deviatoric stress
        σ = s + vol(σᵗ)       # updated stress

        # Compute algorithmic tangent stiffness ``D = \frac{\Delta \sigma }{\Delta \epsilon}``
        κ = H * (state.k + μ) # drag stress
        σₑ = material.σ₀ + κ  # updated yield surface

        δ(i,j) = i == j ? 1.0 : 0.0
        Isymdev(i,j,k,l)  = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
        Q(i,j,k,l) = Isymdev(i,j,k,l) - 3.0 / (2.0*σₑ^2) * s[i,j]*s[k,l]
        b = (3G*μ/σₑ) / (1.0 + 3G*μ/σₑ)

        Dtemp(i,j,k,l) = -2G*b * Q(i,j,k,l) - 9G^2 / (h*σₑ^2) * s[i,j]*s[k,l]
        D = material.Dᵉ + SymmetricTensor{4, 3}(Dtemp)

        # Store outputs in the material state
        Δϵᵖ = 3/2 *μ / σₑ*s            # plastic strain
        state.temp_ϵᵖ = state.ϵᵖ + Δϵᵖ  # plastic strain
        state.temp_k = state.k + μ     # hardening variable
        state.temp_σ = σ               # updated stress
        return state.temp_σ, D
    end
end

function create_values(interpolation)
    # setup quadrature rules
    qr      = QuadratureRule{3,RefTetrahedron}(2)
    face_qr = QuadratureRule{2,RefTetrahedron}(3)

    # create geometric interpolation (use the same as for u)
    interpolation_geom = Lagrange{3,RefTetrahedron,1}()

    # cell and facevalues for u
    cellvalues_u = CellVectorValues(qr, interpolation, interpolation_geom)
    facevalues_u = FaceVectorValues(face_qr, interpolation, interpolation_geom)

    return cellvalues_u, facevalues_u
end;

function create_dofhandler(grid, interpolation)
    dh = DofHandler(grid)
    dim = 3
    push!(dh, :u, dim, interpolation) # add a displacement field with 3 components
    close!(dh)
    return dh
end

function create_bc(dh, grid)
    dbcs = ConstraintHandler(dh)
    # Clamped on the left side
    dofs = [1, 2, 3]
    dbc = Dirichlet(:u, getfaceset(grid, "left"), (x,t) -> [0.0, 0.0, 0.0], dofs)
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end;

function doassemble(cellvalues::CellVectorValues{dim},
                    facevalues::FaceVectorValues{dim}, K::SparseMatrixCSC, grid::Grid,
                    dh::DofHandler, material::J2Plasticity, u, states, t) where {dim}
    r = zeros(ndofs(dh))
    assembler = start_assemble(K, r)
    nu = getnbasefunctions(cellvalues)
    re = zeros(nu)     # element residual vector
    ke = zeros(nu, nu) # element tangent matrix

    for (cell, state) in zip(CellIterator(dh), states)
        fill!(ke, 0)
        fill!(re, 0)
        eldofs = celldofs(cell)
        ue = u[eldofs]
        assemble_cell!(ke, re, cell, cellvalues, facevalues, grid, material,
                       ue, state, t)
        assemble!(assembler, eldofs, re, ke)
    end
    return K, r
end

function assemble_cell!(Ke, re, cell, cellvalues, facevalues, grid, material,
                        ue, state, t)
    n_basefuncs = getnbasefunctions(cellvalues)
    reinit!(cellvalues, cell)

    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
        ∇u = function_gradient(cellvalues, q_point, ue)
        ϵ = symmetric(∇u) # Total strain
        σ, D = compute_stress_tangent(ϵ, material, state[q_point])

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = symmetric(shape_gradient(cellvalues, q_point, i))

            re[i] += (δϵ ⊡ σ) * dΩ # add internal force to residual
            for j in 1:i
                Δϵ = symmetric(shape_gradient(cellvalues, q_point, j))
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)

    # Add traction as a negative contribution to the element residual `re`:
    for face in 1:nfaces(cell)
        if onboundary(cell, face) && (cellid(cell), face) ∈ getfaceset(grid, "right")
            reinit!(facevalues, cell, face)
            for q_point in 1:getnquadpoints(facevalues)
                dΓ = getdetJdV(facevalues, q_point)
                for i in 1:n_basefuncs
                    δu = shape_value(facevalues, q_point, i)
                    re[i] -= (δu ⋅ t) * dΓ
                end
            end
        end
    end

end

function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

function solve()
    # Define material parameters
    E = 200.0e9 # [Pa]
    H = E/20   # [Pa]
    ν = 0.3     # [-]
    σ₀ = 200e6  # [Pa]
    material = J2Plasticity(E, ν, σ₀, H)

    L = 10.0 # beam length [m]
    w = 1.0  # beam width [m]
    h = 1.0  # beam height[m]
    n_timesteps = 100
    u_max = zeros(n_timesteps)
    traction_magnitude = 1.e7 * range(0.5, 1.0, length=n_timesteps)

    # Create geometry, dofs and boundary conditions
    n = 2
    nels = (10n, n, 2n) # number of elements in each spatial direction
    P1 = Vec((0.0, 0.0, 0.0))  # start point for geometry
    P2 = Vec((L, w, h))        # end point for geometry
    grid = generate_grid(Tetrahedron, nels, P1, P2)
    interpolation = Lagrange{3, RefTetrahedron, 1}() # Linear tet with 3 unknowns/node

    dh = create_dofhandler(grid, interpolation) # JuaFEM helper function
    dbcs = create_bc(dh, grid) # create Dirichlet boundary-conditions

    cellvalues, facevalues = create_values(interpolation)

    # Pre-allocate solution vectors, etc.
    n_dofs = ndofs(dh)  # total number of dofs
    u  = zeros(n_dofs)  # solution vector
    plotter = MakiePlotter(dh,u)
    fig = ferriteviewer(plotter)
    display(fig)
    Δu = zeros(n_dofs)  # displacement correction
    r = zeros(n_dofs)   # residual
    K = create_sparsity_pattern(dh); # tangent stiffness matrix

    nqp = getnquadpoints(cellvalues)
    states = [[MaterialState() for _ in 1:nqp] for _ in 1:getncells(grid)]

    # Newton-Raphson loop
    NEWTON_TOL = 1 # 1 N

    for timestep in 1:n_timesteps
        t = timestep # actual time (used for evaluating d-bndc)
        traction = Vec((0.0, 0.0, traction_magnitude[timestep]))
        newton_itr = -1
        update!(dbcs, t) # evaluates the D-bndc at time t
        apply!(u, dbcs)  # set the prescribed values in the solution vector

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
        FerriteVis.update!(plotter,u)
        sleep(0.1)

        # Update all the material states after we have reached equilibrium
        for cell_states in states
            foreach(update_state!, cell_states)
        end
        u_max[timestep] = max(abs.(u)...) # maximum displacement in current timestep
    end

    # ## Postprocessing
    # Only a vtu-file corrsponding to the last time-step is exported.
    #
    # The following is a quick (and dirty) way of extracting average cell data for export.
    mises_values = zeros(getncells(grid))
    κ_values = zeros(getncells(grid))
    for (el, cell_states) in enumerate(states)
        for state in cell_states
            mises_values[el] += vonMises(state.σ)
            κ_values[el] += state.k*material.H
        end
        mises_values[el] /= length(cell_states) # average von Mises stress
        κ_values[el] /= length(cell_states)     # average drag stress
    end
    return u, dh, traction_magnitude
end

u, dh, traction_magnitude = solve();
```

Note that we create `plotter::MakiePlotter` object before the time stepping begins, as well as calling `ferriteviewer` on the `plotter`.
The next function call is crucial to get the live plotting working. `display(fig)` forces the viewer to pop up, even if it's inside a function body.
Now, the only missing piece is the `FerriteVis.update!` of the plotter, which happens directly after the Newton iteration. The result for this code looks like this:

![liveplot](assets/liveplotting.gif)

Since the computational load of one time step is in this example too low, the plotter would just update all the time and likely never display something, so we artificially increate the load of one time step by
`sleep`ing for 0.1s.

If you don't need the full viewer as a live plot, you can of course call instead `solutionplot` (or any other plot/plot combination) with appropriate keyword arguments to only have a specific live plot.
This can be beneficial performancewise.
