# The 3D adaptive heat example (Ferrite's heat_adaptivity tutorial, re-centered
# so the feature is visible from outside), rendered live through FerriteViz:
# one dataset, one solutionplot + adaptive meshplot, following the ForestBWG
# through the AMR sweeps via regrid!.
using GLMakie, Ferrite, FerriteViz, SparseArrays, LinearAlgebra
GLMakie.activate!()
outdir = @__DIR__

# Manufactured Gaussian ring on the sphere ‖x - c‖ = 0.7 around a cube corner,
# so the feature cuts the three adjacent faces.
const C0 = Ferrite.Vec(1.0, 1.0, 1.0)
u_exact(x) = exp(-((norm(x - C0) - 0.7) / 0.12)^2)
f_rhs(x) = -laplace(u_exact, x)   # Tensors.laplace, reexported by Ferrite

function solve_on(grid)
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefHexahedron,1}())
    close!(dh)
    ch = ConstraintHandler(dh)
    add!(ch, Ferrite.ConformityConstraint(:u))            # hanging nodes
    add!(ch, Dirichlet(:u, union((getfacetset(grid, f) for f in
        ("left", "right", "front", "back", "bottom", "top"))...), x -> u_exact(x)))
    close!(ch)
    qr = QuadratureRule{RefHexahedron}(2)
    cv = CellValues(qr, Lagrange{RefHexahedron,1}())
    K = allocate_matrix(dh, ch)
    fvec = zeros(ndofs(dh))
    assembler = start_assemble(K, fvec)
    ke = zeros(8, 8); fe = zeros(8)
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        fill!(ke, 0.0); fill!(fe, 0.0)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, getcoordinates(cell))
            dΩ = getdetJdV(cv, qp)
            for i in 1:8
                fe[i] += f_rhs(x) * shape_value(cv, qp, i) * dΩ
                for j in 1:8
                    ke[i, j] += shape_gradient(cv, qp, j) ⋅ shape_gradient(cv, qp, i) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    apply!(K, fvec, ch)
    u = K \ fvec
    apply!(u, ch)
    (dh, u, cv)
end

# Simple gradient indicator + top-fraction marking (see Ferrite's
# heat_adaptivity tutorial for a proper ZZ estimator with Dörfler marking —
# the visualization does not care how the marks were produced).
function mark_cells(dh, u, cv; frac = 0.2)
    η = zeros(Ferrite.getncells(Ferrite.get_grid(dh)))
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        ue = u[celldofs(cell)]
        h = norm(getcoordinates(cell)[7] - getcoordinates(cell)[1])
        s = 0.0
        for qp in 1:getnquadpoints(cv)
            s += norm(function_gradient(cv, qp, ue))^2 * getdetJdV(cv, qp)
        end
        η[cellid(cell)] = h * sqrt(s)
    end
    n = max(1, round(Int, frac * length(η)))
    return partialsortperm(η, 1:n; rev = true)
end

forest = Ferrite.ForestBWG(generate_grid(Hexahedron, (4, 4, 4)), 6)
grid = Ferrite.creategrid(forest)
dh, u, cv = solve_on(grid)
ds = FEData(dh, u)

fig = Figure(size = (900, 780))
title = Observable("")
ax = Axis3(fig[1, 1]; title, aspect = :data, azimuth = 0.12π, elevation = 0.3π,
           protrusions = 0)
hidedecorations!(ax); hidespines!(ax)
solutionplot!(ax, ds; color = :u, adaptive = true, solution_tol = 2.0e-3, max_depth = 5,
              colormap = :inferno, colorrange = (0.0, 1.0))
meshplot!(ax, ds; adaptive = true, geometry_tol = 1.0e-3, color = (:white, 0.35),
          plotnodes = false)

record(fig, joinpath(outdir, "amr_heat.mp4"); framerate = 1) do io
    for sweep in 0:6
        title[] = "AMR sweep $sweep — $(Ferrite.getncells(grid)) elements, " *
                  "$(length(grid.conformity_info)) hanging nodes"
        recordframe!(io)
        sweep == 6 && break
        marked = mark_cells(dh, u, cv)
        Ferrite.refine!(forest, marked)
        Ferrite.balanceforest!(forest)
        global grid = Ferrite.creategrid(forest)
        global dh, u, cv = solve_on(grid)
        FerriteViz.regrid!(ds, dh, u)     # the live plots follow
        ax.azimuth[] += 0.015π
    end
end
println("final: ", Ferrite.getncells(grid), " elements, ",
        length(grid.conformity_info), " hanging nodes")
