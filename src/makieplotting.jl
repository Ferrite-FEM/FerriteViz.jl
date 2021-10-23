struct MakiePlotter{dim,DH<:Ferrite.AbstractDofHandler} <: AbstractPlotter
    dh::DH
    u::Vector #sorted
    cells_connectivity::Vector
    coords::Matrix{Float64}
    triangle_cells::Vector{NTuple{3,Int}}
end

function MakiePlotter(dh::Ferrite.AbstractDofHandler, u::Vector) 
    cells = Ferrite.getcells(dh.grid)
    dim = Ferrite.getdim(dh.grid)
    coords = [node.x[i] for node in Ferrite.getnodes(dh.grid), i in 1:Ferrite.getdim(dh.grid)] 
    connectivity = getproperty.(cells, :nodes) #TODO, needs interface function getnodes(::AbstractCell)
    triangle_cells = Vector{NTuple{3,Int}}()
    for cell in cells
        append!(triangle_cells, to_triangle(cell))
    end
    return MakiePlotter{dim,typeof(dh)}(dh,u,connectivity,coords,triangle_cells)
end

function reshape_triangles(plotter::MakiePlotter)
    reinterpret(reshape, Int, plotter.triangle_cells)'
end

function Makie.convert_arguments(P::Type{<:Makie.Mesh}, plotter::MakiePlotter)
    return Makie.convert_arguments(P,plotter.coords,reshape_triangles(plotter))
end

@recipe(SolutionPlot) do scene
    Attributes(
        scale_plot=false,
        shading=false,
        field_idx=1,
        process=postprocess,
        colormap=:viridis,
    )
end

function Makie.plot!(SP::SolutionPlot{<:Tuple{<:MakiePlotter}})
    plotter = SP[1][]
    solution = lift((x,y) -> reshape(dof_to_node(plotter.dh, plotter.u; field=x, process=y), Ferrite.getnnodes(plotter.dh.grid)),SP.attributes[:field_idx], SP.attributes[:process])
    return Makie.mesh!(SP,plotter, color=solution, shading=SP.attributes[:shading], scale_plot=SP.attributes[:scale_plot], colormap=SP.attributes[:colormap])
end

@recipe(Wireframe) do scene
    Attributes(
    plotnodes = true,
    strokewidth = 2,
    color=:black,
    markersize=30    
    )
end

function Makie.plot!(WF::Wireframe{<:Tuple{<:MakiePlotter{dim}}}) where dim
    plotter = WF[1][]
    coords = [node.x[i] for node in Ferrite.getnodes(plotter.dh.grid), i in 1:dim] 
    dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
    for cell in Ferrite.getcells(plotter.dh.grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [coords[e,:] for boundary in boundaryentities for e in boundary])
    end
    nodes = lift(x->x ? plotter.coords : zeros(Float32,0,3),WF.attributes[:plotnodes])
    Makie.scatter!(WF,nodes,markersize=WF.attributes[:markersize], color=WF.attributes[:color])
    return Makie.linesegments!(WF,lines,color=WF.attributes[:color], linewidth=WF.attributes[:strokewidth])
end

function Makie.plot!(WF::Wireframe{<:Tuple{<:Ferrite.AbstractGrid{dim}}}) where dim
    grid = WF[1][]
    coords = [Ferrite.getcoordinates(node)[i] for node in Ferrite.getnodes(grid), i in 1:dim] 
    dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
    for cell in Ferrite.getcells(grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [coords[e,:] for boundary in boundaryentities for e in boundary])
    end
    nodes = lift(x->x ? coords : zeros(Float32,0,3),WF.attributes[:plotnodes])
    Makie.scatter!(WF,nodes,markersize=WF.attributes[:markersize], color=WF.attributes[:color])
    return Makie.linesegments!(WF,lines,color=WF.attributes[:color], strokewidth=WF.attributes[:strokewidth])
end

@recipe(DeformedPlot) do scene
    Attributes(
    field_idx = 1,
    process = postprocess,
    scale_plot = false,
    shading = false,
    color = nothing,
    colormap = :viridis,
    scale = 1.0,
    )
end

function Makie.plot!(DP::DeformedPlot{<:Tuple{<:MakiePlotter}})
    plotter = DP[1][]
    @assert Ferrite.getfielddim(plotter.dh,DP.attributes[:field_idx][]) > 1
    nnodes = Ferrite.getnnodes(plotter.dh.grid)
    u_matrix = lift(x->dof_to_node(plotter.dh, plotter.u; field=x, process=identity), DP[:field_idx])
    solution = lift((x,y,z)-> z === nothing ? reshape(dof_to_node(plotter.dh, plotter.u; field=x, process=y), Ferrite.getnnodes(plotter.dh.grid)) : ones(length(solution[]))*x, DP.attributes[:field_idx], DP.attributes[:process], DP.attributes[:color])
    deformed_coords = lift((x,y) -> plotter.coords .+ (x .* y), DP.attributes[:scale], u_matrix)
    return Makie.mesh!(DP, deformed_coords, reshape_triangles(plotter), color=solution, scale_plot=DP.attributes[:scale_plot], shading=DP.attributes[:shading], colormap=DP.attributes[:colormap])
end

@recipe(Surface) do scene
    Attributes(
    field_idx = 1,
    process = postprocess,
    scale_plot = false,
    shading = false,
    colormap = :viridis,
    )
end
 
function Makie.plot!(SF::Surface{<:Tuple{<:MakiePlotter{2}}})
    plotter = SF[1][]
    solution = lift((x,y)->reshape(dof_to_node(plotter.dh, plotter.u; field=x, process=y), Ferrite.getnnodes(plotter.dh.grid)),SF.attributes[:field_idx],SF.attributes[:process])
    points = lift(x->[Makie.Point3f0(coord[1], coord[2], x[idx]) for (idx, coord) in enumerate(eachrow(plotter.coords))],solution)
    return Makie.mesh!(SF,points, reshape_triangles(plotter), color=solution, scale_plot=SF.attributes[:scale_plot], shading=SF.attributes[:shading], colormap=SF.attributes[:colormap])
end

function Makie.arrows(plotter::MakiePlotter, args...; field::Int=1, arrowsize=0.08, normalize=true, kwargs...) where T
    @assert Ferrite.getfielddim(plotter.dh,field) > 1
    solution = dof_to_node(plotter.dh, plotter.u; field=field, process=identity)
    if Ferrite.getdim(plotter.dh.grid) == 2
        Makie.arrows(plotter.coords[:,1], plotter.coords[:,2], solution[:,1], solution[:,2], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    elseif Ferrite.getdim(plotter.dh.grid) == 3
        Makie.arrows(plotter.coords[:,1], plotter.coords[:,2], plotter.coords[:,3], solution[:,1], solution[:,2], solution[:,3], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    end
end

function Makie.arrows!(plotter::MakiePlotter, args...; field::Int=1, arrowsize=0.08, normalize=false, kwargs...) where T
    @assert Ferrite.getfielddim(plotter.dh,field) > 1
    solution = dof_to_node(plotter.dh, plotter.u; field=field, process=identity)
    if Ferrite.getdim(plotter.dh.grid) == 2
        Makie.arrows!(plotter.coords[:,1], plotter.coords[:,2], solution[:,1], solution[:,2], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    elseif Ferrite.getdim(plotter.dh.grid) == 3
        Makie.arrows!(plotter.coords[:,1], plotter.coords[:,2], plotter.coords[:,3], solution[:,1], solution[:,2], solution[:,3], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    end
end

