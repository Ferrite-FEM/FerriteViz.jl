function Makie.mesh(plotter::MakiePlotter, args...; field_idx::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(transfer_solution(plotter; field_idx=field_idx, process=process), num_vertices(plotter))
    return Makie.mesh(plotter.physical_coords,plotter.triangles, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.mesh!(plotter::MakiePlotter, args...; field_idx::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(transfer_solution(plotter; field_idx=field_idx, process=process), num_vertices(plotter))
    return Makie.mesh!(plotter.physical_coords,plotter.triangles, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.surface(plotter::MakiePlotter{2}, args...; field_idx::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(transfer_solution(plotter; field_idx=field_idx, process=process), num_vertices(plotter))
    points = [Makie.Point3f0(coord[1], coord[2], solution[idx]) for (idx, coord) in enumerate(eachrow(plotter.physical_coords))]
    return Makie.mesh(points,plotter.triangles, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.surface!(plotter::MakiePlotter{2}, args...; field_idx::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(transfer_solution(plotter; field_idx=field_idx, process=process), num_vertices(plotter))
    points = [Makie.Point3f0(coord[1], coord[2], solution[idx]) for (idx, coord) in enumerate(eachrow(plotter.physical_coords))]
    return Makie.mesh!(points,plotter.triangles, color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.arrows(plotter::MakiePlotter, args...; field_idx::Int=1, arrowsize=0.08, normalize=true, kwargs...) where T
    @assert Ferrite.getfielddim(plotter.dh,field_idx) > 1
    solution = transfer_solution(plotter; field_idx=field_idx, process=identity)
    if Ferrite.getdim(plotter.dh.grid) == 2
        Makie.arrows(plotter.physical_coords[:,1], plotter.physical_coords[:,2], solution[:,1], solution[:,2], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    elseif Ferrite.getdim(plotter.dh.grid) == 3
        Makie.arrows(plotter.physical_coords[:,1], plotter.physical_coords[:,2], plotter.physical_coords[:,3], solution[:,1], solution[:,2], solution[:,3], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    end
end

function Makie.arrows!(plotter::MakiePlotter, args...; field_idx::Int=1, arrowsize=0.08, normalize=false, kwargs...) where T
    @assert Ferrite.getfielddim(plotter.dh,field_idx) > 1
    solution = transfer_solution(plotter; field_idx=field_idx, process=identity)
    if Ferrite.getdim(plotter.dh.grid) == 2
        Makie.arrows!(plotter.physical_coords[:,1], plotter.physical_coords[:,2], solution[:,1], solution[:,2], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    elseif Ferrite.getdim(plotter.dh.grid) == 3
        Makie.arrows!(plotter.physical_coords[:,1], plotter.physical_coords[:,2], plotter.physical_coords[:,3], solution[:,1], solution[:,2], solution[:,3], args...; arrowsize=arrowsize, normalize=normalize, kwargs...)
    end
end

function warp_by_vector(plotter::MakiePlotter, args...; field_idx::Int=1, scale=1.0, process::Function=postprocess, scale_plot=false, shading=false, color=nothing, kwargs...) where T
    @assert Ferrite.getfielddim(plotter.dh,field_idx) > 1
    u_matrix = transfer_solution(plotter; field_idx=field_idx, process=identity)
    solution = reshape(transfer_solution(plotter; field_idx=field_idx, process=process), num_vertices(plotter))
    plot_color = color===nothing ? solution : color
    return Makie.mesh(plotter.physical_coords .+ (scale .* u_matrix), color=plot_color,plotter.triangles, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function warp_by_vector!(plotter::MakiePlotter, args...; field_idx::Int=1, scale=1.0, process::Function=postprocess, scale_plot=false, shading=false, color=nothing, kwargs...) where T
    @assert Ferrite.getfielddim(plotter.dh,field_idx) > 1
    u_matrix = transfer_solution(plotter; field_idx=field_idx, process=identity)
    solution = reshape(transfer_solution(plotter; field_idx=field_idx, process=process), num_vertices(plotter))
    plot_color = color===nothing ? solution : color
    return Makie.mesh!(plotter.physical_coords .+ (scale .* u_matrix), color=plot_color,plotter.triangles, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function plot_grid(plotter::MakiePlotter{dim},args...;color=:black, strokewidth=3, plotnodes=true, markersize=(dim == 2 ? 20 : 70), kwargs...) where dim
    dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
    cells = Ferrite.getcells(plotter.dh.grid)
    for cell in Ferrite.getcells(plotter.dh.grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [plotter.physical_coords[e,:] for boundary in boundaryentities for e in boundary])
    end
    if plotnodes
        Makie.scatter(plotter.physical_coords, args...; markersize=markersize, color=color, kwargs...)
        Makie.linesegments!(lines,args...;color=color, strokewidth=strokewidth, kwargs...)
        return Makie.current_figure()
    end
    return Makie.linesegments(lines,args...;color=color, strokewidth=strokewidth, kwargs...)
end

function plot_grid!(plotter::MakiePlotter{dim},args...;color=:black, strokewidth=3, plotnodes=true, markersize=(dim == 2 ? 20 : 70), kwargs...) where dim
    dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
    cells = Ferrite.getcells(plotter.dh.grid)
    for cell in Ferrite.getcells(plotter.dh.grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [plotter.physical_coords[e,:] for boundary in boundaryentities for e in boundary])
    end
    if plotnodes
        Makie.scatter!(plotter.physical_coords, args...; markersize=markersize, color=color, kwargs...)
    end
    return Makie.linesegments!(lines,args...;color=color, strokewidth=strokewidth, kwargs...)
end

function plot_grid(grid::Ferrite.AbstractGrid{dim},args...;color=:black, strokewidth=3, plotnodes=true, markersize=(dim == 2 ? 20 : 70), kwargs...) where dim
    dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
    coords = [node.x[i] for node in Ferrite.getnodes(grid), i in 1:Ferrite.getdim(grid)]
    cells = Ferrite.getcells(grid)
    for cell in Ferrite.getcells(grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [coords[e,:] for boundary in boundaryentities for e in boundary])
    end
    if plotnodes
        Makie.scatter(coords, args...; markersize=markersize, color=color, kwargs...)
        Makie.linesegments!(lines,args...;color=color, strokewidth=strokewidth, kwargs...)
        return Makie.current_figure()
    end
    return Makie.linesegments(lines,args...;color=color, strokewidth=strokewidth, kwargs...)
end

function plot_grid!(grid::Ferrite.AbstractGrid{dim},args...;color=:black, strokewidth=3, plotnodes=true, markersize=(dim == 2 ? 20 : 70), kwargs...) where dim
    dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
    coords = [node.x[i] for node in Ferrite.getnodes(grid), i in 1:Ferrite.getdim(grid)]
    cells = Ferrite.getcells(grid)
    for cell in Ferrite.getcells(grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [coords[e,:] for boundary in boundaryentities for e in boundary])
    end
    if plotnodes
        Makie.scatter!(coords, args...; markersize=markersize, color=color, kwargs...)
    end
    return Makie.linesegments!(lines,args...;color=color, strokewidth=strokewidth, kwargs...)
end
