struct MakiePlotter{dim} <: AbstractPlotter
    dh::Ferrite.DofHandler{dim}
    u::Vector #sorted
    cells_connectivity::Vector
    coords::Matrix{Float64}
    triangle_cells::Vector{NTuple{3,Int}}
end

function MakiePlotter(dh::Ferrite.AbstractDofHandler, u::Vector) 
    cells = Ferrite.getcells(dh.grid)
    coords = [node.x[i] for node in Ferrite.getnodes(dh.grid), i in 1:Ferrite.getdim(dh.grid)] 
    connectivity = getproperty.(cells, :nodes) #TODO, needs interface function getnodes(::AbstractCell)
    triangle_cells = Vector{NTuple{3,Int}}()
    for cell in cells
        append!(triangle_cells, to_triangle(cell))
    end
    return MakiePlotter(dh,u,connectivity,coords,triangle_cells)
end

function reshape_triangles(plotter::MakiePlotter)
    reinterpret(reshape, Int, plotter.triangle_cells)'
end

function Makie.mesh(plotter::MakiePlotter, args...; field::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(dof_to_node(plotter.dh, plotter.u; field=field, process=process), Ferrite.getnnodes(plotter.dh.grid))
    return Makie.mesh(plotter.coords, reshape_triangles(plotter), color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.mesh!(plotter::MakiePlotter, args...; field::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(dof_to_node(plotter.dh, plotter.u; field=field, process=process), Ferrite.getnnodes(plotter.dh.grid))
    return Makie.mesh!(plotter.coords, reshape_triangles(plotter), color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.surface(plotter::MakiePlotter{2}, args...; field::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(dof_to_node(plotter.dh, plotter.u; field=field, process=process), Ferrite.getnnodes(plotter.dh.grid))
    points = [Makie.Point3f0(coord[1], coord[2], solution[idx]) for (idx, coord) in enumerate(eachrow(plotter.coords))]
    return Makie.mesh(points, reshape_triangles(plotter), color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.surface!(plotter::MakiePlotter{2}, args...; field::Int=1, process::Function=postprocess, scale_plot=false, shading=false, kwargs...) where T
    solution = reshape(dof_to_node(plotter.dh, plotter.u; field=field, process=process), Ferrite.getnnodes(plotter.dh.grid))
    points = [Makie.Point3f0(coord[1], coord[2], solution[idx]) for (idx, coord) in enumerate(eachrow(plotter.coords))]
    return Makie.mesh!(points, reshape_triangles(plotter), color=solution, args...; scale_plot=scale_plot, shading=shading, kwargs...)
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

function warp_by_vector(plotter::MakiePlotter, args...; field::Int=1, scale=1.0, process::Function=postprocess, scale_plot=false, shading=false, color=nothing, kwargs...) where T
    @assert Ferrite.getfielddim(plotter.dh,field) > 1
    u_matrix = dof_to_node(plotter.dh, plotter.u; field=field, process=identity)
    solution = reshape(dof_to_node(plotter.dh, plotter.u; field=field, process=process), Ferrite.getnnodes(plotter.dh.grid))
    plot_color = color===nothing ? solution : color
    return Makie.mesh(plotter.coords .+ (scale .* u_matrix), color=plot_color, reshape_triangles(plotter), args...; scale_plot=scale_plot, shading=shading, kwargs...)
end

function Makie.poly(plotter::MakiePlotter{2}, args...; color=:transparent, strokecolor=:black, strokewidth=3, kwargs...) where dim
    p = Makie.Polygon[] 
    for cell in plotter.cells_connectivity
        points = Makie.Point2f0[]
        for node in cell
            push!(points, Makie.Point2f0(plotter.coords[node,:]))
        end
        push!(p, Makie.Polygon(points))
    end
    return Makie.poly(p,args...;color=color,strokecolor=strokecolor, strokewidth=strokewidth, kwargs...)
end

function Makie.poly!(plotter::MakiePlotter{2}, args...; color=:transparent, strokecolor=:black, strokewidth=3, kwargs...) where dim
    p = Makie.Polygon[] 
    for cell in plotter.cells_connectivity
        points = Makie.Point2f0[]
        for node in cell
            push!(points, Makie.Point2f0(plotter.coords[node,:]))
        end
        push!(p, Makie.Polygon(points))
    end
    return Makie.poly(p,args...;color=color,strokecolor=strokecolor, strokewidth=strokewidth, kwargs...)
end
