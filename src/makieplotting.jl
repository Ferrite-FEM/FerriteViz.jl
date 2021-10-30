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

const FerriteVisPlots = Union{Type{<:Wireframe},Type{<:SolutionPlot},Type{<:Arrows}}

function Makie.convert_arguments(P::FerriteVisPlots, dh::Ferrite.AbstractDofHandler, u::Vector)
    return (MakiePlotter(dh,u),)
end

@recipe(SolutionPlot) do scene
    Attributes(
    scale_plot=false,
    shading=false,
    field=:default,
    deformation_field=:default,
    process=postprocess,
    colormap=:viridis,
    colorrange=(0,1),
    transparent=false,
    deformation_scale = 1.0,
    )
end

function Makie.plot!(SP::SolutionPlot{<:Tuple{<:MakiePlotter}})
    plotter = SP[1][]
    solution = lift((x,y) -> x===:default ? reshape(dof_to_node(plotter.dh, plotter.u; field=1, process=y), Ferrite.getnnodes(plotter.dh.grid)) : reshape(dof_to_node(plotter.dh, plotter.u; field=Ferrite.find_field(plotter.dh,x), process=y), Ferrite.getnnodes(plotter.dh.grid)),SP[:field], SP[:process])
    u_matrix = lift(x->x===:default ? zeros(0,3) : dof_to_node(plotter.dh, plotter.u; field=Ferrite.find_field(plotter.dh,x), process=identity), SP[:deformation_field])
    coords = lift((x,y,z) -> z===:default ? plotter.coords : plotter.coords .+ (x .* y) , SP[:deformation_scale], u_matrix, SP[:deformation_field])
    mins = @lift(minimum($solution))
    maxs = @lift(maximum($solution))
    SP[:colorrange] = @lift(($mins,$maxs))
    return Makie.mesh!(SP, coords, reshape_triangles(plotter), color=solution, shading=SP[:shading], scale_plot=SP[:scale_plot], colormap=SP[:colormap], transparent=SP[:transparent])
end

@recipe(Wireframe) do scene
    Attributes(
    plotnodes=true,
    strokewidth=2,
    color=theme(scene, :linecolor),
    markersize=30,
    deformation_field=:default,
    visible=true,
    scale=1,
    textsize=15,
    offset=(0.0,0.0),
    nodelabels=false,
    nodelabelcolor=:darkblue,
    celllabels=false,
    celllabelcolor=:darkred
    )
end

function Makie.plot!(WF::Wireframe{<:Tuple{<:MakiePlotter{dim}}}) where dim
    plotter = WF[1][]
    physical_coords = [node.x[i] for node in Ferrite.getnodes(plotter.dh.grid), i in 1:dim] 
    u_matrix = lift(x->x===:default ? zeros(0,3) : dof_to_node(plotter.dh, plotter.u; field=Ferrite.find_field(plotter.dh,x), process=identity), WF[:deformation_field])
    coords = lift((x,y,z) -> z===:default ? physical_coords : physical_coords .+ (x .* y) , WF[:scale], u_matrix, WF[:deformation_field])
    lines = @lift begin
        dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
        for cell in Ferrite.getcells(plotter.dh.grid)
            boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
            append!(lines, [$coords[e,:] for boundary in boundaryentities for e in boundary])
        end
        lines
    end
    nodes = lift((x,y)->x ? y : zeros(Float32,0,3),WF[:plotnodes], coords)
    #plot the nodes
    Makie.scatter!(WF,coords,markersize=WF[:markersize], color=WF[:color], visible=WF[:visible])
    #set up nodelabels
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size($coords,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? [dim < 3 ? Point2f0(row) : Point3f0(row) for row in eachrow($coords)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
    #set up celllabels
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:Ferrite.getncells(plotter.dh.grid)] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,$coords) for cell in Ferrite.getcells(plotter.dh.grid)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
    Makie.text!(WF,nodelabels, position=nodepositions, textsize=WF[:textsize], offset=WF[:offset],color=WF[:nodelabelcolor])
    Makie.text!(WF,celllabels, position=cellpositions, textsize=WF[:textsize], color=WF[:celllabelcolor], align=(:center,:center))
    #plot edges (3D) /faces (2D) of the mesh
    return Makie.linesegments!(WF,lines,color=WF[:color], linewidth=WF[:strokewidth], visible=WF[:visible])
end

function Makie.plot!(WF::Wireframe{<:Tuple{<:Ferrite.AbstractGrid{dim}}}) where dim
    grid = WF[1][]
    coords = [Ferrite.getcoordinates(node)[i] for node in Ferrite.getnodes(grid), i in 1:dim] 
    dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
    for cell in Ferrite.getcells(grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [coords[e,:] for boundary in boundaryentities for e in boundary])
    end
    nodes = lift(x->x ? coords : zeros(Float32,0,3),WF[:plotnodes])
    Makie.scatter!(WF,nodes,markersize=WF[:markersize], color=WF[:color])
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size(coords,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? [dim < 3 ? Point2f0(row) : Point3f0(row) for row in eachrow(coords)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:Ferrite.getncells(grid)] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,coords) for cell in Ferrite.getcells(grid)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
    Makie.text!(WF,nodelabels, position=nodepositions, textsize=WF[:textsize], offset=WF[:offset],color=WF[:nodelabelcolor])
    Makie.text!(WF,celllabels, position=cellpositions, textsize=WF[:textsize], color=WF[:celllabelcolor], align=(:center,:center))
    return Makie.linesegments!(WF,lines,color=WF[:color], strokewidth=WF[:strokewidth])
end

@recipe(Surface) do scene
    Attributes(
    field = :default,
    process = postprocess,
    scale_plot = false,
    shading = false,
    colormap = :viridis,
    )
end
 
function Makie.plot!(SF::Surface{<:Tuple{<:MakiePlotter{2}}})
    plotter = SF[1][]
    field = lift(x->x===:default ? 1 : Ferrite.find_field(plotter.dh,x),SF[:field])
    solution = lift((x,y)->reshape(dof_to_node(plotter.dh, plotter.u; field=x, process=y), Ferrite.getnnodes(plotter.dh.grid)),field,SF[:process])
    points = lift(x->[Makie.Point3f0(coord[1], coord[2], x[idx]) for (idx, coord) in enumerate(eachrow(plotter.coords))],solution)
    return Makie.mesh!(SF,points, reshape_triangles(plotter), color=solution, scale_plot=SF[:scale_plot], shading=SF[:shading], colormap=SF[:colormap])
end

@recipe(Arrows) do scene
    Attributes(
    arrowsize = 0.08,
    normalize = true,
    field = :default,
    color = :default,
    colormap = :viridis,
    process=postprocess,
    lengthscale = 1f0,
    )
end

function Makie.plot!(AR::Arrows{<:Tuple{<:MakiePlotter{dim}}}) where dim
    plotter = AR[1][]
    field = lift(x->x===:default ? 1 : Ferrite.find_field(plotter.dh,x),AR[:field])
    @assert Ferrite.getfielddim(plotter.dh,field[]) > 1
    solution = lift(x->dof_to_node(plotter.dh, plotter.u; field=x, process=identity),field)
    if dim  == 2
        ps = [Point2f0(i) for i in eachrow(plotter.coords)]
        ns = lift(x->[Vec2f(i) for i in eachrow(x)],solution)
        lengths = lift((x,y,z)-> z===:default ? x.(y) : ones(length(y))*z, AR[:process], ns, AR[:color])
    elseif dim  == 3
        ps = [Point3f0(i) for i in eachrow(plotter.coords)]
        ns = lift(x->[Vec3f(i) for i in eachrow(x)],solution)
        lengths = lift((x,y,z)-> z===:default ? x.(y) : ones(length(y))*z, AR[:process], ns, AR[:color])
    else
        error("Arrows plots are only available in dim ≥ 2")
    end
    Makie.arrows!(AR, ps, ns, arrowsize=AR[:arrowsize], normalize=AR[:normalize], colormap=AR[:colormap], color=lengths, lengthscale=AR[:lengthscale])
end

function ferriteviewer(plotter::MakiePlotter{dim}) where dim
    #set up figure and axis, axis either LScene in 3D or Axis if dim < 2
    fig = Figure()
    dim > 2 ? (ax = LScene(fig[1,1])) : (ax = Axis(fig[1,1]))
    #toggles and their labels for switching plot types on/off
    toggles = [Toggle(fig, active=active) for active in [true,true,false]]
    labels = [Label(fig,label) for label in ["mesh", "deformation", "labels"]]
    #setup the deformation_field as Node (Observable)
    deformation_field = Node(Ferrite.getfieldnames(plotter.dh)[1])
    #solutionplot main plot of the viewer
    solutionp = solutionplot!(plotter,colormap=:cividis,deformation_field=@lift $(toggles[2].active) ? $(deformation_field) : :default)

    #setting up various sliders
    markerslider = Slider(fig, range = 0:1:100, startvalue=5)
    strokewidthslider = Slider(fig, range = 0:1:10, startvalue=1)
    markersize = lift(x->x,markerslider.value)
    strokewidth = lift(x->x,strokewidthslider.value)

    #plot the fe-mesh
    wireframep = wireframe!(plotter,markersize=markersize,strokewidth=strokewidth,deformation_field= @lift $(toggles[2].active) ? $(deformation_field) : :default)
    #connect fe-mesh plot to the toggle
    connect!(wireframep.visible,toggles[1].active)
    connect!(wireframep.nodelabels,toggles[3].active)
    connect!(wireframep.celllabels,toggles[3].active)

    #set up dropdown menus for colormap, field, deformation field and processing function
    menu_cm = Menu(fig, options=["cividis", "inferno", "thermal"],label="colormap", direction=:up)
    menu_field = Menu(fig, options=Ferrite.getfieldnames(plotter.dh))
    menu_deformation_field = Menu(fig, options=Ferrite.getfieldnames(plotter.dh))
    menu_process = Menu(fig, options=[x₁,x₂,x₃,("magnitude",l2),("manhatten norm",l1),identity])
    #align all menus as a vgrid under each other
    fig[1,3] = vgrid!(grid!(hcat(toggles,labels), tellheight=false),
                      Label(fig,"nodesize",width=nothing), markerslider,
                      Label(fig,"strokewidth",width=nothing), strokewidthslider,
                      Label(fig,"processing function",width=nothing), menu_process,
                      Label(fig,"field",width=nothing), menu_field,
                      Label(fig, "deformation field",width=nothing),menu_deformation_field,
                      Label(fig, "colormap",width=nothing),menu_cm)
    #add colorbar
    cb = Colorbar(fig[1,2], solutionp)

    #event handling for selecting stuff from the menus
    on(menu_cm.selection) do s
        cb.colormap = s
        solutionp.colormap = s
    end

    on(menu_field.selection) do field
        solutionp.field = field
    end
    
    on(menu_deformation_field.selection) do field
        solutionp.deformation_field = field
        wireframep.deformation_field = field
    end

    on(menu_process.selection) do process_function
        solutionp.process=process_function
    end

    return fig
end
