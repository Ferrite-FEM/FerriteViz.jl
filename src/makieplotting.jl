"""
    FerriteVis.update!(plotter::MakiePlotter, u::Vector)
Updates the Observable `plotter.u` and thereby, triggers the plot to update.
"""
function update!(plotter::MakiePlotter, u::Vector)
    @assert length(plotter.u[]) == length(u)
    plotter.u[] .= u
    Makie.notify(plotter.u)
end

function Makie.convert_arguments(P::Type{<:Makie.Mesh}, plotter::MakiePlotter)
    return Makie.convert_arguments(P,plotter.physical_coords,plotter.triangles)
end

"""
    solutionplot(plotter::MakiePlotter; kwargs...)
    solutionplot(dh::AbstractDofHandler, u::Vector; kwargs...)
    solutionplot!(plotter::MakiePlotter; kwargs...)
    solutionplot!(dh::AbstractDofHandler, u::Vector; kwargs...)
Solutionplot produces the classical contour plot onto the finite element mesh. Most important
keyword arguments are:

- `field::Symbol=:default` representing the field which gets plotted, defaults to the first field in the `dh`.
- `deformation_field::Symbol=:default` field that transforms the mesh by the given deformation, defaults to no deformation
- `process::Function=postprocess` function to construct nodal scalar values from a vector valued problem
- `colormap::Symbol=:cividis`
- `deformation_scale=1.0`
- `shading=false`
- `scale_plot=false`
- `transparent=false`
"""
@recipe(SolutionPlot) do scene
    Attributes(
    scale_plot=false,
    shading=false,
    field=:default,
    deformation_field=:default,
    process=postprocess,
    colormap=:cividis,
    colorrange=(0,1),
    transparent=false,
    deformation_scale = 1.0,
    )
end

function Makie.plot!(SP::SolutionPlot{<:Tuple{<:MakiePlotter}})
    plotter = SP[1][]
    solution = @lift($(SP[:field])===:default ? reshape(transfer_solution(plotter,$(plotter.u); field_idx=1, process=$(SP[:process])), num_vertices(plotter)) : reshape(transfer_solution(plotter,$(plotter.u); field_idx=Ferrite.find_field(plotter.dh,$(SP[:field])), process=$(SP[:process])), num_vertices(plotter)))
    u_matrix = @lift($(SP[:deformation_field])===:default ? zeros(0,3) : transfer_solution(plotter,$(plotter.u); field_idx=Ferrite.find_field(plotter.dh,$(SP[:deformation_field])), process=identity))
    coords = @lift($(SP[:deformation_field])===:default ? plotter.physical_coords : plotter.physical_coords .+ ($(SP[:deformation_scale]) .* $(u_matrix)))
    mins = @lift(minimum($solution))
    maxs = @lift(maximum($solution))
    SP[:colorrange] = @lift(isapprox($mins,$maxs) ? (0,1e-8) : ($mins,$maxs))
    return Makie.mesh!(SP, coords, plotter.triangles, color=solution, shading=SP[:shading], scale_plot=SP[:scale_plot], colormap=SP[:colormap], transparent=SP[:transparent])
end

"""
    cellplot(plotter::MakiePlotter,σ::Vector{T}; kwargs...) where T
    cellplot!(plotter::MakiePlotter,σ::Vector{T}; kwargs...) where T
`cellplot` plots constant scalar data on the cells of the finite element mesh. If `T` is not a number, the keyword argument `process`
can be passed in order to reduce the elements of `σ` to a scalar.

keyword arguments are:

- `deformation_field::Symbol=:default` field that transforms the mesh by the given deformation, defaults to no deformation
- `process::Function=identity` function to construct cell scalar values. Defaults to `identity`, i.e. scalar values.
- `colormap::Symbol=:cividis`
- `deformation_scale=1.0`
- `shading=false`
- `scale_plot=false`
- `transparent=false`
"""
@recipe(CellPlot) do scene
    Attributes(
    scale_plot=false,
    shading=false,
    deformation_field=:default,
    process=identity,
    colormap=:cividis,
    colorrange=(0,1),
    transparent=false,
    deformation_scale = 1.0,
    )
end

function Makie.plot!(CP::CellPlot{<:Tuple{<:MakiePlotter{dim},Vector}}) where dim
    plotter = CP[1][]
    qp_values = CP[2][]
    u_matrix = @lift($(CP[:deformation_field])===:default ? zeros(0,3) : transfer_solution(plotter,$(plotter.u); field_idx=Ferrite.find_field(plotter.dh,$(CP[:deformation_field])), process=identity))
    coords = @lift($(CP[:deformation_field])===:default ? plotter.physical_coords : plotter.physical_coords .+ ($(CP[:deformation_scale]) .* $(u_matrix)))
    mins = minimum(qp_values)
    maxs = maximum(qp_values)
    CP[:colorrange] = isapprox(mins,maxs) ? (0,1e-8) : (mins,maxs) 
    solution =  @lift(reshape(transfer_scalar_celldata(plotter, qp_values; process=$(CP[:process])), num_vertices(plotter)))
    return Makie.mesh!(CP, coords, plotter.triangles, color=solution, shading=CP[:shading], scale_plot=CP[:scale_plot], colormap=CP[:colormap], transparent=CP[:transparent])
end

"""
    wireframe(plotter::MakiePlotter; kwargs...)
    wireframe(dh::AbstractDofHandler, u::Vector; kwargs...)
    wireframe(grid::AbstractGrid; kwargs...)
    wireframe!(plotter::MakiePlotter; kwargs...)
    wireframe!(dh::AbstractDofHandler, u::Vector; kwargs...)
    wireframe!(grid::AbstractGrid; kwargs...)
Plots the finite element mesh, optionally labels it and transforms it if a suitable `deformation_field` is given.

- `plotnodes::Bool=true` plots the nodes as circles/spheres
- `strokewidth::Int=2` how thick faces/edges are drawn
- `color::Symbol=theme(scene,:linecolor)` color of the faces/edges and nodes
- `markersize::Int=30` size of the nodes
- `deformation_field::Symbol=:default` field that transforms the mesh by the given deformation, defaults to no deformation
- `deformation_scale::Number=1.0` scaling of the deformation
- `nodelables=false` global node id labels
- `nodelabelcolor=:darkblue`
- `celllabels=false` global cell id labels
- `celllabelcolor=:darkred`
- `textsize::Int=15` size of the label's text
- `visible=true`
"""
@recipe(Wireframe) do scene
    Attributes(
    plotnodes=true,
    strokewidth=2,
    color=theme(scene, :linecolor),
    markersize=30,
    deformation_field=:default,
    visible=true,
    deformation_scale=1,
    textsize=15,
    offset=(0.0,0.0),
    nodelabels=false,
    nodelabelcolor=:darkblue,
    celllabels=false,
    celllabelcolor=:darkred,
    )
end

function Makie.plot!(WF::Wireframe{<:Tuple{<:MakiePlotter{dim}}}) where dim
    plotter = WF[1][]
    #triangle representation
    # can't use triangle representation, since we don't know by this information which edges to draw
    # further it would be in the incompressible example 2600 nodes vs 15000 in triangle representation
    # u_matrix = @lift($(WF[:deformation_field])===:default ? zeros(0,3) : transfer_solution(plotter; field_idx=Ferrite.find_field(plotter.dh,$(WF[:deformation_field])), process=identity))
    # coords = @lift($(WF[:deformation_field])===:default ? plotter.physical_coords : plotter.physical_coords .+ ($(WF[:scale]) .* $(u_matrix)))
    #original representation
    nodal_u_matrix = @lift($(WF[:deformation_field])===:default ? zeros(0,3) : dof_to_node(plotter.dh, $(WF[1][].u); field=Ferrite.find_field(plotter.dh,$(WF[:deformation_field])), process=identity))
    gridnodes = @lift($(WF[:deformation_field])===:default ? plotter.gridnodes : plotter.gridnodes .+ ($(WF[:deformation_scale]) .* $(nodal_u_matrix))) 
    lines = @lift begin
        dim > 2 ? (lines = Makie.Point3f0[]) : (lines = Makie.Point2f0[])
        for cell in Ferrite.getcells(plotter.dh.grid)
            boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
            append!(lines, [$gridnodes[e,:] for boundary in boundaryentities for e in boundary])
        end
        lines
    end
    nodes = @lift($(WF[:plotnodes]) ? $(gridnodes) : zeros(Float32,0,3))
    #plot the nodes
    Makie.scatter!(WF,gridnodes,markersize=WF[:markersize], color=WF[:color], visible=WF[:visible])
    #set up nodelabels
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size($gridnodes,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? [dim < 3 ? Point2f0(row) : Point3f0(row) for row in eachrow($gridnodes)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
    #set up celllabels
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:Ferrite.getncells(plotter.dh.grid)] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,$gridnodes) for cell in Ferrite.getcells(plotter.dh.grid)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
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
    nodes = @lift($(WF[:plotnodes]) ? coords : zeros(Float32,0,3))
    Makie.scatter!(WF,nodes,markersize=WF[:markersize], color=WF[:color])
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size(coords,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? [dim < 3 ? Point2f0(row) : Point3f0(row) for row in eachrow(coords)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:Ferrite.getncells(grid)] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,coords) for cell in Ferrite.getcells(grid)] : (dim < 3 ? [Point2f0((0,0))] : [Point3f0((0,0,0))])
    Makie.text!(WF,nodelabels, position=nodepositions, textsize=WF[:textsize], offset=WF[:offset],color=WF[:nodelabelcolor])
    Makie.text!(WF,celllabels, position=cellpositions, textsize=WF[:textsize], color=WF[:celllabelcolor], align=(:center,:center))
    return Makie.linesegments!(WF,lines,color=WF[:color], strokewidth=WF[:strokewidth])
end

"""
    surface(plotter::MakiePlotter; kwargs...)
    surface(dh::AbstractDofHandler, u::Vector; kwargs...)
    surface!(plotter::MakiePlotter; kwargs...)
    surface!(dh::AbstractDofHandler, u::Vector; kwargs...)
Uses the given `field` and plots the scalar values as a surface. If it's a vector valued problem, the nodal vector
values are transformed to a scalar based on `process` which defaults to the magnitude. Only availble in `dim=2`.
 
- `field = :default`
- `process = postprocess`
- `scale_plot = false`
- `shading = false`
- `colormap = :cividis`
"""
@recipe(Surface) do scene
    Attributes(
    field = :default,
    process = postprocess,
    scale_plot = false,
    shading = false,
    colormap = :cividis,
    )
end
 
function Makie.plot!(SF::Surface{<:Tuple{<:MakiePlotter{2}}})
    plotter = SF[1][]
    field = @lift($(SF[:field])===:default ? 1 : Ferrite.find_field(plotter.dh,$(SF[:field])))
    solution = @lift(reshape(transfer_solution(plotter, $(plotter.u); field_idx=$(field), process=$(SF[:process])), num_vertices(plotter)))
    points = @lift([Makie.Point3f0(coord[1], coord[2], $(solution)[idx]) for (idx, coord) in enumerate(eachrow(plotter.physical_coords))])
    return Makie.mesh!(SF,points, plotter.triangles, color=solution, scale_plot=SF[:scale_plot], shading=SF[:shading], colormap=SF[:colormap])
end

"""
    arrows(plotter::MakiePlotter; kwargs...)
    arrows(dh::AbstractDofHandler, u::Vector; kwargs...)
    arrows!(plotter::MakiePlotter; kwargs...)
    arrows!(dh::AbstractDofHandler, u::Vector; kwargs...)
At every node position a arrows is drawn, where the arrow tip ends at the node. Only works in `dim >=2`. If a `color` is specified
the arrows are unicolored. Otherwise the color corresponds to the magnitude, or any other scalar value based on the `process` function.

- `arrowsize = 0.08`
- `normalize = true`
- `field = :default`
- `color = :default`
- `colormap = :cividis`
- `process=postprocess`
- `lengthscale = 1f0`
"""
@recipe(Arrows) do scene
    Attributes(
    arrowsize = 0.08,
    normalize = true,
    field = :default,
    color = :default,
    colormap = :cividis,
    process=postprocess,
    lengthscale = 1f0,
    )
end

function Makie.plot!(AR::Arrows{<:Tuple{<:MakiePlotter{dim}}}) where dim
    plotter = AR[1][]
    field = @lift($(AR[:field])===:default ? 1 : Ferrite.find_field(plotter.dh,$(AR[:field])))
    @assert Ferrite.getfielddim(plotter.dh,field[]) > 1
    solution = @lift(transfer_solution(plotter, $(plotter.u); field_idx=$(field), process=identity))
    if dim  == 2
        ps = [Point2f0(i) for i in eachrow(plotter.physical_coords)]
        ns = @lift([Vec2f(i) for i in eachrow($(solution))])
        lengths = @lift($(AR[:color])===:default ? $(AR[:process]).($(ns)) : ones(length($(ns)))*$(AR[:color]))
    elseif dim  == 3
        ps = [Point3f0(i) for i in eachrow(plotter.physical_coords)]
        ns = @lift([Vec3f(i) for i in eachrow($(solution))])
        lengths = @lift($(AR[:color])===:default ? $(AR[:process]).($(ns)) : ones(length($(ns)))*$(AR[:color]))
    else
        error("Arrows plots are only available in dim ≥ 2")
    end
    Makie.arrows!(AR, ps, ns, arrowsize=AR[:arrowsize], normalize=AR[:normalize], colormap=AR[:colormap], color=lengths, lengthscale=AR[:lengthscale])
end

"""
    ferriteviewer(plotter::MakiePlotter)
    ferriteviewer(plotter::MakiePlotter, u_history::Vector{Vector{T}}})
Constructs a viewer with a `solutionplot`, `Colorbar` as well as sliders,toggles and menus to change the current view.
If the second dispatch is called a timeslider is added, in order to step through a set of solutions obtained from a simulation.
"""
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

function ferriteviewer(plotter::MakiePlotter, data::Vector{Vector{T}}) where T
    fig = ferriteviewer(plotter)
    timeslider = labelslider!(fig, "timestep n:", 1:length(data); format = x->"$x", sliderkw = Dict(:snap=>false))
    fig[2,1] = timeslider.layout
    @lift(FerriteVis.update!(plotter,data[$(timeslider.slider.value)]))
    display(fig)
end

####### One Shot Methods ####### 
const FerriteVisPlots = Union{Type{<:Wireframe},Type{<:SolutionPlot},Type{<:Arrows}}

function Makie.convert_arguments(P::FerriteVisPlots, dh::Ferrite.AbstractDofHandler, u::Vector)
    return (MakiePlotter(dh,u),)
end
