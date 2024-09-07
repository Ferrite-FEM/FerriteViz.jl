"""
    FerriteViz.update!(plotter::MakiePlotter, u::Vector)
Updates the Observable `plotter.u` and thereby, triggers the plot to update.
"""
function update!(plotter::MakiePlotter, u::Vector)
    @assert length(plotter.u[]) == length(u)
    plotter.u[] .= u
    Makie.notify(plotter.u)
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
- `colorrange::NTuple{2,<:Number}`: Specify (min, max) of the colorscale. If not given, min and max are calculated automatically from the data. 
- `deformation_scale=1.0`
- `shading=false`
- `scale_plot=false`
- `transparent=false`
- `nan_color::Union{Symbol, <:Colorant}=:red`
"""
@recipe(SolutionPlot) do scene
    Attributes(
    scale_plot=false,
    shading=false,
    field=:default,
    deformation_field=:default,
    process=postprocess,
    colormap=:cividis,
    colorrange=Makie.automatic,
    transparent=false,
    deformation_scale = 1.0,
    nan_color=:red,
    )
end

function Makie.plot!(SP::SolutionPlot{<:Tuple{<:MakiePlotter}})
    plotter = SP[1][]
    solution = @lift begin
        if $(SP[:field]) == :default
            field_name = Ferrite.getfieldnames(plotter.dh)[1]
            reshape(transfer_solution(plotter,$(plotter.u); field_name=field_name, process=$(SP[:process])), num_vertices(plotter))
        else
            reshape(transfer_solution(plotter,$(plotter.u); field_name=$(SP[:field]), process=$(SP[:process])), num_vertices(plotter))
        end
    end
    u_matrix = @lift begin
        if $(SP[:deformation_field])===:default
             Ferrite.getdim(plotter.dh.grid) > 2 ? Point3f[Point3f(0,0,0)] : Point2f[Point2f(0,0)]
        else
            #TODO remove convert
            convert(Vector{Point{Ferrite.getdim(plotter.dh.grid),Float32}},Makie.to_vertices(transfer_solution(plotter,$(plotter.u); field_name=$(SP[:deformation_field]), process=identity)))
        end
    end
    @lift begin
        if $(SP[:deformation_field])===:default
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords
        else
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords .+ ($(SP[:deformation_scale]) .* $(u_matrix))
        end
    end
    return Makie.mesh!(SP, plotter.mesh, color=solution, shading=SP[:shading], scale_plot=SP[:scale_plot], colormap=SP[:colormap],colorrange=SP[:colorrange] , transparent=SP[:transparent], nan_color=SP[:nan_color])
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
- `colorrange::NTuple{2,<:Number}`: Specify (min, max) of the colorscale. If not given, min and max are calculated automatically from the data. 
- `deformation_scale=1.0`
- `shading=false`
- `scale_plot=false`
- `transparent=false`
- `nan_color::Union{Symbol, <:Colorant}=:red`
"""
@recipe(CellPlot) do scene
    Attributes(
    scale_plot=false,
    shading=false,
    deformation_field=:default,
    process=identity,
    colormap=:cividis,
    colorrange=Makie.automatic,
    transparent=false,
    deformation_scale = 1.0,
    nan_color=:red,
    )
end

function Makie.plot!(CP::CellPlot{<:Tuple{<:MakiePlotter{dim},Vector}}) where dim
    plotter = CP[1][]
    qp_values = CP[2][]
    u_matrix = @lift begin
        if $(CP[:deformation_field])===:default
            Point3f[Point3f(0,0,0)]
        else
            convert(Vector{Point{Ferrite.getdim(plotter.dh.grid),Float32}},Makie.to_vertices(transfer_solution(plotter,$(plotter.u); field_name=$(CP[:deformation_field]), process=identity)))
        end
    end
    coords = @lift begin
        if $(CP[:deformation_field])===:default
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords
        else
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords .+ ($(CP[:deformation_scale]) .* $(u_matrix))
        end
    end
    solution =  @lift(reshape(transfer_scalar_celldata(plotter, qp_values; process=$(CP[:process])), num_vertices(plotter)))
    return Makie.mesh!(CP, plotter.mesh, color=solution, shading=CP[:shading], scale_plot=CP[:scale_plot], colormap=CP[:colormap], transparent=CP[:transparent], colorrange=CP[:colorrange], nan_color=CP[:nan_color])
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
- `cellsets=false` Color cells based on their cellset association. If no cellset is found for a cell, the cell is marked blue.
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
    color=theme(scene, :linecolor),
    strokewidth=theme(scene, :linewidth),
    markersize=theme(scene, :markersize),
    deformation_field=:default,
    visible=true,
    deformation_scale=1,
    textsize=15,
    offset=(0.0,0.0),
    nodelabels=false,
    nodelabelcolor=:darkblue,
    celllabels=false,
    celllabelcolor=:darkred,
    cellsets=false,
    depth_shift=-0.0001f0
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
    pointtype = dim > 2 ? Point3f : Point2f
    nodal_u_matrix = @lift begin 
        if $(WF[:deformation_field])===:default
            pointtype[zero(pointtype)]
        else
            convert(Vector{Point{Ferrite.getdim(plotter.dh.grid),Float32}},Makie.to_vertices(dof_to_node(plotter.dh, $(WF[1][].u); field_name=$(WF[:deformation_field]))))
        end
    end
    gridnodes = @lift begin
        if $(WF[:deformation_field])===:default
            plotter.gridnodes
        else
            plotter.gridnodes .+ ($(WF[:deformation_scale]) .* $(nodal_u_matrix))
        end
    end
    lines = @lift begin
        dim > 2 ? (lines = Point3f[]) : (lines = Point2f[])
        for cell in Ferrite.getcells(plotter.dh.grid)
            boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
            append!(lines, [$gridnodes[e] for boundary in boundaryentities for e in boundary])
        end
        lines
    end
    nodes = @lift($(WF[:plotnodes]) ? $(gridnodes) : pointtype[zero(pointtype)])
    #plot cellsets
    cellsets = plotter.dh.grid.cellsets
    cellset_to_value = Dict{String,Int}()
    for (cellsetidx,(cellsetname,cellset)) in enumerate(cellsets)
        cellset_to_value[cellsetname] = cellsetidx
    end
    cellset_u = zeros(Ferrite.getncells(plotter.dh.grid))
    for (cellidx,cell) in enumerate(Ferrite.getcells(plotter.dh.grid))
        for (cellsetname,cellsetvalue) in cellset_to_value
            if cellidx in cellsets[cellsetname]
                cellset_u[cellidx] = cellsetvalue
            end
        end
    end
    u_matrix = @lift begin
        if $(WF[:deformation_field])===:default
            pointtype[zero(pointtype)]
        else
            Makie.to_ndim.(pointtype, Makie.to_vertices(transfer_solution(plotter,$(plotter.u); field_name=$(WF[:deformation_field]), process=identity)), 0f0)
        end
    end
    coords = @lift begin
        if $(WF[:deformation_field])===:default
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords
        else
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords .+ ($(WF[:deformation_scale]) .* $(u_matrix))
        end
    end
    colorrange = isempty(cellset_to_value) ? (0,1) : (0,maximum(values(cellset_to_value)))
    cellset_u =  reshape(transfer_scalar_celldata(plotter, cellset_u; process=identity), num_vertices(plotter))
    Makie.mesh!(WF, plotter.mesh, color=cellset_u, shading=false, scale_plot=false, colormap=:darktest, visible=WF[:cellsets])
    #plot the nodes
    shouldplot = @lift ($(WF[:visible]) && $(WF[:plotnodes]))
    Makie.scatter!(WF,gridnodes,markersize=WF[:markersize], color=WF[:color], visible=shouldplot)
    #set up nodelabels
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size($gridnodes,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? $gridnodes : (dim < 3 ? Point2f[Point2f((0,0))] : Point3f[Point3f((0,0,0))])
    #set up celllabels
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:Ferrite.getncells(plotter.dh.grid)] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,$gridnodes) for cell in Ferrite.getcells(plotter.dh.grid)] : (dim < 3 ? [Point2f((0,0))] : [Point3f((0,0,0))])
    Makie.text!(WF,nodepositions, text=nodelabels, textsize=WF[:textsize], offset=WF[:offset],color=WF[:nodelabelcolor])
    Makie.text!(WF,celllabels, position=cellpositions, textsize=WF[:textsize], color=WF[:celllabelcolor], align=(:center,:center))
    #plot edges (3D) /faces (2D) of the mesh
    Makie.linesegments!(WF,lines,color=WF[:color], linewidth=WF[:strokewidth], visible=WF[:visible], depth_shift=WF[:depth_shift])
end


function Makie.plot!(WF::Wireframe{<:Tuple{<:Ferrite.AbstractGrid{dim}}}) where dim
    grid = WF[1][]
    coords = [Ferrite.getcoordinates(node)[i] for node in Ferrite.getnodes(grid), i in 1:dim] 
    coords = Makie.to_vertices(coords)
    dim > 2 ? (lines = Point3f[]) : (lines = Point2f[])
    for cell in Ferrite.getcells(grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        append!(lines, [coords[e] for boundary in boundaryentities for e in boundary])
    end
    nodes = @lift($(WF[:plotnodes]) ? coords : Point3f[Point3f(0,0,0)])
    shouldplot = @lift ($(WF[:visible]) && $(WF[:plotnodes]))
    Makie.scatter!(WF,nodes,markersize=WF[:markersize], color=WF[:color], visible=shouldplot)
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size(coords,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? coords : (dim < 3 ? Point2f[Point2f((0,0))] : Point3f[Point3f((0,0,0))])
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:Ferrite.getncells(grid)] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,coords) for cell in Ferrite.getcells(grid)] : (dim < 3 ? [Point2f((0,0))] : [Point3f((0,0,0))])
    #cellsetsplot
    if isconcretetype(grid.cells)
        dh = Ferrite.DofHandler(grid)
    else
        dh = Ferrite.MixedDofHandler(grid)
    end
    cellsets = grid.cellsets
    cellset_to_value = Dict{String,Int}()
    for (cellsetidx,(cellsetname,cellset)) in enumerate(cellsets)
        cellset_to_value[cellsetname] = cellsetidx
    end
    cellset_u = zeros(Ferrite.getncells(grid))
    for (cellidx,cell) in enumerate(Ferrite.getcells(grid))
        for (cellsetname,cellsetvalue) in cellset_to_value
            if cellidx in cellsets[cellsetname]
                cellset_u[cellidx] = cellsetvalue
            end
        end
    end
    plotter = MakiePlotter(dh,cellset_u)
    cellset_u =  reshape(transfer_scalar_celldata(plotter, cellset_u; process=identity), num_vertices(plotter))
    colorrange = isempty(cellset_to_value) ? (0,1) : (0,maximum(values(cellset_to_value)))
    Makie.mesh!(WF, plotter.mesh, color=cellset_u, shading=false, scale_plot=false, colormap=:darktest, visible=WF[:cellsets])
    Makie.text!(WF,nodelabels, position=nodepositions, textsize=WF[:textsize], offset=WF[:offset],color=WF[:nodelabelcolor])
    Makie.text!(WF,celllabels, position=cellpositions, textsize=WF[:textsize], color=WF[:celllabelcolor], align=(:center,:center))
    Makie.linesegments!(WF,lines,color=WF[:color], strokewidth=WF[:strokewidth], visible=WF[:visible])
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
- `colorrange=Makie.automatic`
- `nan_color::Union{Symbol, <:Colorant}=:red`
"""
@recipe(Surface) do scene
    Attributes(
    field = :default,
    process = postprocess,
    scale_plot = false,
    shading = false,
    colormap = :cividis,
    colorrange=Makie.automatic,
    nan_color=:red,
    )
end

function Makie.plot!(SF::Surface{<:Tuple{<:MakiePlotter{2}}})
    plotter = SF[1][]
    solution = @lift begin
        if $(SF[:field]) == :default
            field_name = Ferrite.getfieldnames(plotter.dh)[1]
            reshape(transfer_solution(plotter,$(plotter.u); field_name=field_name, process=$(SF[:process])), num_vertices(plotter))
        else
            reshape(transfer_solution(plotter,$(plotter.u); field_name=$(SF[:field]), process=$(SF[:process])), num_vertices(plotter))
        end
    end
    coords = @lift begin
        Point3f[Point3f(coord[1], coord[2], $(solution)[idx]) for (idx, coord) in enumerate(plotter.physical_coords)]
    end
    return Makie.mesh!(SF, coords, plotter.vis_triangles, color=solution, scale_plot=SF[:scale_plot], shading=SF[:shading], colormap=SF[:colormap], colorrange=SF[:colorrange], nan_color=SF[:nan_color])
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
    arrowsize = Makie.Automatic(),
    normalize = true, #TODO: broken
    field = :default,
    color = :default,
    colormap = :cividis,
    process=postprocess,
    lengthscale = 1f0,
    )
end

function Makie.plot!(AR::Arrows{<:Tuple{<:MakiePlotter{dim}}}) where dim
    plotter = AR[1][]
    solution = @lift begin
        if $(AR[:field]) === :default
            field_name = Ferrite.getfieldnames(plotter.dh)[1]
            @assert Ferrite.getfielddim(plotter.dh,field_name) > 1
            transfer_solution(plotter,$(plotter.u); field_name=field_name, process=identity)
        else
            @assert Ferrite.getfielddim(plotter.dh,$(AR[:field])) > 1
            transfer_solution(plotter,$(plotter.u); field_name=$(AR[:field]), process=identity)
        end
    end
    if dim  == 2
        ns = @lift([Vec2f(i) for i in eachrow($(solution))])
        lengths = @lift($(AR[:color])===:default ? $(AR[:process]).($(ns)) : ones(length($(ns)))*$(AR[:color]))
    elseif dim  == 3
        ns = @lift([Vec3f(i) for i in eachrow($(solution))])
        lengths = @lift($(AR[:color])===:default ? $(AR[:process]).($(ns)) : ones(length($(ns)))*$(AR[:color]))
    else
        error("Arrows plots are only available in dim ≥ 2")
    end
    Makie.arrows!(AR, plotter.physical_coords, ns, arrowsize=AR[:arrowsize], colormap=AR[:colormap], color=lengths, lengthscale=AR[:lengthscale])
end

"""
    elementinfo(ip::Interpolation; kwargs...)
    elementinfo(cell::AbstractCell; kwargs...)
    elementinfo(ip::Type{Interpolation}; kwargs...)
    elementinfo(cell::Type{AbstractCell}; kwargs...)

- `plotnodes=true` controls if nodes of element are plotted
- `strokewidth=2` strokwidth of faces/edges
- `color=theme(scene, :linecolor)`
- `markersize=30` size of the nodes
- `textsize=60` textsize of node-, edges- and facelabels
- `nodelabels=true` switch that controls plotting of nodelabels
- `nodelabelcolor=:darkred`
- `nodelabeloffset=(0.0,0.0)` offset of the nodelabel text relative to its associated node
- `facelabels=true` switch that controls plotting of facelabels
- `facelabelcolor=:darkgreen`
- `facelabeloffset=(-40,0)` offset of the facelabel text relative to its associated face middlepoint
- `edgelabels=true` switch that controls plotting of edgelabels
- `edgelabelcolor=:darkblue`
- `edgelabeloffset=(-40,-40)` offset of the edgelabel text relative to its associated edge middlepoint
- `font="Julia Mono"` font of the node-, edge-, and facelabels
"""
@recipe(Elementinfo) do scene
    Attributes(
    plotnodes=true,
    strokewidth=theme(scene, :linewidth),
    color=theme(scene, :linecolor),
    markersize=theme(scene, :markersize),
    textsize=60,
    nodelabels=true,
    nodelabelcolor=:darkred,
    nodelabeloffset=(0.0,0.0),
    facelabels=true,
    facelabelcolor=:darkgreen,
    facelabeloffset=(-40,0),
    edgelabels=true,
    edgelabelcolor=:darkblue,
    edgelabeloffset=(-40,-40),
    font=theme(scene,:font),
    )
end

function Makie.plot!(Ele::Elementinfo{<:Tuple{<:Ferrite.Interpolation{dim,refshape}}}) where {dim,refshape}
    ip = Ele[1][]
    elenodes = Ferrite.reference_coordinates(ip) |> x->reshape(reinterpret(Float64,x),(dim,length(x)))'
    dim > 2 ? (lines = Point3f[]) : (lines = Point2f[])
    facenodes = Ferrite.faces(ip)
    if dim == 2
        append!(lines, [elenodes[e,:] for boundary in facenodes for e in boundary[1:2]]) # 1:2 because higher order node in the middle
    else
        edgenodes = Ferrite.edges(ip)
        order = Ferrite.getorder(ip)
        #TODO remove the index monstrosity below after edges are defined consistently see https://github.com/Ferrite-FEM/Ferrite.jl/issues/520
        append!(lines, [elenodes[e,:] for boundary in edgenodes for e in boundary[1:((refshape == Ferrite.RefCube) ? 1 : (order > 1 ? 2 : 1)):((refshape == Ferrite.RefCube) ? 2 : end)]]) # 1:2 because higher order node in the middle
    end
    boundaryentities = dim == 2 ? facenodes : edgenodes
    #plot element boundary
    Makie.linesegments!(Ele,lines,color=Ele[:color], linewidth=Ele[:strokewidth])
    for (id,face) in enumerate(facenodes)
        idx = 0
        if refshape == Ferrite.RefCube && dim == 3
            idx = 4
        elseif refshape == Ferrite.RefTetrahedron && dim == 3
            idx = 3
        else
            idx = 2
        end
        position = zeros(dim)
        for i in 1:idx
            position += elenodes[face[i],:]
        end
        position ./= idx
        position = dim == 2 ? Point2f(position) : Point3f(position)
        Makie.text!(Ele,"$id", position=position, textsize=Ele[:textsize], offset=Ele[:facelabeloffset],color=Ele[:facelabelcolor],visible=Ele[:facelabels],font=Ele[:font])
    end
    if dim == 3
        for (id,edge) in enumerate(edgenodes)
            position = Point3f((elenodes[edge[1],:] + elenodes[refshape==Ferrite.RefCube ? edge[2] : edge[end],:])*0.5)
            t = Makie.text!(Ele,"$id", position=position, textsize=Ele[:textsize], offset=Ele[:edgelabeloffset],color=Ele[:edgelabelcolor],visible=Ele[:edgelabels],align=(:center,:center),font=Ele[:font])
            # Boundingbox can't switch currently from pixelspace to "coordinate" space in recipes
            #bb = Makie.boundingbox(t)
            #Makie.wireframe!(Ele,bb,space=:pixel)
        end
    end
    #plot the nodes
    Makie.scatter!(Ele,elenodes,markersize=Ele[:markersize], color=Ele[:color], visible=Ele[:plotnodes])
    #set up nodelabels
    nodelabels = @lift $(Ele[:nodelabels]) ? ["$i" for i in 1:size(elenodes,1)] : [""]
    nodepositions = @lift $(Ele[:nodelabels]) ? [dim < 3 ? Point2f(row) : Point3f(row) for row in eachrow(elenodes)] : (dim < 3 ? [Point2f((0,0))] : [Point3f((0,0,0))])
    #set up celllabels
    Makie.text!(Ele,nodelabels, position=nodepositions, textsize=Ele[:textsize], offset=Ele[:nodelabeloffset],color=Ele[:nodelabelcolor],font=Ele[:font])
    #plot edges (3D) /faces (2D) of the mesh
    Makie.linesegments!(Ele,lines,color=Ele[:color], linewidth=Ele[:strokewidth])
end

Makie.convert_arguments(P::Type{<:Elementinfo}, cell::C) where C<:Ferrite.AbstractCell = (Ferrite.default_interpolation(typeof(cell)),)
Makie.convert_arguments(P::Type{<:Elementinfo}, celltype::Type{C}) where C<:Ferrite.AbstractCell = (Ferrite.default_interpolation(celltype),)
Makie.convert_arguments(P::Type{<:Elementinfo}, iptype::Type{IP}) where IP<:Ferrite.Interpolation = (iptype(),)

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
    toggles = [Toggle(fig, active=active) for active in [true,false,false]]
    labels = [Label(fig,label) for label in ["mesh", "deformation", "labels"]]
    #setup the deformation_field as Observable
    deformation_field = Makie.Observable(Ferrite.getfieldnames(plotter.dh)[1])
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
    menu_cm = Menu(fig, options=["cividis", "inferno", "thermal"], direction=:up)
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
    sg = SliderGrid(fig[2,1], (label="timestep n:", range=1:length(data), format = x->"$x"))
    @lift(FerriteViz.update!(plotter,data[$(sg.sliders[1].value)]))
    display(fig)
end

#function initialize_figure(ip::Ferrite.Interpolation{ref_shape,N}) where {ref_shape<:Union{RefQuadrilateral,RefTriangle},N}
function initialize_figure(ip::Ferrite.Interpolation{2,ref_shape,N}) where {ref_shape<:Union{Ferrite.RefCube,Ferrite.RefTetrahedron},N}
    #(ip==CrouzeixRaviart{RefTriangle, 1}()) && error("To-Do: implement")
    rcs = Ferrite.reference_coordinates(ip)
    min_val = Int(min(vcat(rcs...)...))
    rcs = [[c[1]-min_val,c[2]-min_val] for c in rcs]
    scale = round(inv(min(filter(!iszero,vcat(rcs...))...));digits=15)
    max_id = Int(round(max(vcat(rcs...)...)*scale))
    refpos = [Int.([max_id-coord[2]+1, coord[1]+1]) for coord in rcs.*scale]

    fig = Figure()
    ax = [Axis3(fig[pos...];xlabel=L"ξ_x", ylabel=L"ξ_y", zlabel=L"N(ξ)" ,title="id: $(i), node = "*string(round.(rcs[i];digits=2))) for (i,pos) in enumerate(refpos)]
    return fig, ax
end

"""
    show_basis_function(ip::Interpolation{ref_shape,N})

Plots all basis functions of `ip`. Only implemented for 1D and 2D instances of `ip`.
"""
#function show_basis_function(ip::Ferrite.Interpolation{ref_shape,N}) where {ref_shape<:RefTriangle,N}
function show_basis_function(ip::Ferrite.Interpolation{2,ref_shape,N}; meshsize=20, plotnodes=true, strokewidth=1, markersize=10, fontsize=60, nodelabels=true, nodelabelcolor=:darkred, nodelabeloffset=(2.0,2.0), facelabels=false, facelabelcolor=:darkgreen, facelabeloffset=(-40,0), edgelabels=true, edgelabelcolor=:darkblue, edgelabeloffset=(-40,-40), font="Julia Mono") where {ref_shape<:Ferrite.RefTetrahedron,N}
    rcs = Ferrite.reference_coordinates(ip)
    mn = isa(ip,Ferrite.CrouzeixRaviart) ? 0. : min(vcat(collect.(rcs)...)...)
    mx = isa(ip,Ferrite.CrouzeixRaviart) ? 1. : mx = max(vcat(collect.(rcs)...)...)
    y = range(mx,mn,length=meshsize)
    x = range(mn,mx,length=meshsize)

    # initialize axes
    fig,ax = initialize_figure(ip)

    for i in 1:Ferrite.getnbasefunctions(ip)
        local z = [(ix>iy) ? NaN : Ferrite.value(ip,i,Ferrite.Vec{2,Float64}((_x,_y))) for (iy,_y) in enumerate(y), (ix,_x) in enumerate(x)]
        vertices, clist = get_triangulation(ip,x,y,z)
        !isa(ip,Ferrite.CrouzeixRaviart) && elementinfo!(ax[i],ip;plotnodes, strokewidth, markersize, fontsize, nodelabels, nodelabelcolor, nodelabeloffset, facelabels, facelabelcolor, facelabeloffset, edgelabels, edgelabelcolor, edgelabeloffset, font)
        mesh!(ax[i],vertices,clist; shading=false, fxaa=true, transparency=false, color=[vertices[i,3] for i in 1:size(vertices)[1]], colormap=:viridis)
    end
    current_figure()
end

#function show_basis_function(ip::Ferrite.Interpolation{ref_shape,N}) where {ref_shape<:RefQuadrilateral,N}
function show_basis_function(ip::Ferrite.Interpolation{2,ref_shape,N}; meshsize=20, plotnodes=true, strokewidth=1, markersize=10, fontsize=60, nodelabels=true, nodelabelcolor=:darkred, nodelabeloffset=(2.0,2.0), facelabels=false, facelabelcolor=:darkgreen, facelabeloffset=(-40,0), edgelabels=true, edgelabelcolor=:darkblue, edgelabeloffset=(-40,-40), font="Julia Mono") where {ref_shape<:Ferrite.RefCube,N}
    rcs = Ferrite.reference_coordinates(ip)
    xy = range(min(vcat(rcs...)...),max(vcat(rcs...)...),length=meshsize)

    # initialize axes
    fig,ax = initialize_figure(ip)

    for i in 1:Ferrite.getnbasefunctions(ip)
        z = [Ferrite.value(ip,i,Ferrite.Vec{2,Float64}((x,y))) for x in xy, y in xy]
        elementinfo!(ax[i],ip;plotnodes, strokewidth, markersize, fontsize, nodelabels, nodelabelcolor, nodelabeloffset, facelabels, facelabelcolor, facelabeloffset, edgelabels, edgelabelcolor, edgelabeloffset, font)
        Makie.surface!(ax[i],xy,xy,z)
        #mesh!(ax[i],vertices,clist; shading=false, fxaa=true, transparency=false, color=[vertices[i,3] for i in 1:size(vertices)[1]], colormap=:viridis)
    end
    current_figure()
end

#function show_basis_function(ip::Ferrite.Lagrange{RefLine,N}) where {N}
function show_basis_function(ip::Ferrite.Lagrange{1,refshape,N}) where {refshape,N}
    rcs = Ferrite.reference_coordinates(ip)
    x = range(-1,1,length=100)

    # initialize axes
    min_val = Int(min(vcat(rcs...)...))
    rcs = [[c[1]-min_val] for c in rcs]
    scale = round(inv(min(filter(!iszero,vcat(rcs...))...));digits=15)
    max_id = Int(round(max(vcat(rcs...)...)*scale))
    refpos = [Int.([coord[1]+1]) for coord in rcs.*scale]

    fig = Figure()
    ax = [Axis(fig[1,i];xlabel=L"ξ", ylabel=L"N(ξ)" ,title="id: $(i), node = "*string(round.(rcs[i];digits=2))) for i in 1:Ferrite.getnbasefunctions(ip)]
    for i in 1:Ferrite.getnbasefunctions(ip)
        local z = [Ferrite.value(ip,i,Ferrite.Vec{1,Float64}((_x,))) for _x in x]
        lines!(ax[refpos[i]...],x,z)
    end
    current_figure()
end

#function get_triangulation(ip::Ferrite.Interpolation{ref_shape,N},x,_y,_z::Matrix) where {ref_shape,N}
function get_triangulation(ip::Ferrite.Interpolation{2,ref_shape,N},x,y,z::Matrix) where {ref_shape,N}
# ================================
# ===== compute vertice list =====
# ================================
    sz = size(z)
    nvertices_percol = (.!isnan.(z)')*ones(Int,sz[2])
    nvertices = sum(nvertices_percol)
    vertices = zeros(nvertices,3)                           #   vertice numbering according to numbering:
    cnt = 1                                                 #   ┌                   ┐      ┌                   ┐                ┌                   ┐
    for col in 1:sz[2], row in 1:sz[1]                      #   │ 1 NaN NaN NaN NaN │      │ 1 NaN NaN NaN NaN │  not allowed:  │ 1 NaN NaN NaN NaN │
        if !isnan(z[row,col])                               #   │ 2  6  NaN NaN NaN │      │ 2  6  NaN 13  17  │ -connectivity  │ 2  6  NaN 13  16  │
            vertices[cnt,:] = [x[col], y[row], z[row,col]]  #   │ 3  7  10  NaN NaN │  or  │ 3  7  10  14  18  │  list is not   │ 3  7  10  14  17  │
            cnt += 1                                        #   │ 4  8  11  13  NaN │      │ 4  8  11  15  19  │  implemented   │ 4  8  11  15  18  │
        end                                                 #   │ 5  9  12  14  15  │      │ 5  9  12  16  20  │  for this      │ 5  9  12  NaN NaN │
    end                                                     #   └                          └                                    └                    
# =====================================
# ===== compute connectivity list =====
# =====================================
    ntriangles = 0 # number of triangles necessary
    for (i,vertincol) in enumerate(nvertices_percol)
        if i!=length(nvertices_percol)
            ntriangles+=(vertincol-1)*2 # if current collum has same length as next one add 2 triangles per "gap" between two nodes
            if vertincol==nvertices_percol[i+1]-1 # if next col is larger add 1
                ntriangles+=1
            elseif vertincol==nvertices_percol[i+1]+1 # if next col is smaller subtract 1
                ntriangles+=-1
            else
                !(vertincol==nvertices_percol[i+1]) ? error("triangulation not possible") : nothing
            end
        end
    end
    clist = zeros(Int,ntriangles,3)

    cnt = 1
    lastvertidpercol = [sum(nvertices_percol[1:c]) for c in 1:sz[2]]
    for vert in 1:lastvertidpercol[end-1]
        newcol = vert in vcat(1,lastvertidpercol[1:end-1].+1) # vertex first id in a column??
        col = findfirst(i->vert<=i,lastvertidpercol)
        if vert==sum(nvertices_percol[1:col]) # if vert is last vert in collumn...
            if nvertices_percol[col] == 1 # ... check if its the only vertex of the collumn...
                clist[cnt,:] = [vert, vert+1, vert+2] # ...add triangle if necessary...
                cnt+=1                
            end
            continue # ... skip last vertex per column
        end
        # fill clist
        offset = nvertices_percol[col]==nvertices_percol[col+1]-1 ? nvertices_percol[col]+1 : nvertices_percol[col] #set ofset according to next collumn size
        clist[cnt,:] = [vert, vert+1, vert+offset]
        cnt+=1
        if ((nvertices_percol[col]==nvertices_percol[col+1]+1)&&(!newcol)) || ((nvertices_percol[col]==nvertices_percol[col+1]-1)&&(newcol)) # if next column is shorter && not first vertex of collumn or if next column is longer and first vertex in collumn
            clist[cnt,:] = [vert, vert+offset-1, vert+offset]
            cnt+=1
        end
        if (nvertices_percol[col]==nvertices_percol[col+1]) || (nvertices_percol[col]==nvertices_percol[col+1]-1) # if next row is longer or shorter
            clist[cnt,:] = [vert+1, vert+offset+1, vert+offset]
            cnt+=1
        end
    end
(cnt-1)!=ntriangles && error("triangulations failed. too few triangles added to connectivity list")

    return vertices, clist
end

####### One Shot Methods #######
const FerriteVizPlots = Union{Type{<:Wireframe},Type{<:SolutionPlot},Type{<:Arrows},Type{<:Surface}}

function Makie.convert_arguments(P::FerriteVizPlots, dh::Ferrite.AbstractDofHandler, u::Vector)
    return (MakiePlotter(dh,u),)
end
