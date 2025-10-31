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
- `shading=Makie.NoShading`
- `nan_color::Union{Symbol, <:Colorant}=:red`
"""
@recipe(SolutionPlot) do scene
    Attributes(
    shading=Makie.NoShading,
    field=:default,
    deformation_field=:default,
    process=postprocess,
    colormap=:cividis,
    colorrange=Makie.automatic,
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
        grid = Ferrite.get_grid(plotter.dh)
        sdim = Ferrite.getspatialdim(grid)
        if $(SP[:deformation_field])===:default
            sdim > 2 ? Point3f[Point3f(0,0,0)] : Point2f[Point2f(0,0)]
        else
            #TODO remove convert
            convert(Vector{Point{sdim,Float32}},Makie.to_vertices(transfer_solution(plotter,$(plotter.u); field_name=$(SP[:deformation_field]), process=identity)))
        end
    end
    @lift begin
        if $(SP[:deformation_field])===:default
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords
        else
            plotter.physical_coords_mesh[1:end] = plotter.physical_coords .+ ($(SP[:deformation_scale]) .* $(u_matrix))
        end
    end
    return Makie.mesh!(SP, plotter.mesh, color=solution, shading=SP[:shading], colormap=SP[:colormap],colorrange=SP[:colorrange], nan_color=SP[:nan_color])
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
- `shading=Makie.NoShading`
- `nan_color::Union{Symbol, <:Colorant}=:red`
"""
@recipe(CellPlot) do scene
    Attributes(
    shading=Makie.NoShading,
    deformation_field=:default,
    process=identity,
    colormap=:cividis,
    colorrange=Makie.automatic,
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
            grid = Ferrite.get_grid(plotter.dh)
            sdim = Ferrite.getspatialdim(grid)
            convert(Vector{Point{sdim,Float32}},Makie.to_vertices(transfer_solution(plotter,$(plotter.u); field_name=$(CP[:deformation_field]), process=identity)))
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
    return Makie.mesh!(CP, plotter.mesh, color=solution, shading=CP[:shading], colormap=CP[:colormap], colorrange=CP[:colorrange], nan_color=CP[:nan_color])
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
- `linewidth::Int=2` how thick faces/edges are drawn
- `color::Symbol=theme(scene,:linecolor)` color of the faces/edges and nodes
- `markersize::Int=30` size of the nodes
- `deformation_field::Symbol=:default` field that transforms the mesh by the given deformation, defaults to no deformation
- `deformation_scale::Number=1.0` scaling of the deformation
- `cellsets=false` Color cells based on their cellset association. If no cellset is found for a cell, the cell is marked blue.
- `nodelables=false` global node id labels
- `nodelabelcolor=:darkblue`
- `celllabels=false` global cell id labels
- `celllabelcolor=:darkred`
- `fontsize::Int=15` size of the label's text
- `visible=true`
"""
@recipe(Wireframe) do scene
    Attributes(
    plotnodes=true,
    color=theme(scene, :linecolor),
    linewidth=theme(scene, :linewidth),
    markersize=theme(scene, :markersize),
    deformation_field=:default,
    visible=true,
    deformation_scale=1,
    fontsize=15,
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
            grid = Ferrite.get_grid(plotter.dh)
            sdim = Ferrite.getspatialdim(grid)
            # Float32 to be performant on a larger range of GPUs
            convert(Vector{Point{sdim,Float32}},Makie.to_vertices(dof_to_node(plotter.dh, $(WF[1][].u); field_name=$(WF[:deformation_field]))))
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
        grid = Ferrite.get_grid(plotter.dh)
        for cell in Ferrite.getcells(grid)
            boundaryentities = Ferrite.edges(cell)
            append!(lines, [$gridnodes[e] for boundary in boundaryentities for e in boundary])
        end
        lines
    end
    nodes = @lift($(WF[:plotnodes]) ? $(gridnodes) : pointtype[zero(pointtype)])
    #plot cellsets
    grid = Ferrite.get_grid(plotter.dh)
    cellsets = grid.cellsets
    cellset_to_value = Dict{String,Int}()
    for (cellsetidx,(cellsetname,cellset)) in enumerate(cellsets)
        cellset_to_value[cellsetname] = cellsetidx
    end
    cellset_u = zeros(Ferrite.getncells(grid))
    allcells = Ferrite.getcells(grid)
    ncells = Ferrite.getncells(grid)
    for (cellidx,cell) in enumerate(allcells)
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
    Makie.mesh!(WF, plotter.mesh, color=cellset_u, shading=Makie.NoShading, colormap=:darktest, visible=WF[:cellsets])
    #plot the nodes
    shouldplot = @lift ($(WF[:visible]) && $(WF[:plotnodes]))
    Makie.scatter!(WF,gridnodes,markersize=WF[:markersize], color=WF[:color], visible=shouldplot)
    #set up nodelabels
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size($gridnodes,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? $gridnodes : (dim < 3 ? Point2f[Point2f((0,0))] : Point3f[Point3f((0,0,0))])
    #set up celllabels
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:ncells] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,$gridnodes) for cell in allcells] : (dim < 3 ? [Point2f((0,0))] : [Point3f((0,0,0))])
    Makie.text!(WF,nodepositions, text=nodelabels, fontsize=WF[:fontsize], offset=WF[:offset],color=WF[:nodelabelcolor])
    Makie.text!(WF,celllabels, position=cellpositions, fontsize=WF[:fontsize], color=WF[:celllabelcolor], align=(:center,:center))
    #plot edges (3D) /faces (2D) of the mesh
    @show lines, "-----------------------------"
    Makie.linesegments!(WF,lines,color=WF[:color], linewidth=WF[:linewidth], visible=WF[:visible], depth_shift=WF[:depth_shift])
end


function Makie.plot!(WF::Wireframe{<:Tuple{<:Ferrite.AbstractGrid{dim}}}) where dim
    @info WF
    grid   = WF[1][]
    coords = [Ferrite.get_node_coordinate(node)[i] for node in Ferrite.getnodes(grid), i in 1:dim] 
    coords = Makie.to_vertices(coords)
    dim > 2 ? (lines = Point3f[]) : (lines = Point2f[])
    allcells = Ferrite.getcells(grid)
    ncells   = Ferrite.getncells(grid)
    for cell in allcells
        boundaryentities = Ferrite.edges(cell)
        append!(lines, [coords[e] for boundary in boundaryentities for e in boundary])
    end
    nodes = @lift($(WF[:plotnodes]) ? coords : Point3f[Point3f(0,0,0)])
    shouldplot = @lift ($(WF[:visible]) && $(WF[:plotnodes]))
    Makie.scatter!(WF,nodes,markersize=WF[:markersize], color=WF[:color], visible=shouldplot)
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in 1:size(coords,1)] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? coords : (dim < 3 ? Point2f[Point2f((0,0))] : Point3f[Point3f((0,0,0))])
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in 1:ncells] : [""]
    cellpositions = @lift $(WF[:celllabels]) ? [midpoint(cell,coords) for cell in allcells] : (dim < 3 ? [Point2f((0,0))] : [Point3f((0,0,0))])
    #cellsetsplot

    dh = Ferrite.DofHandler(grid)
    cellsets = grid.cellsets
    cellset_to_value = Dict{String,Int}()
    for (cellsetidx,(cellsetname,cellset)) in enumerate(cellsets)
        cellset_to_value[cellsetname] = cellsetidx
    end
    cellset_u = zeros(ncells)
    for (cellidx,cell) in enumerate(allcells)
        for (cellsetname,cellsetvalue) in cellset_to_value
            if cellidx in cellsets[cellsetname]
                cellset_u[cellidx] = cellsetvalue
            end
        end
    end
    plotter = MakiePlotter(dh,cellset_u)
    cellset_u =  reshape(transfer_scalar_celldata(plotter, cellset_u; process=identity), num_vertices(plotter))
    colorrange = isempty(cellset_to_value) ? (0,1) : (0,maximum(values(cellset_to_value)))
    Makie.mesh!(WF, plotter.mesh, color=cellset_u, shading=Makie.NoShading, colormap=:darktest, visible=WF[:cellsets])
    Makie.text!(WF,nodelabels, position=nodepositions, fontsize=WF[:fontsize], offset=WF[:offset],color=WF[:nodelabelcolor])
    Makie.text!(WF,celllabels, position=cellpositions, fontsize=WF[:fontsize], color=WF[:celllabelcolor], align=(:center,:center))
    Makie.linesegments!(WF,lines,color=WF[:color], linewidth=WF[:linewidth], visible=WF[:visible])
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
- `shading = Makie.NoShading`
- `colormap = :cividis`
- `colorrange=Makie.automatic`
- `nan_color::Union{Symbol, <:Colorant}=:red`
"""
@recipe(Surface) do scene
    Attributes(
    field = :default,
    process = postprocess,
    shading = Makie.NoShading,
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
    return Makie.mesh!(SF, coords, plotter.vis_triangles, color=solution, shading=SF[:shading], colormap=SF[:colormap], colorrange=SF[:colorrange], nan_color=SF[:nan_color])
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

function Makie.plot!(AR::Arrows{<:Tuple{<:MakiePlotter{sdim}}}) where sdim
    plotter = AR[1][]
    solution = @lift begin
        if $(AR[:field]) === :default
            field_name = Ferrite.getfieldnames(plotter.dh)[1]
            field_dim = Ferrite.n_components(plotter.dh, field_name)
            @assert field_dim == sdim "Dimension of field $field_name is $field_dim does not match spatial dimension $sdim"
            transfer_solution(plotter,$(plotter.u); field_name=field_name, process=identity)
        else
            field_name = $(AR[:field])
            field_dim = Ferrite.n_components(plotter.dh, field_name)
            @assert field_dim == sdim "Dimension of field $field_name is $field_dim does not match spatial dimension $sdim"
            transfer_solution(plotter,$(plotter.u); field_name=$(AR[:field]), process=identity)
        end
    end
    if sdim  == 2
        ns = @lift([Vec2f(i) for i in eachrow($(solution))])
        lengths = @lift($(AR[:color])===:default ? $(AR[:process]).($(ns)) : ones(length($(ns)))*$(AR[:color]))
    elseif sdim  == 3
        ns = @lift([Vec3f(i) for i in eachrow($(solution))])
        lengths = @lift($(AR[:color])===:default ? $(AR[:process]).($(ns)) : ones(length($(ns)))*$(AR[:color]))
    else
        error("Arrows plots are only available in dimension ≥ 2")
    end
    Makie.arrows!(AR, plotter.physical_coords, ns, arrowsize=AR[:arrowsize], colormap=AR[:colormap], color=lengths, lengthscale=AR[:lengthscale])
end

"""
    elementinfo(ip::Interpolation; kwargs...)
    elementinfo(cell::AbstractCell; kwargs...)
    elementinfo(ip::Type{Interpolation}; kwargs...)
    elementinfo(cell::Type{AbstractCell}; kwargs...)

- `plotnodes=true` controls if nodes of element are plotted
- `linewidth=2` strokwidth of faces/edges
- `color=theme(scene, :linecolor)`
- `markersize=30` size of the nodes
- `fontsize=60` fontsize of node-, edges- and facelabels
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
    linewidth=theme(scene, :linewidth),
    color=theme(scene, :linecolor),
    markersize=theme(scene, :markersize),
    fontsize=60,
    vertexlabels=true,
    vertexlabelcolor=:darkred,
    vertexlabeloffset=(0.0,0.0),
    nodelabels=true,
    nodelabelcolor=:darkred,
    nodelabeloffset=(0.0,20.0),
    facelabels=true,
    facelabelcolor=:darkgreen,
    facelabeloffset=(-40,0),
    edgelabels=true,
    edgelabelcolor=:darkblue,
    edgelabeloffset=(-40,-40),
    font=theme(scene,:font),
    )
end

function Makie.plot!(Ele::Elementinfo{<:Tuple{<:Ferrite.AbstractCell{refshape}}}) where {refshape}
    cell = Ele[1][]
    dim = Ferrite.getrefdim(cell)

    # Draw element outline
    gip = Ferrite.geometric_interpolation(cell)
    elenodes = Ferrite.reference_coordinates(gip) |> x->reshape(reinterpret(Float64,x),(dim,length(x)))'
    dim > 2 ? (lines = Point3f[]) : (lines = Point2f[])
    for edgenodes in Ferrite.edgedof_indices(gip)
        append!(lines, [elenodes[edgenodes[1],:], elenodes[edgenodes[2],:]]) # The convention in Ferrite is that the first two edge nodes are always associated to the vertices
    end
    Makie.linesegments!(Ele,lines,color=Ele[:color], linewidth=Ele[:linewidth])

    # Draw its nodes
    Makie.scatter!(Ele,elenodes,markersize=Ele[:markersize], color=Ele[:color], visible=Ele[:plotnodes])
    nodelabels = @lift $(Ele[:nodelabels]) ? ["N$i" for i in 1:size(elenodes,1)] : [""]
    nodepositions = @lift $(Ele[:nodelabels]) ? [dim < 3 ? Point2f(row) : Point3f(row) for row in eachrow(elenodes)] : (dim < 3 ? [Point2f((0,0))] : [Point3f((0,0,0))])
    Makie.text!(Ele,nodelabels, position=nodepositions, fontsize=Ele[:fontsize], offset=Ele[:nodelabeloffset],color=Ele[:nodelabelcolor],font=Ele[:font])

    # Annotate element vertices
    if dim ≥ 1 && Ele[:vertexlabels][]
        for (id,vertexnodes) in enumerate(Ferrite.vertexdof_indices(gip))
            position = if dim == 3
                Point3f(elenodes[vertexnodes[1],:])
            elseif dim == 2
                Point2f(elenodes[vertexnodes[1],:])
            end
            Makie.text!(Ele,"V$id", position=position, fontsize=Ele[:fontsize], offset=Ele[:edgelabeloffset],color=Ele[:edgelabelcolor],visible=Ele[:edgelabels],font=Ele[:font])
        end
    end

    # Annotate element edges
    if dim ≥ 2 && Ele[:edgelabels][]
        for (id,edgenodes) in enumerate(Ferrite.edgedof_indices(gip))
            position = if dim == 3
                Point3f((elenodes[edgenodes[1],:] + elenodes[edgenodes[2],:])*0.5)
            elseif dim == 2
                Point2f((elenodes[edgenodes[1],:] + elenodes[edgenodes[2],:])*0.5)
            end
            Makie.text!(Ele,"E$id", position=position, fontsize=Ele[:fontsize], offset=Ele[:edgelabeloffset],color=Ele[:edgelabelcolor],visible=Ele[:edgelabels],font=Ele[:font])
        end
    end

    # Annotate element faces
    if dim ≥ 3 && Ele[:facelabels][]
        for face_index ∈ 1:Ferrite.nfaces(cell)
            linear_face = linear_face_cell(cell, face_index)
            gip_face = Ferrite.geometric_interpolation(linear_face)
            vertexdof_list = Ferrite.vertexdof_indices(gip_face)
            refpositions   = Ferrite.reference_coordinates(gip_face)
            facenodes      = [Ferrite.facet_to_element_transformation(refposition, refshape, face_index) for refposition in refpositions] |> x->reshape(reinterpret(Float64,x),(dim,length(x)))'        
            position = zeros(dim)
            for i in 1:length(vertexdof_list)
                position += facenodes[vertexdof_list[i][1],:]
            end
            position ./= length(vertexdof_list)
            position = dim == 2 ? Point2f(position) : Point3f(position)
            Makie.text!(Ele,"F$face_index", position=position, fontsize=Ele[:fontsize], offset=Ele[:facelabeloffset],color=Ele[:facelabelcolor],visible=Ele[:facelabels],font=Ele[:font])
        end
    end
end

function Makie.plot!(Ele::Elementinfo{<:Tuple{<:Ferrite.Interpolation{refshape}}}) where {refshape}
    ip = Ele[1][]
    dim = Ferrite.getrefdim(ip)

    # Draw element boundary
    gip = Ferrite.default_geometric_interpolation(ip)
    geonodes = Ferrite.reference_coordinates(gip) |> x->reshape(reinterpret(Float64,x),(dim,length(x)))'
    dim > 2 ? (lines = Point3f[]) : (lines = Point2f[])
    for edgenodes in Ferrite.edgedof_indices(gip.ip)
        append!(lines, [geonodes[edgenodes[1],:], geonodes[edgenodes[2],:]]) 
    end
    Makie.linesegments!(Ele,lines,color=Ele[:color], linewidth=Ele[:linewidth])

    # Draw the interpolation nodes and its indices
    elenodes = Ferrite.reference_coordinates(ip) |> x->reshape(reinterpret(Float64,x),(dim,length(x)))'
    Makie.scatter!(Ele,elenodes,markersize=Ele[:markersize], color=Ele[:color], visible=Ele[:plotnodes])
    nodelabels = @lift $(Ele[:nodelabels]) ? ["D$i" for i in 1:size(elenodes,1)] : [""]
    nodepositions = @lift $(Ele[:nodelabels]) ? [dim < 3 ? Point2f(row) : Point3f(row) for row in eachrow(elenodes)] : (dim < 3 ? [Point2f((0,0))] : [Point3f((0,0,0))])
    Makie.text!(Ele,nodelabels, position=nodepositions, fontsize=Ele[:fontsize], offset=Ele[:nodelabeloffset],color=Ele[:nodelabelcolor],font=Ele[:font])

    # Annotate element vertices
    if dim ≥ 1 && Ele[:vertexlabels][]
        for (id,vertexnodes) in enumerate(Ferrite.vertexdof_indices(gip.ip))
            position = if dim == 3
                Point3f(elenodes[vertexnodes[1],:])
            elseif dim == 2
                Point2f(elenodes[vertexnodes[1],:])
            end
            Makie.text!(Ele,"V$id", position=position, fontsize=Ele[:fontsize], offset=Ele[:edgelabeloffset],color=Ele[:edgelabelcolor],visible=Ele[:edgelabels],font=Ele[:font])
        end
    end
    # Annotate element edges
    if dim ≥ 2 && Ele[:edgelabels][]
        for (id,edgenodes) in enumerate(Ferrite.edgedof_indices(gip.ip))
            position = if dim == 3
                Point3f((elenodes[edgenodes[1],:] + elenodes[edgenodes[2],:])*0.5)
            elseif dim == 2
                Point2f((elenodes[edgenodes[1],:] + elenodes[edgenodes[2],:])*0.5)
            end
            Makie.text!(Ele,"E$id", position=position, fontsize=Ele[:fontsize], offset=Ele[:edgelabeloffset],color=Ele[:edgelabelcolor],visible=Ele[:edgelabels],font=Ele[:font])
        end
    end

    # Annotate element faces
    if dim ≥ 3 && Ele[:facelabels][]
        for face_index ∈ 1:Ferrite.nfaces(ip)
            linear_face = linear_face_cell(refshape, face_index)
            gip_face = Lagrange{linear_face,1}()
            vertexdof_list = Ferrite.vertexdof_indices(gip_face)
            refpositions   = Ferrite.reference_coordinates(gip_face)
            facenodes      = [Ferrite.facet_to_element_transformation(refposition, refshape, face_index) for refposition in refpositions] |> x->reshape(reinterpret(Float64,x),(dim,length(x)))'        
            position = zeros(dim)
            for i in 1:length(vertexdof_list)
                position += facenodes[vertexdof_list[i][1],:]
            end
            position ./= length(vertexdof_list)
            position = dim == 2 ? Point2f(position) : Point3f(position)
            Makie.text!(Ele,"F$face_index", position=position, fontsize=Ele[:fontsize], offset=Ele[:facelabeloffset],color=Ele[:facelabelcolor],visible=Ele[:facelabels],font=Ele[:font])
        end
    end
end

function Makie.convert_arguments(P::Type{<:Elementinfo}, celltype::Type{C}) where C<:Ferrite.AbstractCell 
    gip = geometric_interpolation(C)
    nnodes = getnbasefunctions(gip)
    nodes  = ntuple(x->1,nnodes)
    return (celltype(nodes),)
end
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
    linewidthslider = Slider(fig, range = 0:1:10, startvalue=1)
    markersize = lift(x->x,markerslider.value)
    linewidth = lift(x->x,linewidthslider.value)

    #plot the fe-mesh
    wireframep = wireframe!(plotter,markersize=markersize,linewidth=linewidth,deformation_field= @lift $(toggles[2].active) ? $(deformation_field) : :default)
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
                      Label(fig,"linewidth",width=nothing), linewidthslider,
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

####### One Shot Methods #######
const FerriteVizPlots = Union{Type{<:Wireframe},Type{<:SolutionPlot},Type{<:Arrows},Type{<:Surface}}
# We default with our axis choice to the spatial dimension of the problem
function Makie.args_preferred_axis(a, b::Union{MakiePlotter{sdim},Grid{sdim}}, args...) where {sdim}
    if sdim ≤ 2
        return Makie.Axis
    else
        return Makie.LScene
    end
end
function Makie.args_preferred_axis(a::Type{<:Elementinfo}, ip_or_cell)
    dim = Ferrite.getrefdim(ip_or_cell)
    if dim ≤ 2
        return Makie.Axis
    else
        return Makie.LScene
    end
end
# Surface plots are special, as they are 2D problems which are deformed into the third dimension
Makie.args_preferred_axis(a::Type{<:Surface}, b::Union{MakiePlotter{sdim},Grid{sdim}}, args...) where {sdim} = Makie.LScene

function Makie.convert_arguments(P::FerriteVizPlots, dh::Ferrite.AbstractDofHandler, u::AbstractVector)
    return (MakiePlotter(dh,u),)
end
