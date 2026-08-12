# Layer 3: thin representations (Makie recipes).
#
# Representations take an FEData and a named data array to color by; field
# selection, deformation and data processing all happen upstream in filters.
#
# The recipes are new-style (declared attribute blocks) and compute derived
# values in the plot's ComputeGraph (`plot.attributes`): FEData's Observables
# enter the graph via `ComputePipeline.add_input!`, transformations are `map!`
# edges, and child plots draw from the resulting graph nodes. The one thing
# that stays Observable-side is following a *named* data array when the name
# attribute changes (see `resolve_color`): which array a plot listens to is a
# structural change, and the graph's dependencies are fixed at registration.

# Attributes shared by the data-colored representations. Colormap-mixin keys
# whose FerriteViz default differs are excluded and re-declared.
function base_fe_attributes()
    return Makie.@DocumentedAttributes begin
        """
        Name of the field / point-data / cell-data array to color by. `:default`
        is the first field of the dof handler, reduced to its magnitude if
        vector-valued; an explicitly named array must be scalar — reduce with
        e.g. `Magnitude()` first. Anything that is not a data name is passed
        through to Makie as a plain color.
        """
        color = :default
        "Sets the colormap that is sampled for numeric colors."
        colormap = :cividis
        "Replacement color for NaN values."
        nan_color = :red
        "Controls if the plot object is shaded by the parent scene's lights."
        shading = Makie.NoShading
        Makie.mixin_colormap_attributes(exclude = (:colormap, :nan_color))...
        Makie.mixin_generic_plot_attributes()...
    end
end

# Backend shim, generalizing the CairoMakie one (Petur Bryde, #146, fixes #118).
#
# The pipeline shares its coordinates and triangles into `ShaderAbstractions.Buffer`s
# so GLMakie can mutate the GPU data in place on `update!` — and only GLMakie:
# CairoMakie's software mesh path expects plain `Vector`s and does not accept a
# `Buffer` for the faces, and WGLMakie renders a `Buffer`-backed mesh once but
# never applies subsequent `Buffer` updates — the browser silently keeps the
# initial geometry (verified on WGLMakie 0.13.13). On every backend but
# GLMakie we therefore draw from the coordinate *Observable* the buffer wraps
# (`ds.coords`) plus the triangles' underlying vector: coordinate and color
# updates then flow through Makie's ordinary observable path (live WGLMakie
# viewers keep updating over the websocket); only CrinkleClip-style *triangle*
# changes need the buffer link and won't propagate outside GLMakie.
function _is_glmakie_backend()
    backend = Makie.current_backend()
    return !ismissing(backend) && nameof(backend) === :GLMakie
end

_buffer_data(x) = x
_buffer_data(x::ShaderAbstractions.Buffer) = ShaderAbstractions.data(x)

# The parent's attributes are forwarded wholesale (the child keeps only the
# ones it understands); explicit keyword arguments override forwarded ones.
function _mesh!(parent, ds::FEData; kwargs...)
    if _is_glmakie_backend()
        return Makie.mesh!(parent, parent.attributes, ds.mesh; kwargs...)
    else
        return Makie.mesh!(parent, parent.attributes, ds.coords, _buffer_data(ds.vis_triangles); kwargs...)
    end
end

function _mesh!(parent, vertices, faces; kwargs...)
    if _is_glmakie_backend()
        return Makie.mesh!(parent, parent.attributes, vertices, faces; kwargs...)
    else
        return Makie.mesh!(parent, parent.attributes, vertices, _buffer_data(faces); kwargs...)
    end
end

# Resolve a recipe's color attribute: a Symbol naming a field / point-data /
# cell-data array (or :default) resolves to the scalar per-vertex array,
# anything else passes through as a plain Makie color. Handles dynamic
# switching (the attribute changing to another name) by rewiring the inner
# listener; `colorattr` is the recipe's :color node (a Computed; a plain
# Observable works too). All listeners are registered to `plot` for cleanup
# on plot deletion.
function resolve_color(plot, ds::FEData, colorattr)
    out = Makie.Observable{Any}()
    listener = Ref{Any}(nothing)
    cache = Dict{Symbol,Makie.Observable}() # avoid re-lifting on repeated switches
    function connect(val)
        if listener[] !== nothing
            Makie.Observables.off(listener[]) # double-off at plot deletion is harmless
            listener[] = nothing
        end
        if val isa Symbol && (val === :default || _data_association(ds, val) !== :none)
            # :default falls back to the magnitude for vector-valued fields
            inner = get!(() -> _scalar_data(ds, val; reduce_default=val === :default), cache, val)
            listener[] = _register_listener!(plot, Makie.on(v -> out[] = v, inner))
            out[] = inner[]
        else
            out[] = val
        end
    end
    _register_listener!(plot, Makie.on(connect, colorattr))
    connect(colorattr[])
    return out
end

# Conversion function making an input node untyped: FEData-driven inputs can
# change type (a named data array one moment, a plain color the next), and an
# input's storage is otherwise locked to the type of its first value. Same
# mechanism as Makie's own polymorphic attributes (e.g. colormap).
_untyped_input(_key, value) = Base.RefValue{Any}(value)

# The resolved color enters the plot's compute graph as the :color_data input,
# so downstream consumers are graph nodes.
function _graph_color!(plot, ds::FEData)
    ComputePipeline.add_input!(_untyped_input, plot.attributes, :color_data, resolve_color(plot, ds, plot.color))
    return plot.color_data
end

function cellset_data(grid::Ferrite.AbstractGrid)
    data = zeros(Float64, Ferrite.getncells(grid))
    for (setidx, (_, cellset)) in enumerate(grid.cellsets)
        for cellidx in cellset
            data[cellidx] = setidx
        end
    end
    return data
end

"""
    solutionplot(ds::FEData; kwargs...)
    solutionplot(dh::AbstractDofHandler, u::Vector; kwargs...)
    solutionplot!(...)

Contour plot of a scalar data array on the finite element mesh.

Deformation is an upstream concern: `solutionplot(ds |> WarpByVector(:u, 2.0))`.
"""
Makie.@recipe SolutionPlot (dataset,) begin
    base_fe_attributes()...
end

function Makie.plot!(SP::SolutionPlot{<:Tuple{<:FEData}})
    ds = SP.dataset[]
    _mesh!(SP, ds, color=_graph_color!(SP, ds))
    return SP
end

"""
    cellplot(ds::FEData, values::Vector{<:Real}; kwargs...)
    cellplot(ds::FEData; color=:name, kwargs...)
    cellplot!(...)

Plot one scalar per cell as constant color on the cells, either passed
directly as a vector or by naming a registered cell-data array (see
[`set_cell_data!`](@ref)). Non-scalar per-cell data (e.g. stress tensors) is
reduced with a filter first (e.g. [`VonMises`](@ref)). Shares the
solutionplot kwargs.
"""
Makie.@recipe CellPlot (dataset, values) begin
    base_fe_attributes()...
end

# The converted arguments must always match the declared arity (Makie
# destructures them into one node per declared name), so the values-less form
# is normalized to an empty sentinel vector.
Makie.convert_arguments(::Type{<:CellPlot}, ds::FEData) = (ds, Float64[])
Makie.convert_arguments(::Type{<:CellPlot}, ds::FEData, values::AbstractVector) =
    (ds, convert(Vector{Float64}, values))

function Makie.plot!(CP::CellPlot{<:Tuple{<:FEData,<:AbstractVector}})
    ds = CP.dataset[]
    if isempty(CP.values[])
        # no per-cell values passed: color by the named cell-data array
        _mesh!(CP, ds, color=_graph_color!(CP, ds))
    else
        Makie.map!(v -> transfer_scalar_celldata(ds, v), CP.attributes, :values, :vertex_color)
        _mesh!(CP, ds, color=CP.vertex_color)
    end
    return CP
end

"""
    meshplot(ds::FEData; kwargs...)
    meshplot(grid::AbstractGrid; kwargs...)
    meshplot!(...)

Plot the finite element mesh (edges and nodes), optionally labeled. The
wireframe is drawn from the dataset's tessellation edges, i.e. from the same
vertices as the surface plots: it follows an upstream [`WarpByVector`](@ref)
(including high-order and discontinuous deformation), is hidden with the cells
a [`CrinkleClip`](@ref) removes, and bends along curved (high-order) cell edges
according to the dataset's subdivision (see [`Refine`](@ref)).

Node markers and labels are drawn at the grid nodes of the visible cells; they
are displaced by a warp only when the warp field is a dof field.
"""
Makie.@recipe MeshPlot (dataset,) begin
    "Plot the nodes as circles/spheres."
    plotnodes = true
    "Color of edges and nodes."
    color = @inherit linecolor
    "Edge line width."
    linewidth = @inherit linewidth
    "Size of the node markers."
    markersize = @inherit markersize
    "Label text size."
    fontsize = 15
    "Offset of the node labels."
    offset = (0.0, 0.0)
    "Show global node id labels."
    nodelabels = false
    "Color of the node id labels."
    nodelabelcolor = :darkblue
    "Show global cell id labels."
    celllabels = false
    "Color of the cell id labels."
    celllabelcolor = :darkred
    "Color cells by their cellset association."
    cellsets = false
    Makie.filter_attributes(Makie.mixin_generic_plot_attributes(); exclude = (:depth_shift,))...
    "Depth shift drawing the wireframe in front of surface plots."
    depth_shift = -0.0001f0
end

function Makie.plot!(WF::MeshPlot{<:Tuple{<:FEData{dim}}}) where {dim}
    ds = WF.dataset[]
    grid = Ferrite.get_grid(ds.dh)
    # Makie only draws 2D/3D points; pad 1D grids with a zero y-coordinate
    pointtype = GeometryBasics.Point{max(dim, 2),Float32}
    topoint(c) = dim == 1 ? pointtype(c[1], 0) : pointtype(c...)
    graph = WF.attributes
    ComputePipeline.add_input!(graph, :ds_coords, ds.coords)
    ComputePipeline.add_input!(graph, :ds_gridnodes, ds.gridnodes)
    # The wireframe is a gather from the tessellation coordinates over the
    # visible cells' edge segments. The index list is static; every dynamic
    # concern (deformation, curved geometry, refinement) is already baked into
    # ds.coords by the upstream pipeline, and clipping into ds.visible.
    edge_indices = _visible_edge_indices(ds)
    Makie.map!(cs -> [topoint(cs[i]) for i in edge_indices], graph, :ds_coords, :edge_lines)
    # cellset coloring (depth_shift is meant for the wireframe, not this mesh)
    cellset_u = cellset_data(grid)
    colorrange = (0, max(1, isempty(cellset_u) ? 1 : maximum(cellset_u)))
    _mesh!(WF, ds, color=transfer_scalar_celldata(ds, cellset_u), shading=Makie.NoShading,
           colormap=:darktest, colorrange=colorrange, visible=WF.cellsets, depth_shift=0.0f0)
    # nodes (of the visible cells)
    visible_nodes = _visible_node_ids(ds)
    Makie.map!(ns -> [topoint(ns[i]) for i in visible_nodes], graph, :ds_gridnodes, :node_positions)
    Makie.map!((v, p) -> v && p, graph, [:visible, :plotnodes], :shownodes)
    Makie.scatter!(WF, WF.node_positions, markersize=WF.markersize, color=WF.color, visible=WF.shownodes)
    # labels (global ids, restricted to the visible cells and their nodes)
    visible_cells = findall(ds.visible)
    Makie.map!(graph, [:nodelabels, :node_positions], [:nodelabel_text, :nodelabel_positions]) do nl, positions
        nl ? (["$i" for i in visible_nodes], positions) : ([""], pointtype[zero(pointtype)])
    end
    Makie.map!(graph, [:celllabels, :ds_gridnodes], [:celllabel_text, :celllabel_positions]) do cl, ns
        if cl
            (["$i" for i in visible_cells],
             [topoint(midpoint(Ferrite.getcells(grid, i), ns)) for i in visible_cells])
        else
            ([""], [zero(pointtype)])
        end
    end
    Makie.text!(WF, WF.nodelabel_positions, text=WF.nodelabel_text, fontsize=WF.fontsize,
                offset=WF.offset, color=WF.nodelabelcolor)
    Makie.text!(WF, WF.celllabel_positions, text=WF.celllabel_text, fontsize=WF.fontsize,
                color=WF.celllabelcolor, align=(:center, :center))
    # edges (3D) / faces (2D) of the mesh
    return Makie.linesegments!(WF, WF.edge_lines, color=WF.color, linewidth=WF.linewidth,
                               visible=WF.visible, depth_shift=WF.depth_shift)
end

Makie.convert_arguments(::Type{<:MeshPlot}, grid::Ferrite.AbstractGrid) = (FEData(Ferrite.DofHandler(grid), Float64[]),)

"""
    surfaceplot(ds::FEData{2}; kwargs...)
    surfaceplot!(ds::FEData{2}; kwargs...)

Plot a scalar data array of a 2D problem as a surface, with the value as the
z-coordinate. `color=:default` names the array (same resolution rules as
[`solutionplot`](@ref), but it must be a data array).
"""
Makie.@recipe SurfacePlot (dataset,) begin
    base_fe_attributes()...
end

function Makie.plot!(SF::SurfacePlot{<:Tuple{<:FEData{2}}})
    ds = SF.dataset[]
    color_data = _graph_color!(SF, ds)
    # eager, so construction throws a plain error instead of a wrapped
    # ResolveException from inside the graph edge below
    color_data[] isa AbstractVector || error("surfaceplot needs a data array as `color`, got $(color_data[])")
    ComputePipeline.add_input!(SF.attributes, :ds_coords, ds.coords)
    Makie.map!(SF.attributes, [:ds_coords, :color_data], :positions) do coords, sol
        sol isa AbstractVector || error("surfaceplot needs a data array as `color`, got $sol")
        [Point3f(coords[i][1], coords[i][2], sol[i]) for i in eachindex(coords)]
    end
    return _mesh!(SF, SF.positions, ds.vis_triangles, color=SF.color_data)
end

"""
    arrowplot(ds::FEData; kwargs...)
    arrowplot!(ds::FEData; kwargs...)

Draw an arrow at every tessellation vertex for a vector-valued data array
(only for spatial dim ≥ 2).
"""
Makie.@recipe ArrowPlot (dataset,) begin
    "Name of the vector data array."
    field = :default
    """
    Scalar data array name to color by, or a plain color; `:default` colors
    by the vector magnitude.
    """
    color = :default
    "Sets the colormap that is sampled for numeric colors."
    colormap = :cividis
    "Normalize arrow lengths."
    normalize = false
    "Scale arrow lengths."
    lengthscale = 1.0f0
    Makie.mixin_colormap_attributes(exclude = (:colormap,))...
    Makie.mixin_generic_plot_attributes()...
end

# One Makie.Vec per row of a data matrix, per frame for arrowplot. Function
# barrier: `dim` captured in the map! closure is a plain `Int`, which would
# make `Vec{dim,Float32}` a dynamic type application on every row.
function _row_vectors(A::AbstractMatrix, ::Val{dim}) where {dim}
    return [Makie.Vec{dim,Float32}(ntuple(j -> Float32(A[i, j]), Val(dim))) for i in 1:size(A, 1)]
end

# Follow the named color array like `resolve_color`, but emit `nothing` when
# the color attribute is not a data name (arrowplot resolves :default to the
# arrow magnitude, a graph node, so the fallback lives graph-side).
function _named_color_data(plot, ds::FEData, colorattr)
    out = Makie.Observable{Any}(nothing)
    listener = Ref{Any}(nothing)
    cache = Dict{Symbol,Makie.Observable}()
    function connect(val)
        if listener[] !== nothing
            Makie.Observables.off(listener[]) # double-off at plot deletion is harmless
            listener[] = nothing
        end
        if val isa Symbol && val !== :default && _data_association(ds, val) !== :none
            inner = get!(() -> _scalar_data(ds, val), cache, val)
            listener[] = _register_listener!(plot, Makie.on(v -> out[] = v, inner))
            out[] = inner[]
        else
            out[] = nothing
        end
    end
    _register_listener!(plot, Makie.on(connect, colorattr))
    connect(colorattr[])
    return out
end

function Makie.plot!(AR::ArrowPlot{<:Tuple{<:FEData{dim}}}) where {dim}
    dim >= 2 || error("arrowplot is only available for spatial dim ≥ 2")
    ds = AR.dataset[]
    graph = AR.attributes
    # vector data: the name switching is Observable-side, the resolved array
    # enters the graph
    fname = Makie.lift(f -> _resolve_name(ds, f), ComputePipeline.get_observable!(graph, :field))
    vecdata = _switching_point_data(ds, fname; owner=AR)
    # eager, so construction throws a plain error instead of a wrapped
    # ResolveException from inside the graph edge below
    size(vecdata[], 2) == dim || error("arrowplot needs a $dim-component vector array, :$(fname[]) has $(size(vecdata[], 2))")
    ComputePipeline.add_input!(graph, :vector_data, vecdata)
    Makie.map!(graph, [:vector_data], :directions) do A
        size(A, 2) == dim || error("arrowplot needs a $dim-component vector array, :$(fname[]) has $(size(A, 2))")
        _row_vectors(A, Val(dim))
    end
    Makie.map!(d -> LinearAlgebra.norm.(d), graph, :directions, :magnitude)
    ComputePipeline.add_input!(_untyped_input, graph, :named_color_data, _named_color_data(AR, ds, AR.color))
    Makie.map!(graph, [:color, :named_color_data, :magnitude], :arrow_color) do c, named, mag
        Base.RefValue{Any}(c === :default ? mag : (named === nothing ? c : named))
    end
    arrows! = dim == 2 ? Makie.arrows2d! : Makie.arrows3d!
    return arrows!(AR, AR.attributes, ds.coords, AR.directions, color=AR.arrow_color)
end

"""
    elementinfo(ip::Interpolation; kwargs...)
    elementinfo(cell::AbstractCell; kwargs...)
    elementinfo(ip::Type{Interpolation}; kwargs...)
    elementinfo(cell::Type{AbstractCell}; kwargs...)

Plot the reference element with vertex/edge/face annotations; for a cell the
geometry nodes are labeled "N", for an interpolation the dofs are labeled "D".
"""
Makie.@recipe Elementinfo (element,) begin
    "Plot the nodes."
    plotnodes = true
    "Stroke width of faces/edges."
    linewidth = @inherit linewidth
    "Color of edges and nodes."
    color = @inherit linecolor
    "Size of the node markers."
    markersize = @inherit markersize
    "Fontsize of the labels."
    fontsize = 60
    "Show the node labels."
    nodelabels = true
    "Color of the node labels."
    nodelabelcolor = :darkred
    "Offset of the node labels."
    nodelabeloffset = (0.0, 20.0)
    "Show the vertex labels."
    vertexlabels = true
    "Color of the vertex labels."
    vertexlabelcolor = :darkred
    "Offset of the vertex labels."
    vertexlabeloffset = (0.0, 0.0)
    "Show the edge labels."
    edgelabels = true
    "Color of the edge labels."
    edgelabelcolor = :darkblue
    "Offset of the edge labels."
    edgelabeloffset = (-40, -40)
    "Show the face labels."
    facelabels = true
    "Color of the face labels."
    facelabelcolor = :darkgreen
    "Offset of the face labels."
    facelabeloffset = (-40, 0)
    "Font of the labels."
    font = @inherit font
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(Ele::Elementinfo{<:Tuple{<:Ferrite.AbstractCell{refshape}}}) where {refshape}
    cell = Ele.element[]
    gip = Ferrite.geometric_interpolation(typeof(cell))
    _draw_reference_element!(Ele, gip, gip, refshape, "N")
end

function Makie.plot!(Ele::Elementinfo{<:Tuple{<:Ferrite.Interpolation{refshape}}}) where {refshape}
    ip = Ele.element[]
    gip = Ferrite.default_geometric_interpolation(ip)
    gip isa VectorizedInterpolation && (gip = gip.ip)
    _draw_reference_element!(Ele, ip, gip, refshape, "D")
end

# Shared drawing for both Elementinfo variants: `node_ip` provides the labeled
# node positions, `gip` (scalar) the element geometry/boundary structure.
function _draw_reference_element!(Ele, node_ip::Ferrite.Interpolation, gip::Ferrite.ScalarInterpolation, refshape, nodeprefix::String)
    dim = Ferrite.getrefdim(gip)
    topoint = dim > 2 ? Point3f : Point2f
    geocoords = [topoint(ξ...) for ξ in Ferrite.reference_coordinates(gip)]

    # element boundary
    lines = topoint[]
    for edgenodes in Ferrite.edgedof_indices(gip)
        # by convention the first two edge dofs are the vertices
        push!(lines, geocoords[edgenodes[1]], geocoords[edgenodes[2]])
    end
    Makie.linesegments!(Ele, lines, color=Ele.color, linewidth=Ele.linewidth)

    # nodes (geometry nodes or dofs) and their labels
    elenodes = [topoint(ξ...) for ξ in Ferrite.reference_coordinates(node_ip)]
    Makie.scatter!(Ele, elenodes, markersize=Ele.markersize, color=Ele.color, visible=Ele.plotnodes)
    Makie.map!(Ele.attributes, [:nodelabels], [:nodelabel_text, :nodelabel_positions]) do nl
        nl ? (["$nodeprefix$i" for i in 1:length(elenodes)], elenodes) : ([""], [zero(topoint)])
    end
    Makie.text!(Ele, Ele.nodelabel_positions, text=Ele.nodelabel_text, fontsize=Ele.fontsize,
                offset=Ele.nodelabeloffset, color=Ele.nodelabelcolor, font=Ele.font)

    # vertex annotations
    if Ele.vertexlabels[]
        for (id, vertexnodes) in enumerate(Ferrite.vertexdof_indices(gip))
            Makie.text!(Ele, "V$id", position=geocoords[vertexnodes[1]], fontsize=Ele.fontsize,
                        offset=Ele.vertexlabeloffset, color=Ele.vertexlabelcolor, font=Ele.font)
        end
    end
    # edge annotations
    if dim ≥ 2 && Ele.edgelabels[]
        for (id, edgenodes) in enumerate(Ferrite.edgedof_indices(gip))
            position = (geocoords[edgenodes[1]] + geocoords[edgenodes[2]]) * 0.5
            Makie.text!(Ele, "E$id", position=position, fontsize=Ele.fontsize, offset=Ele.edgelabeloffset,
                        color=Ele.edgelabelcolor, font=Ele.font)
        end
    end
    # face annotations
    if dim ≥ 3 && Ele.facelabels[]
        for (id, face) in enumerate(Ferrite.reference_faces(refshape))
            position = sum(geocoords[collect(face)]) / length(face)
            Makie.text!(Ele, "F$id", position=position, fontsize=Ele.fontsize, offset=Ele.facelabeloffset,
                        color=Ele.facelabelcolor, font=Ele.font)
        end
    end
    return Ele
end

function Makie.convert_arguments(P::Type{<:Elementinfo}, celltype::Type{C}) where {C<:Ferrite.AbstractCell}
    gip = Ferrite.geometric_interpolation(C)
    nnodes = Ferrite.getnbasefunctions(gip)
    return (celltype(ntuple(x -> 1, nnodes)),)
end
Makie.convert_arguments(P::Type{<:Elementinfo}, iptype::Type{IP}) where {IP<:Ferrite.Interpolation} = (iptype(),)

####### One Shot Methods #######
const FerriteVizPlots = Union{Type{<:MeshPlot},Type{<:SolutionPlot},Type{<:ArrowPlot},Type{<:SurfacePlot}}

# We default with our axis choice to the spatial dimension of the problem
function Makie.args_preferred_axis(a, b::Union{FEData{sdim},Grid{sdim}}, args...) where {sdim}
    return sdim ≤ 2 ? Makie.Axis : Makie.LScene
end
function Makie.args_preferred_axis(a::Type{<:Elementinfo}, ip_or_cell)
    return Ferrite.getrefdim(ip_or_cell) ≤ 2 ? Makie.Axis : Makie.LScene
end
# Surface plots are special, as they are 2D problems which are deformed into the third dimension
Makie.args_preferred_axis(a::Type{<:SurfacePlot}, b::Union{FEData{sdim},Grid{sdim}}, args...) where {sdim} = Makie.LScene

function Makie.convert_arguments(P::FerriteVizPlots, dh::Ferrite.AbstractDofHandler, u::AbstractVector)
    return (FEData(dh, collect(u)),)
end
