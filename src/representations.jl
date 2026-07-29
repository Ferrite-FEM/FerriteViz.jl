# Layer 3: thin representations (Makie recipes).
#
# Representations take an FEData and a named data array to color by; field
# selection, deformation and data processing all happen upstream in filters.

function base_fe_attributes()
    return Makie.Attributes(
        shading=Makie.NoShading,
        colormap=:cividis,
        colorrange=Makie.automatic,
        nan_color=:red,
    )
end

# Backend shim for CairoMakie (Petur Bryde, #146, fixes #118).
#
# The pipeline shares its coordinates and triangles into `ShaderAbstractions.Buffer`s
# so GL/WGLMakie can mutate the GPU data in place on `update!`. CairoMakie's
# software mesh path, however, expects plain `Vector`s and does not accept a
# `Buffer` for the faces. When CairoMakie is the active backend we therefore
# draw from the buffers' underlying vectors instead of the buffer-backed
# `GeometryBasics.Mesh` (CairoMakie renders a static frame, so losing the live
# buffer link is inconsequential — `data()` still returns the vector kept in
# sync with the coordinate observable).
function _is_cairomakie_backend()
    backend = Makie.current_backend()
    return !ismissing(backend) && nameof(backend) === :CairoMakie
end

_buffer_data(x) = x
_buffer_data(x::ShaderAbstractions.Buffer) = ShaderAbstractions.data(x)

function _mesh!(parent, ds::FEData; kwargs...)
    if _is_cairomakie_backend()
        return Makie.mesh!(parent, _buffer_data(ds.coords_buffer), _buffer_data(ds.vis_triangles); kwargs...)
    else
        return Makie.mesh!(parent, ds.mesh; kwargs...)
    end
end

function _mesh!(parent, vertices, faces; kwargs...)
    if _is_cairomakie_backend()
        return Makie.mesh!(parent, vertices, _buffer_data(faces); kwargs...)
    else
        return Makie.mesh!(parent, vertices, faces; kwargs...)
    end
end

# Resolve a recipe's color attribute: a Symbol naming a field / point-data /
# cell-data array (or :default) resolves to the scalar per-vertex array,
# anything else passes through as a plain Makie color. Handles dynamic
# switching (the attribute changing to another name) by rewiring the inner
# listener; `colorattr` may be an Observable or a recipe-attribute Computed.
# All listeners are registered to `plot` for cleanup on plot deletion.
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

- `color=:default`: name of the field / point-data / cell-data array to color
  by. `:default` is the first field of the dof handler, reduced to its
  magnitude if vector-valued; an explicitly named array must be scalar —
  reduce with e.g. `Magnitude()` first. Anything that is not a data name is
  passed through to Makie as a plain color.
- `colormap=:cividis`
- `colorrange`: (min, max) of the colorscale, automatic by default.
- `shading=Makie.NoShading`
- `nan_color=:red`

Deformation is an upstream concern: `solutionplot(ds |> WarpByVector(:u, 2.0))`.
"""
@recipe(SolutionPlot) do scene
    attrs = base_fe_attributes()
    attrs[:color] = :default
    attrs
end

function Makie.plot!(SP::SolutionPlot{<:Tuple{<:FEData}})
    ds = SP[1][]
    solution = resolve_color(SP, ds, SP[:color])
    return _mesh!(SP, ds, color=solution, shading=SP[:shading], colormap=SP[:colormap],
                  colorrange=SP[:colorrange], nan_color=SP[:nan_color])
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
@recipe(CellPlot) do scene
    attrs = base_fe_attributes()
    attrs[:color] = :default
    attrs
end

function Makie.plot!(CP::CellPlot{<:Tuple{<:FEData,<:AbstractVector}})
    ds = CP[1][]
    solution = Makie.lift(v -> transfer_scalar_celldata(ds, v), CP[2])
    return _mesh!(CP, ds, color=solution, shading=CP[:shading], colormap=CP[:colormap],
                  colorrange=CP[:colorrange], nan_color=CP[:nan_color])
end

function Makie.plot!(CP::CellPlot{<:Tuple{<:FEData}})
    ds = CP[1][]
    solution = resolve_color(CP, ds, CP[:color])
    return _mesh!(CP, ds, color=solution, shading=CP[:shading], colormap=CP[:colormap],
                  colorrange=CP[:colorrange], nan_color=CP[:nan_color])
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
according to the dataset's `edge_resolution` (see [`FEData`](@ref)).

Node markers and labels are drawn at the grid nodes of the visible cells; they
are displaced by a warp only when the warp field is a dof field.

- `plotnodes=true` plot the nodes as circles/spheres
- `linewidth` edge line width
- `color` color of edges and nodes
- `markersize` size of the nodes
- `cellsets=false` color cells by their cellset association
- `nodelabels=false` global node id labels
- `nodelabelcolor=:darkblue`
- `celllabels=false` global cell id labels
- `celllabelcolor=:darkred`
- `fontsize=15` label text size
- `visible=true`
"""
@recipe(MeshPlot) do scene
    Makie.Attributes(
        plotnodes=true,
        color=theme(scene, :linecolor),
        linewidth=theme(scene, :linewidth),
        markersize=theme(scene, :markersize),
        visible=true,
        fontsize=15,
        offset=(0.0, 0.0),
        nodelabels=false,
        nodelabelcolor=:darkblue,
        celllabels=false,
        celllabelcolor=:darkred,
        cellsets=false,
        depth_shift=-0.0001f0,
    )
end

function Makie.plot!(WF::MeshPlot{<:Tuple{<:FEData{dim}}}) where {dim}
    ds = WF[1][]
    grid = Ferrite.get_grid(ds.dh)
    # Makie only draws 2D/3D points; pad 1D grids with a zero y-coordinate
    pointtype = GeometryBasics.Point{max(dim, 2),Float32}
    topoint(c) = dim == 1 ? pointtype(c[1], 0) : pointtype(c...)
    # The wireframe is a gather from the tessellation coordinates over the
    # visible cells' edge segments. The index list is static; every dynamic
    # concern (deformation, curved geometry, refinement) is already baked into
    # ds.coords by the upstream pipeline, and clipping into ds.visible.
    edge_indices = _visible_edge_indices(ds)
    lines = Makie.lift(cs -> [topoint(cs[i]) for i in edge_indices], ds.coords)
    # cellset coloring
    cellset_u = cellset_data(grid)
    colorrange = (0, max(1, isempty(cellset_u) ? 1 : maximum(cellset_u)))
    _mesh!(WF, ds, color=transfer_scalar_celldata(ds, cellset_u), shading=Makie.NoShading,
           colormap=:darktest, colorrange=colorrange, visible=WF[:cellsets])
    # nodes (of the visible cells)
    visible_nodes = _visible_node_ids(ds)
    gridnodes = Makie.lift(ns -> [topoint(ns[i]) for i in visible_nodes], ds.gridnodes)
    shouldplot = @lift($(WF[:visible]) && $(WF[:plotnodes]))
    Makie.scatter!(WF, gridnodes, markersize=WF[:markersize], color=WF[:color], visible=shouldplot)
    # labels (global ids, restricted to the visible cells and their nodes)
    visible_cells = findall(ds.visible)
    nodelabels = @lift $(WF[:nodelabels]) ? ["$i" for i in visible_nodes] : [""]
    nodepositions = @lift $(WF[:nodelabels]) ? $gridnodes : pointtype[zero(pointtype)]
    celllabels = @lift $(WF[:celllabels]) ? ["$i" for i in visible_cells] : [""]
    cellpositions = @lift $(WF[:celllabels]) ?
                    [topoint(midpoint(Ferrite.getcells(grid, i), $(ds.gridnodes))) for i in visible_cells] :
                    [zero(pointtype)]
    Makie.text!(WF, nodepositions, text=nodelabels, fontsize=WF[:fontsize], offset=WF[:offset], color=WF[:nodelabelcolor])
    Makie.text!(WF, cellpositions, text=celllabels, fontsize=WF[:fontsize], color=WF[:celllabelcolor], align=(:center, :center))
    # edges (3D) / faces (2D) of the mesh
    return Makie.linesegments!(WF, lines, color=WF[:color], linewidth=WF[:linewidth], visible=WF[:visible], depth_shift=WF[:depth_shift])
end

Makie.convert_arguments(::Type{<:MeshPlot}, grid::Ferrite.AbstractGrid) = (FEData(Ferrite.DofHandler(grid), Float64[]),)

"""
    surfaceplot(ds::FEData{2}; kwargs...)
    surfaceplot!(ds::FEData{2}; kwargs...)

Plot a scalar data array of a 2D problem as a surface, with the value as the
z-coordinate. `color=:default` names the array (same resolution rules as
[`solutionplot`](@ref), but it must be a data array).
"""
@recipe(SurfacePlot) do scene
    attrs = base_fe_attributes()
    attrs[:color] = :default
    attrs
end

function Makie.plot!(SF::SurfacePlot{<:Tuple{<:FEData{2}}})
    ds = SF[1][]
    solution = resolve_color(SF, ds, SF[:color])
    solution[] isa AbstractVector || error("surfaceplot needs a data array as `color`, got $(solution[])")
    positions = Makie.lift(ds.coords, solution) do coords, sol
        [Point3f(coords[i][1], coords[i][2], sol[i]) for i in eachindex(coords)]
    end
    return _mesh!(SF, positions, ds.vis_triangles, color=solution, shading=SF[:shading],
                  colormap=SF[:colormap], colorrange=SF[:colorrange], nan_color=SF[:nan_color])
end

"""
    arrowplot(ds::FEData; kwargs...)
    arrowplot!(ds::FEData; kwargs...)

Draw an arrow at every tessellation vertex for a vector-valued data array
(only for spatial dim ≥ 2).

- `field=:default` name of the vector data array
- `color=:default` scalar data array name to color by, or a plain color;
  `:default` colors by the vector magnitude
- `normalize=false` normalize arrow lengths
- `lengthscale=1f0` scale arrow lengths
- `colormap=:cividis`
"""
@recipe(ArrowPlot) do scene
    Makie.Attributes(
        field=:default,
        color=:default,
        colormap=:cividis,
        normalize=false,
        lengthscale=1f0,
    )
end

function Makie.plot!(AR::ArrowPlot{<:Tuple{<:FEData{dim}}}) where {dim}
    dim >= 2 || error("arrowplot is only available for spatial dim ≥ 2")
    ds = AR[1][]
    fname = Makie.lift(f -> _resolve_name(ds, f), AR[:field])
    vecdata = _switching_point_data(ds, fname; owner=AR)
    directions = Makie.lift(vecdata) do A
        size(A, 2) == dim || error("arrowplot needs a $dim-component vector array, :$(fname[]) has $(size(A, 2))")
        [Makie.Vec{dim,Float32}(view(A, i, :)...) for i in 1:size(A, 1)]
    end
    arrowcolor = Makie.Observable{Any}()
    listener = Ref{Any}(nothing)
    cache = Dict{Symbol,Makie.Observable}()
    magnitude = Makie.lift(d -> LinearAlgebra.norm.(d), directions)
    function connect_color(c)
        if listener[] !== nothing
            Makie.Observables.off(listener[])
            listener[] = nothing
        end
        if c === :default
            listener[] = _register_listener!(AR, Makie.on(v -> arrowcolor[] = v, magnitude))
            arrowcolor[] = magnitude[]
        elseif c isa Symbol && _data_association(ds, c) !== :none
            inner = get!(() -> _scalar_data(ds, c), cache, c)
            listener[] = _register_listener!(AR, Makie.on(v -> arrowcolor[] = v, inner))
            arrowcolor[] = inner[]
        else
            arrowcolor[] = c
        end
    end
    _register_listener!(AR, Makie.on(connect_color, AR[:color]))
    connect_color(AR[:color][])
    arrows! = dim == 2 ? Makie.arrows2d! : Makie.arrows3d!
    return arrows!(AR, ds.coords, directions, color=arrowcolor, colormap=AR[:colormap],
                   normalize=AR[:normalize], lengthscale=AR[:lengthscale])
end

"""
    elementinfo(ip::Interpolation; kwargs...)
    elementinfo(cell::AbstractCell; kwargs...)
    elementinfo(ip::Type{Interpolation}; kwargs...)
    elementinfo(cell::Type{AbstractCell}; kwargs...)

Plot the reference element with vertex/edge/face annotations; for a cell the
geometry nodes are labeled "N", for an interpolation the dofs are labeled "D".

- `plotnodes=true` plot the nodes
- `linewidth` stroke width of faces/edges
- `color`
- `markersize` size of the nodes
- `fontsize=60` fontsize of the labels
- `nodelabels=true`, `nodelabelcolor=:darkred`, `nodelabeloffset=(0.0,20.0)`
- `vertexlabels=true`, `vertexlabelcolor=:darkred`, `vertexlabeloffset=(0.0,0.0)`
- `edgelabels=true`, `edgelabelcolor=:darkblue`, `edgelabeloffset=(-40,-40)`
- `facelabels=true`, `facelabelcolor=:darkgreen`, `facelabeloffset=(-40,0)`
- `font="Julia Mono"`
"""
@recipe(Elementinfo) do scene
    Makie.Attributes(
        plotnodes=true,
        linewidth=theme(scene, :linewidth),
        color=theme(scene, :linecolor),
        markersize=theme(scene, :markersize),
        fontsize=60,
        vertexlabels=true,
        vertexlabelcolor=:darkred,
        vertexlabeloffset=(0.0, 0.0),
        nodelabels=true,
        nodelabelcolor=:darkred,
        nodelabeloffset=(0.0, 20.0),
        facelabels=true,
        facelabelcolor=:darkgreen,
        facelabeloffset=(-40, 0),
        edgelabels=true,
        edgelabelcolor=:darkblue,
        edgelabeloffset=(-40, -40),
        font=theme(scene, :font),
    )
end

function Makie.plot!(Ele::Elementinfo{<:Tuple{<:Ferrite.AbstractCell{refshape}}}) where {refshape}
    cell = Ele[1][]
    gip = Ferrite.geometric_interpolation(typeof(cell))
    _draw_reference_element!(Ele, gip, gip, refshape, "N")
end

function Makie.plot!(Ele::Elementinfo{<:Tuple{<:Ferrite.Interpolation{refshape}}}) where {refshape}
    ip = Ele[1][]
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
    Makie.linesegments!(Ele, lines, color=Ele[:color], linewidth=Ele[:linewidth])

    # nodes (geometry nodes or dofs) and their labels
    elenodes = [topoint(ξ...) for ξ in Ferrite.reference_coordinates(node_ip)]
    Makie.scatter!(Ele, elenodes, markersize=Ele[:markersize], color=Ele[:color], visible=Ele[:plotnodes])
    nodelabels = @lift $(Ele[:nodelabels]) ? ["$nodeprefix$i" for i in 1:length(elenodes)] : [""]
    nodepositions = @lift $(Ele[:nodelabels]) ? elenodes : [zero(topoint)]
    Makie.text!(Ele, nodepositions, text=nodelabels, fontsize=Ele[:fontsize], offset=Ele[:nodelabeloffset],
                color=Ele[:nodelabelcolor], font=Ele[:font])

    # vertex annotations
    if Ele[:vertexlabels][]
        for (id, vertexnodes) in enumerate(Ferrite.vertexdof_indices(gip))
            Makie.text!(Ele, "V$id", position=geocoords[vertexnodes[1]], fontsize=Ele[:fontsize],
                        offset=Ele[:vertexlabeloffset], color=Ele[:vertexlabelcolor], font=Ele[:font])
        end
    end
    # edge annotations
    if dim ≥ 2 && Ele[:edgelabels][]
        for (id, edgenodes) in enumerate(Ferrite.edgedof_indices(gip))
            position = (geocoords[edgenodes[1]] + geocoords[edgenodes[2]]) * 0.5
            Makie.text!(Ele, "E$id", position=position, fontsize=Ele[:fontsize], offset=Ele[:edgelabeloffset],
                        color=Ele[:edgelabelcolor], font=Ele[:font])
        end
    end
    # face annotations
    if dim ≥ 3 && Ele[:facelabels][]
        for (id, face) in enumerate(Ferrite.reference_faces(refshape))
            position = sum(geocoords[collect(face)]) / length(face)
            Makie.text!(Ele, "F$id", position=position, fontsize=Ele[:fontsize], offset=Ele[:facelabeloffset],
                        color=Ele[:facelabelcolor], font=Ele[:font])
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
