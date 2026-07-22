# Composable viewer.
#
# A viewer is described declaratively with Makie's SpecApi: a `layout(ds, state)`
# hook returns a `GridLayoutSpec`, driven by an `Observable` of view state that a
# set of pluggable `Control`s feed. Makie diffs successive specs and updates only
# what changed, so switching field/colormap/labels re-uses the existing plots
# (their `resolve_color` already handles reactive color switching).
#
# Two reactivity mechanisms coexist:
#   * structural view state (field, colormap, labels, which panels) rebuilds the
#     spec — cheap, since the same `ds` object is reused and only attributes diff;
#   * data streaming (FerriteViz.update!, deformation scale) mutates the shared
#     GPU buffers in place and never rebuilds the spec.

const S = Makie.SpecApi

##############
# Spec helpers
##############

"""
    solutionplotspec(ds; kwargs...) -> Makie.PlotSpec

`PlotSpec` for a [`solutionplot`](@ref), for use inside a [`layout`](@ref
ferriteviewer) hook / `Makie.SpecApi`. Analogous helpers exist for every
representation: [`meshplotspec`](@ref), [`cellplotspec`](@ref),
[`surfaceplotspec`](@ref), [`arrowplotspec`](@ref).
"""
solutionplotspec(ds; kwargs...) = Makie.PlotSpec(SolutionPlot, ds; kwargs...)
"""    meshplotspec(ds; kwargs...) -> Makie.PlotSpec"""
meshplotspec(ds; kwargs...) = Makie.PlotSpec(MeshPlot, ds; kwargs...)
"""    cellplotspec(ds; kwargs...) -> Makie.PlotSpec"""
cellplotspec(ds; kwargs...) = Makie.PlotSpec(CellPlot, ds; kwargs...)
"""    surfaceplotspec(ds; kwargs...) -> Makie.PlotSpec"""
surfaceplotspec(ds; kwargs...) = Makie.PlotSpec(SurfacePlot, ds; kwargs...)
"""    arrowplotspec(ds; kwargs...) -> Makie.PlotSpec"""
arrowplotspec(ds; kwargs...) = Makie.PlotSpec(ArrowPlot, ds; kwargs...)

"""
    panelspec(plots...; colorbar=nothing, dim=2, axis=(;), colorbar_attributes=(;))

Assemble one or more `PlotSpec`s into a `GridLayoutSpec` of a single axis (an
`Axis` for `dim ≤ 2`, an `LScene` otherwise), optionally with a linked
`Colorbar`. Pass the `PlotSpec` to colorbar-link as `colorbar` (typically one of
`plots`). `axis` is a `NamedTuple` of axis attributes.
"""
function panelspec(plots...; colorbar=nothing, dim::Int=2, axis=(;), colorbar_attributes=(;))
    specs = collect(Makie.PlotSpec, plots)
    axblock = dim > 2 ? S.LScene(; plots=specs, axis...) : S.Axis(; plots=specs, axis...)
    colorbar === nothing && return S.GridLayout([axblock])
    return S.GridLayout([axblock S.Colorbar(colorbar; _colorbar_kw(colorbar, colorbar_attributes)...)])
end

# A `Colorbar` linked to a `PlotSpec` derives its colormap from the plot: when
# the plot spec doesn't carry an explicit `colormap`, Makie reads the recipe's
# default via `lookup_default`, which our `Attributes`-based recipes (whose
# defaults are Observables) don't support. Forward the plot's colormap — or fall
# back to the representation default — so the colorbar is always self-sufficient.
function _colorbar_kw(colorbar, user_attributes)
    kw = Dict{Symbol,Any}(pairs(user_attributes))
    if colorbar isa Makie.PlotSpec
        get!(kw, :colormap, get(colorbar.kwargs, :colormap, :cividis))
        haskey(colorbar.kwargs, :colorrange) && get!(kw, :colorrange, colorbar.kwargs[:colorrange])
    end
    return kw
end

###########
# Controls
###########

"""
    ControlResult(content; structural=[], dynamic=[], placement=:column)

Return value of a [`Control`](@ref)'s builder. `content` is a vector of
layoutables; `placement` decides where they go: `:column` (the right-hand
controls column, default) or `:below` (each on its own full-width row under the
plot — e.g. a [`TimeSlider`](@ref)). `structural` state pairs
(`name => Observable`) trigger a spec rebuild when they change; `dynamic` state
pairs are made available to the pipeline/layout but do not (they stream through
the shared observables instead, e.g. a deformation scale).
"""
struct ControlResult
    content::Vector{Any}
    structural::Vector{Pair{Symbol,Makie.Observable}}
    dynamic::Vector{Pair{Symbol,Makie.Observable}}
    placement::Symbol
end
function ControlResult(content::Vector; structural=Pair{Symbol,Makie.Observable}[],
                       dynamic=Pair{Symbol,Makie.Observable}[], placement::Symbol=:column)
    placement in (:column, :below) || error("placement must be :column or :below, got :$placement")
    return ControlResult(convert(Vector{Any}, content), structural, dynamic, placement)
end

"""
    Control(make)

A pluggable viewer control. `make(fig, ds)` builds its widget(s) and returns a
[`ControlResult`](@ref) declaring the widgets and the view-state observables they
contribute. Built-in controls: [`FieldMenu`](@ref), [`ProcessMenu`](@ref),
[`ColormapMenu`](@ref), [`LabelsToggle`](@ref), [`DeformationToggle`](@ref),
[`TimeSlider`](@ref).
"""
struct Control
    make::Function
end

_labeled(fig, label, w) = Any[Label(fig, label, width=nothing), w]
_toggle_row(fig, label, tog) = Any[grid!([tog Label(fig, label, halign=:left)], tellheight=false)]

"""    FieldMenu(; label="field")

Control selecting which dof field to color by (state `:field`)."""
FieldMenu(; label="field") = Control() do fig, ds
    fields = collect(Ferrite.getfieldnames(ds.dh))
    menu = Menu(fig, options=fields)  # Menu defaults to the first option
    ControlResult(_labeled(fig, label, menu); structural=[:field => menu.selection])
end

"""    ProcessMenu(; label="processing", options=["magnitude","x₁","x₂","x₃"])

Control choosing the scalar reduction of the colored field (state `:process`)."""
ProcessMenu(; label="processing", options=["magnitude", "x₁", "x₂", "x₃"]) = Control() do fig, ds
    menu = Menu(fig, options=options)
    ControlResult(_labeled(fig, label, menu); structural=[:process => menu.selection])
end

"""    ColormapMenu(; label="colormap", options=["cividis","inferno","thermal"])

Control choosing the colormap (state `:colormap`, a `Symbol`)."""
ColormapMenu(; label="colormap", options=["cividis", "inferno", "thermal"]) = Control() do fig, ds
    menu = Menu(fig, options=options, direction=:up)
    ControlResult(_labeled(fig, label, menu); structural=[:colormap => Makie.lift(Symbol, menu.selection)])
end

"""    WireframeToggle(; label="wireframe", active=true)

Control toggling the mesh wireframe overlay (state `:wireframe`)."""
WireframeToggle(; label="wireframe", active=true) = Control() do fig, ds
    tog = Toggle(fig, active=active)
    ControlResult(_toggle_row(fig, label, tog); structural=[:wireframe => tog.active])
end

"""    LabelsToggle(; label="labels", active=false)

Control toggling node/cell id labels (state `:labels`)."""
LabelsToggle(; label="labels", active=false) = Control() do fig, ds
    tog = Toggle(fig, active=active)
    ControlResult(_toggle_row(fig, label, tog); structural=[:labels => tog.active])
end

"""    DeformationToggle(; label="deformation", active=false)

Control toggling warp-by-vector deformation (dynamic state `:deform_scale`; the
scale streams through the shared observables without a spec rebuild)."""
DeformationToggle(; label="deformation", active=false) = Control() do fig, ds
    tog = Toggle(fig, active=active)
    scale = Makie.lift(a -> a ? 1.0 : 0.0, tog.active)
    ControlResult(_toggle_row(fig, label, tog); dynamic=[:deform_scale => scale])
end

"""    TimeSlider(u_history; label="timestep")

Control stepping [`FerriteViz.update!`](@ref) through a solution history. Placed
as a full-width slider below the plot (`placement=:below`). Purely
side-effecting: it streams new solutions into the pipeline, contributing no
structural state."""
function TimeSlider(data::AbstractVector{<:AbstractVector}; label="timestep")
    Control() do fig, ds
        sg = SliderGrid(fig, (label=label, range=1:length(data), format=x -> "$x"))
        Makie.on(i -> update!(ds, data[i]), sg.sliders[1].value)
        ControlResult(Any[sg]; placement=:below)
    end
end

###############
# Default hooks
###############

_has_deformable(ds::FEData{dim}) where {dim} =
    any(f -> Ferrite.n_components(ds.dh, f) == dim, Ferrite.getfieldnames(ds.dh))

"""
    default_controls(ds) -> Vector{Control}

The default control set: field/process/colormap menus, a wireframe toggle, a
labels toggle and (if a field can deform the mesh) a deformation toggle. Override
via the `controls` keyword of [`ferriteviewer`](@ref).
"""
function default_controls(ds::FEData)
    ctrls = Control[FieldMenu(), ProcessMenu(), ColormapMenu(), WireframeToggle(), LabelsToggle()]
    _has_deformable(ds) && push!(ctrls, DeformationToggle())
    return ctrls
end

"""
    default_pipeline(ds, state) -> FEData

Build the dataset the default layout plots. Applies an observable-scale
[`WarpByVector`](@ref) when a [`DeformationToggle`](@ref) contributed a
`:deform_scale` (so the toggle deforms the mesh in place).
"""
function default_pipeline(ds::FEData{dim}, state) where {dim}
    haskey(state, :deform_scale) || return ds
    deformable = filter(f -> Ferrite.n_components(ds.dh, f) == dim, collect(Ferrite.getfieldnames(ds.dh)))
    isempty(deformable) && return ds
    return apply(WarpByVector(Makie.Observable(first(deformable)), state[:deform_scale]), ds)
end

"""
    default_layout(ds, state) -> Makie.GridLayoutSpec

The default single-panel layout: a [`solutionplot`](@ref) colored by the selected
field/process, a [`meshplot`](@ref) wireframe (shown/labeled per the wireframe and
labels toggles), and a linked colorbar. Pass your own `layout(ds, state)` to
[`ferriteviewer`](@ref) for any composition of panels and representations.
"""
function default_layout(ds::FEData{dim}, state) where {dim}
    field = hasproperty(state, :field) ? state.field : default_field(ds)
    process = hasproperty(state, :process) ? state.process : "magnitude"
    colormap = hasproperty(state, :colormap) ? state.colormap : :cividis
    labels = hasproperty(state, :labels) ? state.labels : false
    wireframe = hasproperty(state, :wireframe) ? state.wireframe : true
    color = _viewer_color_array!(ds, field, process)
    sol = solutionplotspec(ds; color=color, colormap=colormap)
    msh = meshplotspec(ds; visible=wireframe, nodelabels=labels, celllabels=labels)
    return panelspec(sol, msh; colorbar=sol, dim=dim)
end

# Register (once) the derived scalar array for a field/process combination and
# return its name.
function _viewer_color_array!(ds::FEData{dim}, field::Symbol, process::AbstractString) where {dim}
    ncomps = size(point_data(ds, field)[], 2)
    process == "magnitude" && ncomps == 1 && return field
    valfun = process == "magnitude" ? LinearAlgebra.norm :
             process == "x₁" ? (x -> x[1]) :
             process == "x₂" ? (x -> x[min(2, end)]) :
             (x -> x[min(3, end)])
    name = Symbol(field, :_, process)
    if !haskey(ds.point_data, name)
        ds.point_data[name] = Makie.lift(A -> _rows_to_matrix(valfun, A, dim), point_data(ds, field))
    end
    return name
end

#########
# Scaffold
#########

"""
    ferriteviewer(ds::FEData; layout=default_layout, controls=default_controls(ds), pipeline=default_pipeline)
    ferriteviewer(ds::FEData, u_history::Vector{<:Vector}; kwargs...)

Interactive viewer composed declaratively with `Makie.SpecApi`. `controls` is a
vector of [`Control`](@ref)s feeding a view-state observable; `pipeline(ds,
state)` derives the dataset to plot (default applies deformation); `layout(ds,
state)` returns the `GridLayoutSpec` rendered on every structural state change.
All three are overridable — the defaults reproduce a single solutionplot panel
with a colorbar. The second form appends a [`TimeSlider`](@ref) for a solution
history.

Because plots pass the same `ds` and named data, structural changes reuse the
existing plots (only attributes diff), while [`FerriteViz.update!`](@ref) and
deformation stream through the shared GPU buffers without rebuilding.
"""
function ferriteviewer(ds::FEData;
                       layout=default_layout,
                       controls::AbstractVector{Control}=default_controls(ds),
                       pipeline=default_pipeline)
    fig = Figure()

    column = Any[]           # right-hand controls column
    below = Any[]            # full-width rows under the plot (e.g. a TimeSlider)
    state_pairs = Pair{Symbol,Makie.Observable}[]
    structural_keys = Symbol[]
    for ctrl in controls
        r = ctrl.make(fig, ds)
        append!(r.placement === :below ? below : column, r.content)
        append!(state_pairs, r.structural)
        append!(state_pairs, r.dynamic)
        append!(structural_keys, first.(r.structural))
    end
    state = Dict{Symbol,Makie.Observable}(state_pairs)

    plotds = pipeline(ds, state)

    allkeys = Tuple(keys(state))
    snapshot() = NamedTuple{allkeys}(Tuple(Makie.to_value(state[k]) for k in allkeys))
    view_state = Makie.Observable(snapshot())
    for k in structural_keys
        Makie.on(_ -> view_state[] = snapshot(), state[k])
    end

    spec = Makie.lift(s -> layout(plotds, s), view_state)
    # non-mutating `plot` into an empty slot: GridLayoutSpec is figure-level
    # (`args_preferred_axis == FigureOnly`), so it must create the nested layout
    # itself rather than plot into a pre-made axis (which `plot!` would require).
    Makie.plot(fig[1, 1], spec)
    isempty(column) || (fig[1, 2] = vgrid!(column...))
    # `below` controls span the plot column, stacked on their own rows underneath
    for (k, c) in enumerate(below)
        fig[1 + k, 1] = c
    end
    return fig
end

function ferriteviewer(ds::FEData, data::AbstractVector{<:AbstractVector};
                       layout=default_layout,
                       controls::AbstractVector{Control}=default_controls(ds),
                       pipeline=default_pipeline)
    return ferriteviewer(ds; layout, controls=Control[controls; TimeSlider(data)], pipeline)
end
