# API Reference

## Data source

```@docs
FerriteViz.FEData
FerriteViz.update!
FerriteViz.point_data
FerriteViz.cell_data
FerriteViz.set_point_data!
FerriteViz.set_cell_data!
```

## Filters

Filters transform an [`FEData`](@ref) into a new one and compose with `|>`:

```julia
ds = FEData(dh, u)
solutionplot(ds |> WarpByVector(:u, 2.0) |> Gradient(:u) |> VonMises(); color=:vonMises)
```

```@docs
FerriteViz.apply
FerriteViz.WarpByVector
FerriteViz.Gradient
FerriteViz.CrinkleClip
FerriteViz.ClipPlane
FerriteViz.Refine
FerriteViz.FirstOrderRefinement
FerriteViz.Component
FerriteViz.Magnitude
FerriteViz.Norm1
FerriteViz.VonMises
FerriteViz.Deviator
FerriteViz.Threshold
FerriteViz.Derive
FerriteViz.vonmises
```

## Representations

```@docs
FerriteViz.solutionplot
FerriteViz.cellplot
FerriteViz.meshplot
FerriteViz.arrowplot
FerriteViz.surfaceplot
FerriteViz.elementinfo
```

## Composable viewer

The viewer is assembled declaratively with `Makie.SpecApi`: pluggable
[`FerriteViz.Control`](@ref)s feed a view-state observable, a `layout(ds, state)`
hook returns a `GridLayoutSpec`, and spec helpers build the panels.

```@docs
FerriteViz.ferriteviewer
FerriteViz.Control
FerriteViz.ControlResult
FerriteViz.default_controls
FerriteViz.default_layout
FerriteViz.default_pipeline
FerriteViz.FieldMenu
FerriteViz.ProcessMenu
FerriteViz.ColormapMenu
FerriteViz.WireframeToggle
FerriteViz.LabelsToggle
FerriteViz.DeformationToggle
FerriteViz.TimeSlider
FerriteViz.panelspec
FerriteViz.solutionplotspec
FerriteViz.meshplotspec
FerriteViz.cellplotspec
FerriteViz.surfaceplotspec
FerriteViz.arrowplotspec
```
