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
FerriteViz.ferriteviewer
```
