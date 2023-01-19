# FerriteViz.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]
### Added
 - Docs and helper for gradient field visualization based on interpolation [#51][github-51].
   Currently only useful in 2D, because we have no clip plane feature to introspect the interior 
   of 3D problems.
 - Manufactured heat problem to test correctness of gradient field computation and as a 
   helper to generate scalar-valued solutions with different ansatz [#51][github-51].

### Modified
 - Incompressible elasticity solver now takes the Ansatz functions and the actual material 
   parameters instead of the poisson number the Ansatz functions [#51][github-51].


[github-51](https://github.com/Ferrite-FEM/Ferrite.jl/pull/51)
