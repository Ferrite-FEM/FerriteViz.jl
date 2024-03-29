# FerriteViz.jl changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.2] - 2023-11-10
### Added
 - uniform refinement for high-order solutions ([#97][github-97])
 - dependabot for GitHub actions ([#101][github-101])
 - attempt to increase internal machinery test coverage ([#104][github-104])

### Modified
 - `README.md` improvements with example gifs ([#96][github-96])
 - CI trigger only for PRs and master ([#105][github-105])
 - update docs to Documenter v1 ([#106][github-106])
 - update Makie in docs to v0.19.12 ([#109][github-109])

### Fixed
 - 0 `ntriangles` for empty domains ([#92][github-92])
 - correct link for plasticity example ([#93][github-93])
 - colorbar for 0 values and `ferriteviewer` deformation default changed to false ([#95][github-95])

## [0.2.1] - 2023-05-24
### Added
 - Basic culling where all faces of all boundary elements are rendered ([#56][github-56]).
 - Citation file ([#65](github-65))
 - Support for MixedDofHandler ([#70][github-70])
### Modified
 - Removed unnecessary extra dispatches for three-dimensional case ([#56][github-56]).
 - function barrier for `transfer_solution` such that its closer to type groundedness ([#68][github-68]).
 - `MakiePlotter` holds now `ShaderAbstractions.Buffer`s ([#69][github-69])
    - triangles are now stored in `Buffer`s with Observables
    - triangle coords are now `Buffers`s with Observables
- replace overcomplicated ternary operators by `begin end` expressions ([#69][github-69])
- remove unused functions ([#69][github-69])
- default linear rendering of high order triangles ([#83][github-83])
- keyword argument `copy_fields` added to `interpolate_gradient_field` ([#83][github-83])
### Fixed
 - Renamed `Crincle` to `Crinkle` ([#56][github-56]).
 - wireframe plot could not selectively disable the plotting of the nodes ([#83][github-83])
 - let CI error if example block errors ([#71][github-71])
 - removed bug in `transfer_solution` from ([#70][github-70]) in ([#89][github-89])
 - fix JSServe documentation issue ([#85][github-85])

## [0.2.0] - 2023-03-06
### Added
 - Functionality to obtain a first-order refined mesh and the corresponding
   dof handler and solution to approximately visualize high order solutions ([#57][github-57]).
 - Subtitles for the tutorial to find useful stuff faster ([#57][github-57]).
 - Crincle clip in 3D ([#56][github-56]), which basically removes all elements above some surface,
   which can be described by a function, from visualization.
 - `ClipPlane` to describe planar surfaces in 3D ([#56][github-56]).
 - Docs and helper for gradient field visualization based on interpolation ([#51][github-51]).
   Currently only useful in 2D, because we have no clip plane feature to introspect the interior
   of 3D problems.
 - Manufactured heat problem to test correctness of gradient field computation and as a
   helper to generate scalar-valued solutions with different ansatz ([#51][github-51]).

### Modified
 - Incompressible elasticity solver now takes the Ansatz functions and the actual material
   parameters instead of the poisson number the Ansatz functions ([#51][github-51]).

### Fixed
 - Visualization of non-conforming solution fields in 3D ([#59][github-59]).
 - An unknown bug has been fixed, which computes the colorbar `(min,max)` wrong. Now the `max` is
   set to be `1.01` of `min` guaranteeing that the value is larger than `min` if close to zero ([#51][github-51]).
 - Update Makie dependencies to fix some visualization bugs ([#51][github-51]).

[github-51]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/51
[github-56]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/56
[github-57]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/57
[github-59]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/59
[github-65]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/65
[github-63]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/63
[github-68]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/68
[github-69]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/69
[github-70]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/70
[github-71]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/71
[github-83]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/83
[github-85]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/85
[github-89]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/89
[github-92]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/92
[github-93]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/93
[github-95]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/95
[github-96]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/96
[github-97]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/97
[github-101]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/101
[github-104]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/104
[github-105]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/105
[github-106]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/106
[github-109]: https://github.com/Ferrite-FEM/FerriteViz.jl/pull/109

[Unreleased]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.1...v0.2.0
[0.2.0]: https://github.com/Ferrite-FEM/FerriteViz.jl/compare/v0.2.0...v0.1.4
