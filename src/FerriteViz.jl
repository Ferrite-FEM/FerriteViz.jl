module FerriteViz

using Makie
using Tensors
using Ferrite
import Ferrite: get_grid, getrefshape
import GeometryBasics
import ShaderAbstractions
import LinearAlgebra

abstract type AbstractPlotter end
abstract type AbstractFilter end

include("tessellation.jl")
include("dataset.jl")
include("gradient.jl")
include("filters.jl")
include("representations.jl")

export FEData, apply
export set_point_data!, set_cell_data!
export WarpByVector, Gradient, CrinkleClip, Refine, FirstOrderRefinement
export Component, Magnitude, Norm1, VonMises, Deviator, Threshold, Derive
export ClipPlane, vonmises
export ferriteviewer
export solutionplot, solutionplot!, cellplot, cellplot!, meshplot, meshplot!,
       surfaceplot, surfaceplot!, arrowplot, arrowplot!, elementinfo, elementinfo!

end
