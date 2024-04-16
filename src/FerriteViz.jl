module FerriteViz

using Makie
using StaticArrays
using Tensors
using Ferrite
import Ferrite: get_grid, getrefshape
import GeometryBasics
import ShaderAbstractions
import LinearAlgebra

abstract type AbstractPlotter end

include("utils.jl")
include("makieplotting.jl")
include("lor_tools.jl")

export MakiePlotter
export ferriteviewer
export for_discretization

end
