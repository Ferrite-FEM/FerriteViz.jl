module FerriteViz

using Makie
using Tensors
import Ferrite
import LinearAlgebra

abstract type AbstractPlotter end

include("utils.jl")
include("makieplotting.jl")

export MakiePlotter
export ferriteviewer

end
