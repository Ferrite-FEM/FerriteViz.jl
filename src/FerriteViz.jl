module FerriteViz

using Makie
using Tensors
import Ferrite
import LinearAlgebra

abstract type AbstractPlotter end

include("utils.jl")
include("makieplotting.jl")
include("lor_tools.jl")

export MakiePlotter
export ferriteviewer

end
