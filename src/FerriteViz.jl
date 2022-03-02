module FerriteViz

using Makie
using Tensors
using PGFPlotsX
import Ferrite
import LinearAlgebra

abstract type AbstractPlotter end

include("utils.jl")
include("makieplotting.jl")
include("plotsplotting.jl")

export MakiePlotter
export ferriteviewer

end
