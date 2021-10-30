module FerriteVis

using Makie
import Ferrite
import LinearAlgebra

abstract type AbstractPlotter end

include("utils.jl")
include("makieplotting.jl")

export MakiePlotter
export ferriteviewer

end
