module FerriteVis

import Makie
import Ferrite

abstract type AbstractPlotter end

include("utils.jl")
include("makieplotting.jl")

export MakiePlotter
export warp_by_vector, warp_by_vector!, plot_grid, plot_grid!

end
