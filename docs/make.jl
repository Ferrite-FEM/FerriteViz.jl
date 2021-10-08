push!(LOAD_PATH,"../src/")
using Documenter, FerriteVis, GLMakie

makedocs(sitename="FerriteVis",
         modules=[FerriteVis],
         authors="Maximilian Köhler",
         pages=["Home"=> "index.md",
                "Tutorial" => "tutorial.md",]
)

deploydocs(repo = "github.com/koehlerson/FerriteVis.jl.git",
           push_preview=true,)
