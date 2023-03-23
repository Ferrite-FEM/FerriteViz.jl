push!(LOAD_PATH,"../src/")
using Documenter, FerriteViz, WGLMakie

makedocs(sitename="FerriteViz",
         modules=[FerriteViz],
         authors="Maximilian Köhler",
         pages=["Home"=> "index.md",
                "Tutorial" => "tutorial.md",
                "Advanced Topics" => "atopics.md",
                "API Reference" => "api.md",],
         strict=:example_block,
)

deploydocs(repo = "github.com/Ferrite-FEM/FerriteViz.jl.git",
           push_preview=true,
           forcepush=true,)
