push!(LOAD_PATH,"../src/")
using Documenter, FerriteViz, WGLMakie

makedocs(sitename="FerriteViz",
         modules=[FerriteViz],
         authors="Maximilian Köhler",
         format=Documenter.HTML(
              prettyurls=false,
              assets = ["assets/custom.css", "assets/favicon.ico"],
         ),
         pages=["Home"=> "index.md",
                "Tutorial" => "tutorial.md",
                "Advanced Topics" => "atopics.md",
                "api.md",
                "devdocs.md",
                ],
)

deploydocs(repo = "github.com/Ferrite-FEM/FerriteViz.jl.git",
           push_preview=true,
           forcepush=true,)
