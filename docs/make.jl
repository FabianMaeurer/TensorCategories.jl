using Documenter, DocumenterLaTeX, JuCat

#DocMeta.setdocmeta!(JuCat, :DocTestSetup, :(using JuCat); recursive=true)

makedocs(
    sitename = "JuCat.jl",
    modules = [JuCat],
    format = Documenter.HTML(
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "Exampes" => [
            "Vector Spaces" => "VectorSpaces.md",
            "Representations" => "Representations.md",
            "Coherent Sheaves" => "CoherentSheaves.md"
        ],
        "Multitensor Categories" => "Multitensor.md",
        "showcase.md"
    ],
)

makedocs(
    sitename = "JuCat.jl",
    modules = [JuCat],
    format = LaTeX(platform = "none"),
    pages = [
            "Home" => "index.md",
            "Exampes" => [
                "Vector Spaces" => "VectorSpaces.md",
                "Representations" => "Representations.md",
                "Coherent Sheaves" => "CoherentSheaves.md"
            ],
            "Multitensor Categories" => "Multitensor.md",
            "showcase.md"
        ],
)

deploydocs(
    repo   = "github.com/FabianMaeurer/JuCat.jl.git",
)
