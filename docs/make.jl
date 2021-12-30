using Documenter, JuCat

makedocs(
    sitename = "JuCat.jl",
    modules = [JuCat],
    format = Documenter.HTML(
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "Implemented Examples" => [
            "Vector Spaces" => "VectorSpaces.md",
            "Representations of Finite Groups" => "Representations.md"
        ],
        "showcase.md"
    ],
)

deploydocs(
    repo   = "github.com/FabianMaeurer/JuCat.jl.git",
)
