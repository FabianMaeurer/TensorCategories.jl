using Documenter, JuCat

DocMeta.setdocmeta!(JuCat, :DocTestSetup, :(using JuCat,Oscar); recursive=true)

makedocs(
    sitename = "JuCat.jl",
    modules = [JuCat],
    format = Documenter.HTML(
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Vector Spaces" => "VectorSpaces.md",
            "Representations" => "Representations.md"
        ],
        "showcase.md"
    ],
)

deploydocs(
    repo   = "github.com/FabianMaeurer/JuCat.jl.git",
)
