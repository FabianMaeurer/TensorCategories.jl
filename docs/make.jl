using Documenter

makedocs(
    sitename = "JuCat.jl",
    modules = [JuCat],
    pages = Any[
        "index.md",
        "types.md",
        "Showcase" => "showcase.md"
    ],
)
