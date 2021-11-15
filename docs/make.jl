using Documenter, JuCat

makedocs(
    sitename = "JuCat.jl",
    modules = [JuCat],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
        assets = ["assets/favicon.ico"],
        analytics = "UA-136089579-2",
        highlights = ["yaml"],
    ),
    pages = [
        "Home" => "index.md",
        "types.md",
        "showcase.md"
    ],
)

deploydocs(
    repo = "github.com/FabianMaeurer/JuCat.jl",
)
