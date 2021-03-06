using Documenter, TensorCategories

#DocMeta.setdocmeta!(JuCat, :DocTestSetup, :(using JuCat); recursive=true)

makedocs(
    sitename = "TensorCategories.jl",
    modules = [TensorCategories],
    format = Documenter.HTML(
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "Exampes" => [
            "Vector Spaces" => "VectorSpaces.md",
            "Representations" => "Representations.md",
            "Coherent Sheaves" => "CoherentSheaves.md",
            "Ring Categories" => "RingCategories.md"
        ],
        "Multitensor Categories" => "Multitensor.md",
        "The Center Construction" => "Center.md",
    ],
)

# makedocs(
#     sitename = "JuCat.jl",
#     modules = [JuCat],
#     format = LaTeX(platform = "none"),
#     pages = [
#             "Home" => "index.md",
#             "Exampes" => [
#                 "Vector Spaces" => "VectorSpaces.md",
#                 "Representations" => "Representations.md",
#                 "Coherent Sheaves" => "CoherentSheaves.md"
#             ],
#             "Multitensor Categories" => "Multitensor.md",
#             "showcase.md"
#         ],
# )

deploydocs(
    repo   = "github.com/FabianMaeurer/TensorCategories.jl.git",
)
