using Documenter, TensorCategories, Oscar, DocumenterCitations

bib = CitationBibliography("MyBib.bib")

makedocs(
    prettyurls = !("local" in ARGS),
    bib,
    sitename = "TensorCategories.jl",
    modules = [TensorCategories],
    format = Documenter.HTML(
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "Introduction" => "Interface.md",
        "Concrete Examples" => [
            "Vector Spaces" => "VectorSpaces.md",
            "Representations" => "Representations.md",
            "Coherent Sheaves" => "CoherentSheaves.md"
        ],
        "Fusion Categories from 6j Symbols" => [
            "Idea" => "RingCategory.md",
            "Examples" => "RingCatExamples.md"
        ],
        "ℤ₊-Rings" => [
            "ℤ₊-Rings" => "ZPlusRings.md"
        ],
       #"Multitensor Category Interface" => "Multitensor.md",
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
