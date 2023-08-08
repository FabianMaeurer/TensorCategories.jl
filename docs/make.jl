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
        "Category Interface" => [
            "Philosophy" => "Interface/Philosophy.md",
            "Categories" => "Interface/BasicInterface.md",
            "Abelian Categories" => "Interface/Abelian Categories.md",
            "Monoidal Categories" => "Interface/MonoidalCategories.md",
            "Tensor Categories" => "Interface/TensorCategories.md",
            "Optimisations" => "Interface/AdvancedInterface.md"
        ],
        # "Concrete Examples" => [
        #     "Vector Spaces" => "VectorSpaces.md",
        #     "Representations" => "Representations.md",
        #     "Coherent Sheaves" => "CoherentSheaves.md"
        # ],
        # "Fusion Categories from 6j Symbols" => [
        #     "Idea" => "SixJCategory.md",
        #     "Examples" => "RingCatExamples.md"
        # ],
        # "ℤ₊-Rings" => [
        #     "ℤ₊-Rings" => "ZPlusRings.md"
        # ],
       #"Multitensor Category Interface" => "Multitensor.md",
        "The Center Construction" => "Center.md",
    ],
)


deploydocs(
    repo   = "github.com/FabianMaeurer/TensorCategories.jl.git",
)
