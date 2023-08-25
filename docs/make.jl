using Documenter, TensorCategories, Oscar, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "MyBib.bib"))

makedocs(
    bib,
    sitename = "TensorCategories.jl",
    modules = [TensorCategories],
    format = Documenter.HTML(
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
        prettyurls = !("local" in ARGS)
    ),
    pages = [
        "Home" => "index.md",
        "Category Interface" => [
            "Philosophy" => "Interface/Philosophy.md",
            "Categories" => "Interface/Categories.md",
            "Abelian Categories" => "Interface/AbelianCategories.md",
            "Monoidal Categories" => "Interface/MonoidalCategories.md",
            "Linear Categories" => "Interface/LinearCategories.md",
            "Tensor Categories" => "Interface/TensorCategories.md",
            "Optimisations" => "Interface/AdvancedInterface.md"
        ],
        "Examples" => [
            "Graded Vector Spaces" => "ConcreteExamples/VectorSpaces.md",
            "Group Representations" => "ConcreteExamples/Representations.md",
            "Equivariant Coherent Sheaves" => "ConcreteExamples/CoherentSheaves.md"
        ],
        "6j-Symbols" => [
            "SixJCategories/SixJCategories.md",
            "Examples" => "SixJCategories/6JExamples.md"
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
        "References" => "References.md"
    ],
)


deploydocs(
    repo   = "github.com/FabianMaeurer/TensorCategories.jl.git",
)
