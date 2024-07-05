using Documenter, TensorCategories, Oscar, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "MyBib.bib"))

makedocs(
    plugins = [bib],
    sitename = "TensorCategories.jl",
    modules = [TensorCategories],
    warnonly = true,
    format = Documenter.HTML(
        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
        prettyurls = !("local" in ARGS)
    ),
    pages = [
        "Home" => "index.md",
        "Implementing Categories" => [
            "Philosophy" => "Interface/Philosophy.md",
            "Categories" => "Interface/Categories.md",
            "Abelian Categories" => "Interface/AbelianCategories.md",
            "Monoidal Categories" => "Interface/MonoidalCategories.md",
            "Linear Categories" => "Interface/LinearCategories.md",
            "Tensor Categories" => [
                "Framework" => "Interface/TensorCategories.md",
                "6j-Symbols" => "SixJCategories/SixJCategories.md"
            ],
            "Genericity" => "Interface/Generic.md",
            "Basic Constructions" => "Interface/BasicConstructions.md"
        ],
        "Examples" => [
            "Graded Vector Spaces" => "ConcreteExamples/VectorSpaces.md",
            "Group Representations" => "ConcreteExamples/Representations.md",
            "Equivariant Coherent Sheaves" => "ConcreteExamples/CoherentSheaves.md",
            "Fusion Categories with 6j-Symbols" => "SixJCategories/6JExamples.md",
            "Representations of  ``U_q(\\mathfrak{sl}_2(K))``" => "ConcreteExamples/UqSl2.md"
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
        "Categorical Constructions" => [
            "The Drinfeld Center" => "Constructions/Center.md",
            "The Drinfeld Centralizer" => "Constructions/Centralizer.md",
            "Internal Module Categories" => "Constructions/ModuleCategories.md"
        ],
        "References" => "References.md"
    ],
)


deploydocs(
    repo   = "github.com/FabianMaeurer/TensorCategories.jl.git",
)
