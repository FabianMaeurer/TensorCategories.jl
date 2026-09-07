using Documenter, TensorCategories, Oscar, DocumenterCitations

DocMeta.setdocmeta!(TensorCategories, :DocTestSetup,
    :(using TensorCategories, Oscar); recursive=true, warn=false)
bib = CitationBibliography(joinpath(@__DIR__, "src", "MyBib.bib");
    style=:authoryear)

function check_rendered_example_outputs(build_dir)
    raw_output = r"</code></pre>\s*[^<\s]"
    for (root, _, files) in walkdir(build_dir)
        for file in files
            endswith(file, ".html") || continue
            path = joinpath(root, file)
            occursin(raw_output, read(path, String)) || continue
            error("unformatted example output in $(relpath(path, build_dir))")
        end
    end
end

# Home links to the repository's legal notices. Stage these without maintaining
# duplicate source copies or changing the Home page.
legal_files = ("LICENSE", "COPYRIGHT")
try
    for name in legal_files
        cp(joinpath(@__DIR__, "..", name), joinpath(@__DIR__, "src", name); force=true)
    end
    makedocs(
        root = @__DIR__,
        plugins = [bib],
        sitename = "TensorCategories.jl",
        modules = [TensorCategories],
        checkdocs = :none,
        warnonly = false,
        format = Documenter.HTML(
            edit_link = "master",
            canonical = "https://TensorCategories.github.io/TensorCategories.jl/stable/",
            prettyurls = !("local" in ARGS),
            collapselevel = 1,
            mathengine = MathJax3(Dict(
                :tex => Dict(
                    "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
                    "tags" => "all",
                    "packages" => ["base", "ams", "autoload"],
                ),
            )),
        ),
        pages = [
            "Home" => "index.md",
            "Basics" => [
                "Julia and OSCAR" => "Basics/Julia.md",
                "Implementing categories" => "Interface/Categories.md",
                "Linear and abelian categories" => "Interface/LinearCategories.md",
                "Coefficient fields and numeric computations" => "Basics/BaseFields.md",
                "Simple objects and finite length" => "Interface/SimpleObjects.md",
                "Idempotents and Krull–Schmidt categories" => "Interface/KaroubianCategories.md",
                "Finite and semisimple categories" => "Interface/SemisimpleCategories.md",
                "Splitting and scalar extension" => "Interface/SplittingFields.md",
                "Functors and natural transformations" => "Interface/Functors.md",
            ],
            "Monoidal categories" => [
                "Tensor products and associators" => "Interface/MonoidalCategories.md",
                "Implementing a monoidal category" => "Implementing/MonoidalCategory.md",
                "Rigidity" => "Interface/RigidCategories.md",
                "Ring and tensor categories" => "Interface/RingCategories.md",
                "Fusion and multifusion categories" => "Interface/TensorCategories.md",
                "Grothendieck rings" => "Interface/GrothendieckRings.md",
                "Computing with fusion rings" => "Interface/FusionRings.md",
                "Pivotal and spherical structures" => "Interface/PivotalCategories.md",
                "Braided and symmetric categories" => "Interface/BraidedCategories.md",
                "Unitarity and modularity" => "Interface/UnitaryModularCategories.md",
                "Monoidal functors" => "Interface/AdvancedInterface.md",
                "Fiber functors" => "Interface/FiberFunctors.md",
            ],
            "Skeletal fusion categories" => [
                "Skeletal models and conventions" => "F-symbols/SkeletalFusion.md",
                "First computations" => "Introduction/Introduction.md",
                "Working with fusion data" => "F-symbols/WorkedExamples.md",
                "Anyons and CFT terminology" => "Basics/Physics.md",
                "Numerical fusion categories" => "F-symbols/Numerical.md",
                "Data exchange" => "F-symbols/Data.md",
            ],
            "Drinfeld centers" => [
                "Half-braidings and computation" => "Constructions/Center.md",
                "The Ising center and splitting" => "Introduction/Ising.md",
                "Relative centers" => "Constructions/Centralizer.md",
            ],
            "Further constructions" => [
                "Products and scalar extension" => "Interface/BasicConstructions.md",
                "Algebras and internal modules" => "Constructions/ModuleCategories.md",
                "Group actions" => "Constructions/GroupActions.md",
                "Interface checklist" => "Interface/Generic.md",
            ],
            "Catalogue" => [
                "Inventory" => "F-symbols/Examples.md",
                "AnyonWiki" => "F-symbols/AnyonWiki.md",
                "E₆" => "F-symbols/E6.md",
                "Equivariant sheaves and convolution" => "ConcreteExamples/CoherentSheaves.md",
                "Fibonacci" => "F-symbols/Fibonacci.md",
                "Group representations" => "ConcreteExamples/Representations.md",
                "Haagerup" => "F-symbols/Haagerup.md",
                "sl₂, Verlinde, and dihedral models" => "ConcreteExamples/UqSl2.md",
                "SU(3)₃ subcategory" => "F-symbols/SU3_3.md",
                "Tambara–Yamagami and Ising" => "F-symbols/TambaraYamagami.md",
                "Trivial fusion category" => "F-symbols/Trivial.md",
                "Vector spaces and gradings" => "ConcreteExamples/VectorSpaces.md",
                "Vercleyen–Slingerland data" => "F-symbols/VercleyenSlingerland.md",
            ],
            "API reference" => [
                "Index" => "API.md",
                "Category interface" => "API/Framework.md",
                "Tensor structure" => "API/TensorStructure.md",
                "Concrete categories" => "API/Categories.md",
                "Fusion data and databases" => "API/FusionData.md",
                "Constructions" => "API/Constructions.md",
                "Utilities" => "API/Utilities.md",
            ],
            "Project" => [
                "Developing" => "Project/Developing.md",
                "Citations" => "Project/Citations.md",
                "Further literature" => "Project/FurtherLiterature.md",
            ],
            "References" => "References.md",
        ],
    )
    check_rendered_example_outputs(joinpath(@__DIR__, "build"))
finally
    for name in legal_files
        rm(joinpath(@__DIR__, "src", name); force=true)
    end
end

if "deploy" in ARGS
    deploydocs(
        repo = "github.com/TensorCategories/TensorCategories.jl.git",
        devbranch = "master",
    )
end
