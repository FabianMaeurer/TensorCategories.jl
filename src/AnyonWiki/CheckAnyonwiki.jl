#=----------------------------------------------------------
    Check the properties of AnyonWiki categories
----------------------------------------------------------=#

labels = ["pentagon axiom", "pivotal", "spherical"]
braided_labels = ["hexagon_axiom", "modular"]

cyclic_n = parse.(Int, readlines("src/AnyonWiki/AnyonWikiData/cyclopos.txt"))

for n ∈ cyclic_n
    attrs = TensorCategories.load_anyonwiki_attributes(n)
    C = TensorCategories.anyonwiki(n)

    if ! pentagon_axiom(C)
        
        println("n = $n: pentagon axiom failed")
        continue
    end

    if ! is_pivotal(C) 
        println("n = $n: pivotal axiom failed")
        continue 
    end

    if ! (is_spherical(C) == attrs[:spherical])
        println("n = $n: spherical in anonwiki mismatched: $(attrs[:spherical]) ")
        continue
    end

    if attrs[:braided]
        
        if ! hexagon_axiom(C) 
            println("n = $n: hexagon axiom failed")
            continue
        end

        if attrs[:spherical] && ! (is_modular(C) == attrs[:modular])
            println("n = $n: modularity mismatched: $(attrs[:modular])")
            continue 
        end
    end
end