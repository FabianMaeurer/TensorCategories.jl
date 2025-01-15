#=----------------------------------------------------------
    Induction to the G-Equivariantization 
----------------------------------------------------------=#

function equivariant_induction(X::Object, T::GTensorAction, parent::Equivariantization = equivariantization(parent(X), T))

    G = group(T)

    object_IX, incl, _ = direct_sum([T(g)(X) for g ∈ G])

    structure_maps = Morphism[]

    for g ∈ G
        _,_,proj = direct_sum([T(g*h)(X) for h ∈ G])

        permutation = [findfirst(==(g*h), elements(G)) for h ∈ G]

        matched_incl = incl[permutation]

        u = sum([i ∘ monoidal_structure(T,g,h)(X) ∘ p for (h,i,p) ∈ zip(elements(G), matched_incl, proj)])

        push!(structure_maps, u)
    end

    return EquivariantObject(parent, object_IX, structure_maps)
end


function equivariant_induction_adjunction(X::Object, Y::EquivariantObject, T::GTensorAction, IX::EquivariantObject = equivariant_induction(X,T, parent(Y)))

    H = Hom(X, object(Y))

    G = group(T)

    _,_,proj = direct_sum([T(g)(X) for g ∈ G])

    base = Morphism[sum([u ∘ T(g)(f) ∘ p for (g,u,p) ∈ zip(elements(G), Y.structure_maps, proj)]) for f ∈ H]

    base = [EquivariantMorphism(IX, Y, m) for m ∈ base]
    HomSpace(IX, Y, base)
end