#=----------------------------------------------------------
    Checks for equivariance 
----------------------------------------------------------=#

function is_equivariant(X::EquivariantObject, Y::EquivariantObject, f::Morphism)

    T = gaction(parent(X))
    G = group(T)
    for (g,u,v) ∈ zip(elements(G), X.structure_maps, Y.structure_maps)
        f ∘ u != v ∘ T(g)(f) && return false
    end
    
    true
end

is_equivariant(f::EquivariantMorphism) = is_equivariant(domain(f), codomain(f), morphism(f))


function is_equivariant(X::Object, T::GTensorAction, structure::Vector{<:Morphism})

    G = group(T)

    elems = elements(G)

    for (g,u) ∈ zip(elems,structure), (h,v) ∈ zip(elems, structure)
        i = findfirst(==(g*h), elems)

        uv = structure[i]
        γ = monoidal_structure(T,g,h)(id(X))
        (u ∘ T(g)(v) != uv ∘ γ) && return false        
    end

    true
end
    
is_equivariant(X::EquivariantObject) = is_equivariant(object(X), gaction(parent(X)), X.structure_maps)