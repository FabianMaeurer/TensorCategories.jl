"""
    is_half_braiding(X::Object, half_braiding::Vector{<:Morphism})

TBW
"""
function is_half_braiding(Z::Object, half_braiding::Vector{<:Morphism})
    simple_objects = simples(parent(Z))
    n = length(simple_objects)
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        Xᵢ,Xⱼ,Xₖ = simple_objects[[i,j,k]]
        γᵢ,γⱼ,γₖ = half_braiding[[i,j,k]]

        for t ∈ basis(Hom(Xₖ, Xᵢ⊗Xⱼ))
            left = associator(Xᵢ,Xⱼ,Z) ∘ (t⊗id(Z)) ∘ γₖ
            right = (id(Xᵢ)⊗γⱼ) ∘ associator(Xᵢ,Z,Xⱼ) ∘ (γᵢ⊗id(Xⱼ)) ∘ inv(associator(Z,Xᵢ,Xⱼ)) ∘ (id(Z)⊗t)
            if left != right
                @show matrix(left)
                @show matrix(right)
                return false
            end
        end
    end
    return true
end