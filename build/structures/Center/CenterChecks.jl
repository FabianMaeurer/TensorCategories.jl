"""
    is_half_braiding(X::CategoryObject, half_braiding::Vector{<:CategoryMorphism})

TBW
"""
function is_half_braiding(Z::CategoryObject, half_braiding::Vector{<:CategoryMorphism})
    simple_objects = simples(parent(Z))
    n = length(simple_objects)
    flag = false
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        Xᵢ,Xⱼ,Xₖ = simple_objects[[i,j,k]]
        γᵢ,γⱼ,γₖ = half_braiding[[i,j,k]]

        for t ∈ basis(Hom(Xₖ, Xᵢ⊗Xⱼ))
            left = associator(Xᵢ,Xⱼ,Z) ∘ (t⊗id(Z)) ∘ γₖ
            right = (id(Xᵢ)⊗γⱼ) ∘ associator(Xᵢ,Z,Xⱼ) ∘ (γᵢ⊗id(Xⱼ)) ∘ inv_associator(Z,Xᵢ,Xⱼ) ∘ (id(Z)⊗t)
            if left != right
                @show (i,j,k)
                @show matrix(left)
                @show matrix(right)
                flag = true
            end
        end
    end
    return !flag 
end

is_central(X::CenterCategoryObject) = is_half_braiding(object(X), half_braiding(X))