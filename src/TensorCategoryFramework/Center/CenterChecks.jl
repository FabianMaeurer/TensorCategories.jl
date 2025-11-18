"""
    is_half_braiding(X::Object, half_braiding::Vector{<:Morphism})

TBW
"""
function is_half_braiding(Z::Object, half_braiding::Vector{<:Morphism}, log::Bool = false)
    if typeof(base_ring(Z)) <: Union{ArbField, ComplexField, AcbField}
        return is_half_braiding_numeric(Z, half_braiding, log)
    end

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
                if log == true
                    @show (i,j,k)
                    @show matrix(left)
                    @show matrix(right)
                    @show matrix(left-right)
                else
                    return false
                end
                flag = true
            end
        end
    end
    return !flag 
end

function is_half_braiding_numeric(Z::Object, half_braiding::Vector{<:Morphism}, log::Bool = false)
    simple_objects = simples(parent(Z))
    n = length(simple_objects)
    flag = false
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        Xᵢ,Xⱼ,Xₖ = simple_objects[[i,j,k]]
        γᵢ,γⱼ,γₖ = half_braiding[[i,j,k]]

        for t ∈ basis(Hom(Xₖ, Xᵢ⊗Xⱼ))
            left = associator(Xᵢ,Xⱼ,Z) ∘ (t⊗id(Z)) ∘ γₖ
            right = (id(Xᵢ)⊗γⱼ) ∘ associator(Xᵢ,Z,Xⱼ) ∘ (γᵢ⊗id(Xⱼ)) ∘ inv_associator(Z,Xᵢ,Xⱼ) ∘ (id(Z)⊗t)
            if !overlaps(matrix(left), matrix(right)) 
                if log == true
                    @show (i,j,k)
                    @show matrix(left)
                    @show matrix(right)
                    @show matrix(left-right)
                else
                    return false
                end
                flag = true
            end
        end
    end
    return !flag 
end

function is_relative_braiding(Z::Object, half_braiding::Vector{<:Morphism}, simple_objects::Vector{<:Object}, log::Bool = false)
    
    n = length(simple_objects)
    flag = false
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        Xᵢ,Xⱼ,Xₖ = simple_objects[[i,j,k]]
        γᵢ,γⱼ,γₖ = half_braiding[[i,j,k]]

        for t ∈ basis(Hom(Xₖ, Xᵢ⊗Xⱼ))
            left = associator(Xᵢ,Xⱼ,Z) ∘ (t⊗id(Z)) ∘ γₖ
            right = (id(Xᵢ)⊗γⱼ) ∘ associator(Xᵢ,Z,Xⱼ) ∘ (γᵢ⊗id(Xⱼ)) ∘ inv_associator(Z,Xᵢ,Xⱼ) ∘ (id(Z)⊗t)
            if left != right 
                if log == true
                    @show (i,j,k)
                    @show matrix(left)
                    @show matrix(right)
                    @show matrix(left-right)
                else
                    return false
                end
                flag = true
            end
        end
    end
    return !flag 
end

is_central(X::CentralizerObject, log = false) = is_relative_braiding(object(X), half_braiding(X), parent(X).subcategory_simples, log)

is_central(X::CenterObject, log = false) = is_half_braiding(object(X), half_braiding(X), log)

function is_central(f::Morphism, X::CenterObject, Y::CenterObject)
    if typeof(base_ring(f)) <: Union{ArbField, ComplexField, AcbField}
        return is_central_numeric(f, X,Y)
    end

    C = category(parent(X))
    S = simples(C)

    for (s,γ,δ) ∈ zip(S, half_braiding(X), half_braiding(Y))
        if (id(s)⊗f) ∘ γ != δ ∘ (f ⊗ id(s))
            return false
        end
    end
    true
end

function is_central_numeric(f::Morphism, X::CenterObject, Y::CenterObject)
    C = category(parent(X))
    S = simples(C)

    for (s,γ,δ) ∈ zip(S, half_braiding(X), half_braiding(Y))
        if !overlaps(matrix((id(s)⊗f) ∘ γ), matrix(δ ∘ (f ⊗ id(s))))
            # @show matrix((id(s)⊗f) ∘ γ)- matrix(δ ∘ (f ⊗ id(s)))
            return false
        end
    end
    true
end

is_central(f::CenterMorphism) = is_central(morphism(f), domain(f), codomain(f))