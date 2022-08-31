
#-------------------------------------------------------------------------------
#   Induction
#-------------------------------------------------------------------------------

function induction(X::Object, simples::Vector = simples(parent(X)))
    @assert issemisimple(parent(X)) "Requires semisimplicity"
    Z = dsum([s⊗X⊗dual(s) for s ∈ simples])

    γ = Vector{Morphism}(undef, length(simples))

    for i ∈ 1:length(simples)
        γ[i] = zero_morphism(zero(parent(X)), dsum([simples[i]⊗s⊗X⊗dual(s) for s ∈ simples]))

        for S ∈ simples
            dom_i = S⊗X⊗dual(S)⊗simples[i]
            @show γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ simples
                # Set up basis and dual basis
                basis, basis_dual = dual_basis(Hom(S, simples[i]⊗T), Hom(dual(S), dual(T)⊗dual(simples[i])))

                # Correct dual basis to right (co)domain via Hom(U⊗V,W) ≃ Hom(U,W⊗V∗)
                basis_dual = [(id(dual(T))⊗ev(simples[i])) ∘ associator(dual(T),dual(simples[i]),simples[i]) ∘ (f⊗id(simples[i])) for f ∈ basis_dual]

                if length(basis) == 0 
                    γ_i_temp = vertical_dsum(γ_i_temp, zero_morphism(dom_i, simples[i]⊗T⊗X⊗dual(T))) 
                else
                    γ_i_temp = vertical_dsum(γ_i_temp, (dim(S)*sum([f⊗id(X)⊗g for (f,g) ∈ zip(basis, basis_dual)])))
                end
            end
            γ[i] = horizontal_dsum(γ[i], γ_i_temp)
        end
        iso1 = isisomorphic(Z⊗simples[i], domain(γ[i]))[2] #inv(decompose_morphism(domain(γ[i])))∘decompose_morphism(Z⊗simples[i])
        iso2 = isisomorphic(codomain(γ[i]), simples[i]⊗Z)[2] #inv(decompose_morphism(simples[i]⊗Z))∘decompose_morphism(codomain(γ[i]))

        γ[i] = iso2 ∘ γ[i] ∘ iso1

    end
        
    return CenterObject(CenterCategory(base_ring(X),parent(X)),Z,γ)
end


#-------------------------------------------------------------------------------
#   Pairing and Dual Basis
#-------------------------------------------------------------------------------

function pairing(f::Morphism, g::Morphism)
    A,B = domain(f), codomain(g)
    return ev(B)∘(f⊗g)∘coev(A)
end

function dual_basis(V::AbstractHomSpace, W::AbstractHomSpace)
    dual_basis = []
    F = base_ring(V)
    n = length(basis(V))
    m = length(basis(W))
    for i ∈ 1:m
        M = matrix(F, n,m, [pairing(basis(V)[j], basis(W)[k]) for j ∈ 1:n, k ∈ 1:m])
        b = matrix(F,m,1, [i == j ? 1 : 0 for j ∈ 1:m])
        coeffs = solve(M,b)
        push!(dual_basis, sum(coeffs .* basis(W)))
    end
    return basis(V), dual_basis
end


