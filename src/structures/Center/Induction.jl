
#-------------------------------------------------------------------------------
#   Induction
#-------------------------------------------------------------------------------

function induction(X::Object, simples::Vector = simples(parent(X)))
    @assert issemisimple(parent(X)) "Requires semisimplicity"
    Z = dsum([s⊗X⊗dual(s) for s ∈ simples])
    a = associator
    γ = Vector{Morphism}(undef, length(simples))

    for i ∈ 1:length(simples)
        W = simples[i]
        γ[i] = zero_morphism(zero(parent(X)), dsum([W⊗((s⊗X)⊗dual(s)) for s ∈ simples]))
  
        for S ∈ simples
            dom_i = ((S⊗X)⊗dual(S))⊗simples[i]
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ simples
                # Set up basis and dual basis
                basis, basis_dual = dual_basis(Hom(S, W⊗T), Hom(dual(S), dual(T)⊗dual(W)))

                # Correct dual basis to right (co)domain via Hom(U⊗V,W) ≃ Hom(U,W⊗V∗)
                basis_dual = [(id(dual(T))⊗ev(W)) ∘ a(dual(T),dual(W),W) ∘ (f⊗id(W)) for f ∈ basis_dual]

                if length(basis) == 0 
                    γ_i_temp = vertical_dsum(γ_i_temp, zero_morphism(dom_i, W⊗T⊗X⊗dual(T))) 
                else
                    component_iso = sum([(id(W)⊗inv(a(T,X,dual(T)))) ∘ a(W,T,X⊗dual(T)) ∘ a(W⊗T,X,dual(T)) ∘ (f⊗id(X)⊗g) ∘ a(S⊗X,dual(S),W) for (f,g) ∈ zip(basis, basis_dual)])
                   
                    γ_i_temp = vertical_dsum(γ_i_temp, (component_iso))
                end
            end
            γ[i] = horizontal_dsum(γ[i], dim(S)*γ_i_temp)
        end
        # distribution Before
        distr_before = distribute_left([s⊗X⊗dual(s) for s ∈ simples],W)
        # distribution After
        distr_after = distribute_right(W,[s⊗X⊗dual(s) for s ∈ simples])

        γ[i] = inv(distr_after) ∘ γ[i] ∘ distr_before
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


