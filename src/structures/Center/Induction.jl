
#-------------------------------------------------------------------------------
#   Induction
#-------------------------------------------------------------------------------

function induction(X::Object, simples::Vector = simples(parent(X)))
    @assert issemisimple(parent(X)) "Requires semisimplicity"
    Z = dsum([dual(s)⊗X⊗s for s ∈ simples])
    a = associator
    γ = Vector{Morphism}(undef, length(simples))
    C = parent(X)
    for i ∈ 1:length(simples)
        W = simples[i]
        γ[i] = zero_morphism(zero(C), dsum([((dual(s)⊗X)⊗s)⊗W for s ∈ simples]))

        for S ∈ simples
            dom_i = W⊗((dual(S)⊗X)⊗S)
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ simples
                # Set up basis and dual basis
                basis, basis_dual = dual_basis(Hom(one(C), (dual(S)⊗W)⊗T), Hom(one(C), dual(T)⊗(dual(W)⊗S)))
 
                # Correct basis
                basis = [(ev(S)⊗id(W⊗T)) ∘ (inv(a(S,dual(S),W))⊗id(T)) ∘ inv(a(S,dual(S)⊗W,T)) ∘ (id(S)⊗g) for g ∈ basis]
                # Correct dual basis to right (co)domain via Hom(U⊗V,W) ≃ Hom(U,U*⊗W)
                basis_dual = [(id(T)⊗ev(dual(W))) ∘ a(dual(T),dual(W),dual(dual(W))) ∘ (id(dual(T))⊗(id(dual(W))⊗ev(S))⊗id(dual(dual(W)))) ∘ (id(dual(T))⊗a(dual(W),S,dual(S))⊗id(dual(dual(W)))) ∘ (a(dual(T),dual(W)⊗S,dual(S))⊗id(dual(dual(W)))) ∘ (f⊗id(dual(S))⊗spherical(W)) for f ∈ basis_dual]

                if length(basis) == 0 
                    γ_i_temp = vertical_dsum(γ_i_temp, zero_morphism(dom_i, ((dual(T)⊗X)⊗T)⊗W)) 
                else
                    component_iso = sum([inv(a(dual(T)⊗X,T,W)) ∘ (g⊗id(X)⊗f) ∘ (inv(a(W,dual(S),X))⊗id(S)) ∘ inv(a(W,dual(S)⊗X,S)) for (f,g) ∈ zip(basis, basis_dual)])

                    γ_i_temp = vertical_dsum(γ_i_temp, component_iso)
                end
            end
            γ[i] = horizontal_dsum(γ[i], γ_i_temp)
        end
        # distribution Before
        distr_before = distribute_right(W,[dual(s)⊗X⊗s for s ∈ simples])
        # distribution After
        distr_after = distribute_left([dual(s)⊗X⊗s for s ∈ simples],W)

        γ[i] = inv(distr_after) ∘ γ[i] ∘ distr_before
    end
        
    return CenterObject(CenterCategory(base_ring(X),parent(X)),Z,γ)
end


#-------------------------------------------------------------------------------
#   Pairing and Dual Basis
#-------------------------------------------------------------------------------

function pairing(f::Morphism, g::Morphism)
    A,B = domain(f), codomain(g)
    @show codomain(f) == dual(codomain(g))
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


