
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
        γ[i] = zero_morphism(zero(parent(X)), dsum([((dual(s)⊗X)⊗s)⊗W for s ∈ simples]))

        for S ∈ simples
            dom_i = simples[i]⊗((dual(S)⊗X)⊗S)
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ simples
                # Set up basis and dual basis
                basis, basis_dual = dual_basis(Hom(S, T⊗W), Hom(dual(S), dual(W)⊗dual(T)))

                # Correct dual basis to right (co)domain via Hom(U⊗V,W) ≃ Hom(U,W⊗V∗)
                basis_dual = [(ev(dual(W))⊗id(dual(T))) ∘ inv(a(dual(dual(W)),dual(W),dual(T))) ∘ ((spherical(W)∘id(W))⊗f) for f ∈ basis_dual]

                if length(basis) == 0 
                    γ_i_temp = vertical_dsum(γ_i_temp, zero_morphism(dom_i, ((dual(T)⊗X)⊗T)⊗W)) 
                else
                    component_iso = sum([inv(a(dual(T)⊗X,T,W)) ∘ (g⊗id(X)⊗f) ∘ (inv(a(W,dual(S),X))⊗id(S)) ∘ inv(a(W,dual(S)⊗X,S)) for (f,g) ∈ zip(basis, basis_dual)])

                    γ_i_temp = vertical_dsum(γ_i_temp, (dim(S)*component_iso))
                end
            end
            γ[i] = horizontal_dsum(γ[i], γ_i_temp)
        end
        # distribution Before
        distr_before = distribute_right(W,[dual(s)⊗X⊗s for s ∈ simples])
        @show matrix(distr_before)
        # distribution After
        distr_after = distribute_left([dual(s)⊗X⊗s for s ∈ simples],W)
        @show matrix(distr_after)

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


