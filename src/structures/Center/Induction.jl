
#-------------------------------------------------------------------------------
#   Induction
#-------------------------------------------------------------------------------

function induction(X::Object, simples::Vector = simples(parent(X)); parent_category::CenterCategory = Center(parent(X)))
    @assert issemisimple(parent(X)) "Requires semisimplicity"
    Z = dsum([s⊗X⊗dual(s) for s ∈ simples])
    a = associator
    γ = Vector{Morphism}(undef, length(simples))
    C = parent(X)
    for i ∈ 1:length(simples)
        W = simples[i]
        γ[i] = zero_morphism(zero(parent(X)), dsum([W⊗((s⊗X)⊗dual(s)) for s ∈ simples]))
  
        for S ∈ simples
            dom_i = ((S⊗X)⊗dual(S))⊗simples[i]
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ simples
                #@show S,T,W
                # Set up basis and dual basis
                
                #basis, basis_dual = dual_basis(Hom(S, W⊗T), Hom(dual(S), dual(W⊗T)))
                _basis, basis_dual = basis(Hom(S, W⊗T)), basis(Hom(dual(S)⊗W, dual(T)))

                if length(_basis) == 0 
                    γ_i_temp = vertical_dsum(γ_i_temp, zero_morphism(dom_i, W⊗((T⊗X)⊗dual(T)))) 
                    continue
                end

                corrections = base_ring(X).([(id(W)⊗ev(dual(T))) ∘ (id(W)⊗(spherical(T)⊗id(dual(T)))) ∘ a(W,T,dual(T)) ∘ (f⊗g) ∘ a(S,dual(S),W) ∘ (coev(S)⊗id(W)) for (f,g) ∈ zip(_basis,basis_dual)])
                
                #@show "component_iso"
                component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(inv(dim(W))*inv(k)*g)) ∘ a(S⊗X,dual(S),W) for (k,f,g) ∈ zip(corrections,_basis, basis_dual)])
                #@show "sum"
                γ_i_temp = vertical_dsum(γ_i_temp, (component_iso))
            end
            γ[i] = horizontal_dsum(γ[i], dim(S)*γ_i_temp)
        end
      

        # distribution Before
        distr_before = distribute_left([s⊗X⊗dual(s) for s ∈ simples],W)
        # distribution After
        distr_after = distribute_right(W,[s⊗X⊗dual(s) for s ∈ simples])

        γ[i] = inv(distr_after) ∘ γ[i] ∘ distr_before 
    end

    return CenterObject(parent_category,Z,γ)
end

function induction_restriction(X::Object, simples::Vector = simples(parent(X)))
    @assert issemisimple(parent(X)) "Requires semisimplicity"
    Z = dsum([s⊗X⊗dual(s) for s ∈ simples])
end
# function induction(X::Object, simples::Vector = simples(parent(X)))
#     @assert issemisimple(parent(X)) "Requires semisimplicity"
#     Z = dsum([s⊗X⊗dual(s) for s ∈ simples])
#     a = associator
#     γ = Vector{Morphism}(undef, length(simples))
#     C = parent(X)
#     for i ∈ 1:length(simples)
#         W = simples[i]
#         γ[i] = zero_morphism(zero(parent(X)), dsum([W⊗((s⊗X)⊗dual(s)) for s ∈ simples]))
  
#         for S ∈ simples
#             dom_i = ((S⊗X)⊗dual(S))⊗simples[i]
#             γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
#             for  T ∈ simples
#                 @show dim(S),dim(T),dim(W)
#                 # Set up basis and dual basis
#                 basis, basis_dual = dual_basis(Hom(one(C), (dual(S)⊗W)⊗T), Hom(one(C), dual((dual(S)⊗W)⊗T)))
#                 @show "dual_transform"
#                 # Isomorphism (XY)* ≃ Y*X*
#                 dual_transform = inv(a(dual(T),dual(W),S)) ∘ (id(dual(T))⊗dual_monoidal_structure(dual(S),W)) ∘ dual_monoidal_structure(dual(S)⊗W,T)
#                 basis_dual = [dual_transform ∘ f for f ∈ basis_dual]

#                 @show "basis correction"
#                 # Correct basis
#                 basis = [(ev(dual(S))⊗id(W)⊗id(T)) ∘ (spherical(S)⊗id(dual(S))⊗id(W)⊗id(T)) ∘ (inv(a(S,dual(S),W))⊗id(T)) ∘ inv(a(S,dual(S)⊗W,T)) ∘ (id(S)⊗g) for g ∈ basis]
#                 # Correct dual basis to right (co)domain via Hom(U⊗V,W) ≃ Hom(U,U*⊗W)

#                 basis_dual = [(id(dual(T))⊗ev(W)) ∘ a(dual(T),dual(W),W) ∘ (id(dual(T)⊗dual(W))⊗((ev(dual(S))∘(spherical(S)⊗id(dual(S))))⊗id(W))) ∘ (id(dual(T)⊗dual(W))⊗inv(a(S,dual(S),W))) ∘ a(dual(T)⊗dual(W),S,dual(S)⊗W) ∘ (f⊗id(dual(S)⊗W)) for f ∈ basis_dual]

#                 if length(basis) == 0 
#                     γ_i_temp = vertical_dsum(γ_i_temp, zero_morphism(dom_i, W⊗T⊗X⊗dual(T))) 
#                 else
#                     @show "component_iso"
#                     component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗g) ∘ a(S⊗X,dual(S),W) for (f,g) ∈ zip(basis, basis_dual)])
#                     @show "sum"
#                     γ_i_temp = vertical_dsum(γ_i_temp, (component_iso))
#                 end
#             end
#             γ[i] = horizontal_dsum(γ[i], dim(S)*γ_i_temp)
#         end
#         # distribution Before
#         distr_before = distribute_left([s⊗X⊗dual(s) for s ∈ simples],W)
#         # distribution After
#         distr_after = distribute_right(W,[s⊗X⊗dual(s) for s ∈ simples])

#         γ[i] = inv(distr_after) ∘ γ[i] ∘ distr_before
#     end
        
#     return CenterObject(CenterCategory(base_ring(X),parent(X)),Z,γ)
# end


#-------------------------------------------------------------------------------
#   Pairing and Dual Basis
#-------------------------------------------------------------------------------

function pairing(f::Morphism, g::Morphism)
    A,B = domain(f), codomain(g)
    return ev(B)∘(f⊗g)∘coev(A)
end

# function pairing(f::Morphism, g::Morphism,S,W,T)
#     return (id(W)⊗ev(dual(T))) ∘ (id(W)⊗(spherical(T)⊗id(dual(T)))) ∘ associator(W,T,dual(T)) ∘ (f⊗g) ∘ inv(associator(S,dual(S),W)) ∘ (coev(S)⊗id(W))
# end

# function dual_basis(S::Object, W::Object, T::Object)
#     base = basis(Hom(S,W⊗T))
#     base2 = basis(Hom(dual(S)⊗W,dual(T)))
#     dual_base = []
#     F = base_ring(S)
#     n = length(base)
#     for i ∈ 1:n
#         M = matrix(F, n,n, [pairing(base[j], base2[k], S,W,T) for j ∈ 1:n, k ∈ 1:n])
#         b = matrix(F,n,1, [i == j ? 1 : 0 for j ∈ 1:n])
#         coeffs = solve(M,b)
#         push!(dual_base, sum(coeffs .* base2))
#     end
#     @show  M = matrix(F, n,n, [pairing(base[j], dual_base[k], S,W,T) for j ∈ 1:n, k ∈ 1:n])
#     return base, dual_base
# end
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


