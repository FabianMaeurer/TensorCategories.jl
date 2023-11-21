
#-------------------------------------------------------------------------------
#   Induction
#-------------------------------------------------------------------------------

function induction(X::Object, simples::Vector = simples(parent(X)); parent_category::CenterCategory = Center(parent(X)))
    @assert is_semisimple(parent(X)) "Requires semisimplicity"
    Z = direct_sum([s⊗X⊗dual(s) for s ∈ simples])[1]
    a = associator
    γ = Vector{Morphism}(undef, length(simples))
    C = parent(X)
    for i ∈ 1:length(simples)
        W = simples[i]
        γ[i] = zero_morphism(zero(parent(X)), direct_sum([W⊗((s⊗X)⊗dual(s)) for s ∈ simples])[1])
  
        for S ∈ simples
            dom_i = ((S⊗X)⊗dual(S))⊗simples[i]
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ simples
                #@show S,T,W
                # Set up basis and dual basis
                
                #_basis, basis_dual = dual_basis(Hom(S, W⊗T), Hom(dual(S), dual(W⊗T)))
                _basis, basis_dual = adjusted_dual_basis(Hom(S, W⊗T), Hom(dual(S)⊗W, dual(T)), S, W, T)
                


                if length(_basis) == 0 
                    γ_i_temp = vertical_direct_sum(γ_i_temp, zero_morphism(dom_i, W⊗((T⊗X)⊗dual(T)))) 
                    continue
                end

                #corrections = base_ring(X).([(id(W)⊗ev(dual(T))) ∘ (id(W)⊗(spherical(T)⊗id(dual(T)))) ∘ a(W,T,dual(T)) ∘ (f⊗g) ∘ a(S,dual(S),W) ∘ (coev(S)⊗id(W)) for (f,g) ∈ zip(_basis,basis_dual)])
                
            
                #@show "component_iso"
                #component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(inv(dim(W))*inv(k)*g)) ∘ a(S⊗X,dual(S),W) for (k,f,g) ∈ zip(corrections,_basis, basis_dual)])

                component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(g)) ∘ a(S⊗X,dual(S),W) for (f,g) ∈ zip(_basis, basis_dual)])

                #@show "sum"
#γ_i_temp = vertical_direct_sum(γ_i_temp, sqrt(dim(T))*(component_iso))
                γ_i_temp = vertical_direct_sum(γ_i_temp, (component_iso))
            end
            #γ[i] = horizontal_direct_sum(γ[i], sqrt(dim(S))*γ_i_temp)
            γ[i] = horizontal_direct_sum(γ[i], (dim(S))*γ_i_temp)

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
    @assert is_semisimple(parent(X)) "Requires semisimplicity"
    Z = direct_sum([s⊗X⊗dual(s) for s ∈ simples])[1]
end


function end_of_induction(X::Object, IX = induction(X))
    B = basis(Hom(X,object(IX)))

    ind_B = [induction_adjunction(f, IX, IX) for f ∈ B]

    return CenterHomSpace(IX,IX,ind_B, VectorSpaces(base_ring(X)))
end

function induction_adjunction(f::Morphism, Y::CenterObject, IX = induction(domain(f)))

    ind_f = [(dim(xi))*compose(
        (id(xi) ⊗ f) ⊗ id(dual(xi)),
        associator(xi, object(Y), dual(xi)),
        id(xi) ⊗ half_braiding(Y, dual(xi)),
        inv_associator(xi, dual(xi), object(Y)),
        (ev(dual(xi)) ∘ (spherical(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
    ) for xi ∈ simples(parent(f))]

    Morphism(IX, Y, horizontal_direct_sum(ind_f))
    # Morphism(IX,Y, horizontal_direct_sum([sqrt(dim(xi))*((ev(dual(xi)) ∘(spherical(xi)⊗id(dual(xi))))⊗id(object(IX))) ∘ (id(xi)⊗half_braiding(Y,dual(xi))) ∘ associator(xi,object(Y),dual(xi)) ∘ ((id(xi)⊗f)⊗id(dual(xi))) for xi in simples(parent(X))]))
end

function inverse_induction_adjunction(f::CenterMorphism, X::Object)
    @assert object(domain(f)) == induction_restriction(X)

    C = parent(X)
    simpls = simples(C)
    one_indices = [int_dim(Hom(one(C), s)) > 0 for s ∈ simpls]

    P₀ = direct_sum([i == true ? id(s)⊗id(X)⊗id(s) : zero_morphism(zero(C), s⊗X⊗dual(s)) for (i,s) ∈ zip(one_indices, simpls)])

    return morphism(f) ∘ P₀
end



#-------------------------------------------------------------------------------
#   Pairing and Dual Basis
#-------------------------------------------------------------------------------

function pairing(f::Morphism, g::Morphism)
    A,B = domain(f), codomain(g)
    return ev(B)∘(f⊗g)∘coev(A)
end

function adjusted_pairing(f::Morphism, g::Morphism, S::Object, W::Object, T::Object)
    # The pairing applied to morphisms 
    # f: S → W⊗T, g:S*⊗W, T*
    # Correspondes to the pairing intriduced in 
    # https://doi.org/10.48550/arXiv.1010.1222 after 
    # natural Isomorphisms
    dim(W)*(id(W)⊗ev(dual(T))) ∘ (id(W)⊗(spherical(T)⊗id(dual(T)))) ∘ associator(W,T,dual(T)) ∘ (f⊗g) ∘ associator(S,dual(S),W) ∘ (coev(S)⊗id(W))
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

function adjusted_dual_basis(V::AbstractHomSpace, U::AbstractHomSpace, S::Object, W::Object, T::Object)
    dual_basis = []
    F = base_ring(V)
    n = length(basis(V))
    m = length(basis(U))
    for i ∈ 1:m
        M = matrix(F, n,m, [adjusted_pairing(basis(V)[j], basis(U)[k], S, W, T) for j ∈ 1:n, k ∈ 1:m])
        b = matrix(F,m,1, [i == j ? 1 : 0 for j ∈ 1:m])
        coeffs = solve(M,b)
        push!(dual_basis, sum(coeffs .* basis(U)))
    end
    return basis(V), dual_basis
end

function partial_induction(X::Object, ind_simples::Vector, simples::Vector = simples(parent(X)); parent_category::CenterCategory = Center(parent(X)))
    @assert is_semisimple(parent(X)) "Requires semisimplicity"
    Z = direct_sum([s⊗X⊗dual(s) for s ∈ ind_simples])[1]
    a = associator
    γ = Vector{Morphism}(undef, length(simples))
    C = parent(X)
    for i ∈ 1:length(simples)
        W = simples[i]
        γ[i] = zero_morphism(zero(parent(X)), direct_sum([W⊗((s⊗X)⊗dual(s)) for s ∈ ind_simples])[1])
  
        for S ∈ ind_simples
            dom_i = ((S⊗X)⊗dual(S))⊗simples[i]
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ ind_simples
                #@show S,T,W
                # Set up basis and dual basis
                
                #basis, basis_dual = dual_basis(Hom(S, W⊗T), Hom(dual(S), dual(W⊗T)))
                _basis, basis_dual = basis(Hom(S, W⊗T)), basis(Hom(dual(S)⊗W, dual(T)))

                if length(_basis) == 0 
                    γ_i_temp = vertical_direct_sum(γ_i_temp, zero_morphism(dom_i, W⊗((T⊗X)⊗dual(T)))) 
                    continue
                end

                corrections = base_ring(X).([(id(W)⊗ev(dual(T))) ∘ (id(W)⊗(spherical(T)⊗id(dual(T)))) ∘ a(W,T,dual(T)) ∘ (f⊗g) ∘ a(S,dual(S),W) ∘ (coev(S)⊗id(W)) for (f,g) ∈ zip(_basis,basis_dual)])
                
                #@show "component_iso"
                component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(inv(dim(W))*inv(k)*g)) ∘ a(S⊗X,dual(S),W) for (k,f,g) ∈ zip(corrections,_basis, basis_dual)])
                #@show "sum"
                γ_i_temp = vertical_direct_sum(γ_i_temp, (component_iso))
            end
            γ[i] = horizontal_direct_sum(γ[i], dim(S)*γ_i_temp)
        end
      

        # distribution Before
        distr_before = distribute_left([s⊗X⊗dual(s) for s ∈ ind_simples],W)
        # distribution After
        distr_after = distribute_right(W,[s⊗X⊗dual(s) for s ∈ ind_simples])

        γ[i] = inv(distr_after) ∘ γ[i] ∘ distr_before 
    end

    return CenterObject(parent_category,Z,γ)
end
