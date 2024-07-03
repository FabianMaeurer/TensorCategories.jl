function relative_induction(X::Object, ind_simples::Vector{<:Object}, simples::Vector = simples(parent(X)); parent_category::CentralizerCategory = centralizer(parent(X), ind_simples))

    if isdefined(parent_category, :inductions) && X ∈ collect(keys(parent_category.inductions)) 
        return parent_category.inductions[X]
    end

    ind_simples = parent_category.subcategory_simples

    @assert is_semisimple(parent(X)) "Requires semisimplicity"
    Z = direct_sum([s⊗X⊗dual(s) for s ∈ ind_simples])[1]
    a = associator
    γ = Vector{Morphism}(undef, length(ind_simples))

    C = parent(X)

    for i ∈ 1:length(ind_simples)
        W = ind_simples[i]
        γ[i] = zero_morphism(zero(parent(X)), direct_sum([W⊗((s⊗X)⊗dual(s)) for s ∈ ind_simples])[1])
  
        for S ∈ ind_simples
            dom_i = ((S⊗X)⊗dual(S))⊗W
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ ind_simples
                #@show S,T,W
                # Set up basis and dual basis
                
                #_basis, basis_dual = dual_basis(Hom(S, W⊗T), Hom(dual(S), dual(W⊗T)))
                
                #_basis, basis_dual = adjusted_dual_basis(Hom(S, W⊗T), Hom(dual(S)⊗W, dual(T)), S, W, T)

                 _,_,ic,p = direct_sum_decomposition(W⊗T, simples)
                _basis = [f for f ∈ ic if domain(f) == S]
                dual_basis = [f for f ∈ p if codomain(f) == S]
                #_,_,_,p = direct_sum_decomposition(dual(S)⊗W, dual.(simples))

                basis_dual = transform_dual_basis(dual_basis, S,W,T)

                if length(_basis) == 0 
                    γ_i_temp = vertical_direct_sum(γ_i_temp, zero_morphism(dom_i, W⊗((T⊗X)⊗dual(T)))) 
                    continue
                end


                #_basis, basis_dual = adjusted_dual_basis(_basis, basis_dual, S, W, T)


                #corrections = base_ring(X).([(id(W)⊗ev(dual(T))) ∘ (id(W)⊗(spherical(T)⊗id(dual(T)))) ∘ a(W,T,dual(T)) ∘ (f⊗g) ∘ a(S,dual(S),W) ∘ (coev(S)⊗id(W)) for (f,g) ∈ zip(_basis,basis_dual)])
                
            
                #@show "component_iso"
                #component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(inv(dim(W))*inv(k)*g)) ∘ a(S⊗X,dual(S),W) for (k,f,g) ∈ zip(corrections,_basis, basis_dual)])

                component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(g)) ∘ a(S⊗X,dual(S),W) for (f,g) ∈ zip(_basis, basis_dual)])

                #@show "sum"
                γ_i_temp = vertical_direct_sum(γ_i_temp, (component_iso))
            end
            γ[i] = horizontal_direct_sum(γ[i], γ_i_temp)
        end
      

        # distribution Before
        distr_before = distribute_left([s⊗X⊗dual(s) for s ∈ ind_simples],W)
        # distribution After
        distr_after = distribute_right(W,[s⊗X⊗dual(s) for s ∈ ind_simples])

        γ[i] = inv(distr_after) ∘ γ[i] ∘ distr_before 
    end

    if !is_split_semisimple(C)
        # factor out such that the left and right tensor product coincide.
        r = direct_sum([horizontal_direct_sum([image(f⊗id(X)⊗id(dual(b)) - id(b)⊗id(X)⊗dual(f))[2] for f in End(b)]) for b in simples])

        @show Z,r = cokernel(r)
        ir = right_inverse(r)

        γ = [(id(b)⊗r) ∘ γᵢ ∘ (ir ⊗ id(b)) for (γᵢ,b) ∈ zip(γ, simples)]
    end

    IX = CentralizerObject(parent_category,Z,γ)

    add_induction!(parent_category, X, IX)

    return IX
end

function end_of_induction(X::Object, simpls::Vector{<:Object}, IX = relative_induction(X, S))
    @assert is_split_semisimple(parent(X))

    B = basis(Hom(X,object(IX)))
    
    m = [zero_morphism(X,X) for _ in 1:length(simpls)]
    
    @threads for i ∈ 1:length(simpls)
        xi = simpls[i]
        dxi = dual(xi)
        m[i] = (dim(xi))*compose(
        inv(half_braiding(IX, xi)) ⊗ id(dxi),
        associator(object(IX), xi, dual(xi)),
        id(object(IX)) ⊗ (ev(dxi) ∘ (spherical(xi) ⊗ id(dxi)))
        ) 
    end

    m = horizontal_direct_sum(m)

    ind_B = [m ∘ induction_restriction(f, simpls) for f ∈ B]

    ind_B = [morphism(IX,IX, f) for f ∈ ind_B]

    return HomSpace(IX,IX,ind_B)
end

function induction_right_adjunction(H::AbstractHomSpace, Y::CentralizerObject, IX = relative_induction(codomain(H[1]), parent(Y).subcategory_simples, parent_category = parent(Y)))

    simpls = parent(Y).subcategory_simples

    duals = dual.(simpls)

    ind_f = [compose(
        id(object(Y))⊗coev(xi),
        inv_associator(object(Y),xi,dxi),
        half_braiding(Y,xi) ⊗ id(dxi)
    ) for (xi, dxi) ∈ zip(simpls,duals)]

    # ind_f = [(dim(xi))*compose(
    #     (id(xi) ⊗ f) ⊗ id(dual(xi)),
    #     associator(xi, object(Y), dual(xi)),
    #     id(xi) ⊗ half_braiding(Y, dual(xi)),
    #     inv_associator(xi, dual(xi), object(Y)),
    #     (ev(dual(xi)) ∘ (spherical(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
    # ) for xi ∈ simples(parent(f))]

    base = [morphism(Y, IX, induction_restriction(f, simpls) ∘ vertical_direct_sum(ind_f)) for f ∈ H]

    HomSpace(Y, IX, base)
end

function induction_adjunction(H::AbstractHomSpace, Y::CentralizerObject, IX = relative_induction(domain(H), parent(Y).subcategory_simples, parent_category = parent(Y)))
    @assert is_split_semisimple(parent(H[1]))

    simpls = parent(Y).subcategory_simples

    ind_f = [dim(xi)*compose(
        inv(half_braiding(Y, xi)) ⊗ id(dual(xi)),
        associator(object(Y), xi, dual(xi)),
        id(object(Y)) ⊗ (ev(dual(xi)) ∘ (spherical(xi) ⊗ id(dual(xi))))
    ) for xi ∈ simpls]

    # ind_f = [(dim(xi))*compose(
    #     (id(xi) ⊗ f) ⊗ id(dual(xi)),
    #     associator(xi, object(Y), dual(xi)),
    #     id(xi) ⊗ half_braiding(Y, dual(xi)),
    #     inv_associator(xi, dual(xi), object(Y)),
    #     (ev(dual(xi)) ∘ (spherical(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
    # ) for xi ∈ simples(parent(f))]

    mors = [morphism(IX, Y, horizontal_direct_sum(ind_f) ∘ induction_restriction(f,simpls)) for f ∈ H]

    return HomSpace(IX, Y, mors)
    # morphism(IX,Y, horizontal_direct_sum([sqrt(dim(xi))*((ev(dual(xi)) ∘(spherical(xi)⊗id(dual(xi))))⊗id(object(IX))) ∘ (id(xi)⊗half_braiding(Y,dual(xi))) ∘ associator(xi,object(Y),dual(xi)) ∘ ((id(xi)⊗f)⊗id(dual(xi))) for xi in simples(parent(X))]))
end