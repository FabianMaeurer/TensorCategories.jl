
#-------------------------------------------------------------------------------
#   Induction
#-------------------------------------------------------------------------------

function induction(X::Object, simples::Vector = simples(parent(X)); parent_category::CenterCategory = center(parent(X)))

    if isdefined(parent_category, :inductions) && X ∈ collect(keys(parent_category.inductions)) 
        return parent_category.inductions[X]
    end

    @assert is_semisimple(parent(X)) "Requires semisimplicity"
    Z = direct_sum([s⊗X⊗dual(s) for s ∈ simples])[1]
    a = associator
    γ = Vector{Morphism}(undef, length(simples))

    C = parent(X)

    @threads for i ∈ 1:length(simples)
        W = simples[i]
        γ[i] = zero_morphism(zero(parent(X)), direct_sum([W⊗((s⊗X)⊗dual(s)) for s ∈ simples])[1])
  
        for S ∈ simples
            dom_i = ((S⊗X)⊗dual(S))⊗W
            γ_i_temp = zero_morphism(dom_i, zero(parent(X)))
            for  T ∈ simples
                #@show S,T,W
                # Set up basis and dual basis
                
                #_basis, basis_dual = dual_basis(Hom(S, W⊗T), Hom(dual(S), dual(W⊗T)))
                
                #_basis, basis_dual = adjusted_dual_basis(Hom(S, W⊗T), Hom(dual(S)⊗W, dual(T)), S, W, T)

                 _,_,ic,p = direct_sum_decomposition(W⊗T, simples)
                _basis = [f for f ∈ ic if domain(f) == S]
                dual_basis = [f for f ∈ p if codomain(f) == S]
                
                #_,_,_,p = direct_sum_decomposition(dual(S)⊗W, dual.(simples))

                if length(_basis) == 0 
                    γ_i_temp = vertical_direct_sum(γ_i_temp, zero_morphism(dom_i, W⊗((T⊗X)⊗dual(T)))) 
                    continue
                end

                basis_dual = transform_dual_basis(dual_basis, S,W,T)

                #_basis, basis_dual = adjusted_dual_basis(_basis, basis_dual, S, W, T)


                #corrections = base_ring(X).([(id(W)⊗ev(dual(T))) ∘ (id(W)⊗(pivotal(T)⊗id(dual(T)))) ∘ a(W,T,dual(T)) ∘ (f⊗g) ∘ a(S,dual(S),W) ∘ (coev(S)⊗id(W)) for (f,g) ∈ zip(_basis,basis_dual)])
                
            
                #@show "component_iso"
                #component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(inv(dim(W))*inv(k)*g)) ∘ a(S⊗X,dual(S),W) for (k,f,g) ∈ zip(corrections,_basis, basis_dual)])

                component_iso = sum([a(W,T⊗X,dual(T)) ∘ (a(W,T,X)⊗id(dual(T))) ∘ ((f⊗id(X))⊗(g)) ∘ a(S⊗X,dual(S),W) for (f,g) ∈ zip(_basis, basis_dual)])

                #@show "sum"
                γ_i_temp = vertical_direct_sum(γ_i_temp, (component_iso))
            end
            γ[i] = horizontal_direct_sum(γ[i], γ_i_temp)
        end
      

        # distribution Before
        distr_before = distribute_left([s⊗X⊗dual(s) for s ∈ simples],W)
        # distribution After
        distr_after = distribute_right(W,[s⊗X⊗dual(s) for s ∈ simples])

        γ[i] = inv(distr_after) ∘ γ[i] ∘ distr_before 
    end

    # if !is_split_semisimple(C)
    #     # factor out such that the left and right tensor product coincide.
    #     r = direct_sum([horizontal_direct_sum([image(f⊗id(X)⊗id(dual(b)) - id(b)⊗id(X)⊗dual(f))[2] for f in End(b)]) for b in simples])

    #     Z,r = cokernel(r)
    #     ir = right_inverse(r)

    #     γ = [(id(b)⊗r) ∘ γᵢ ∘ (ir ⊗ id(b)) for (γᵢ,b) ∈ zip(γ, simples)]
    # end

    IX = CenterObject(parent_category,Z,γ)

    add_induction!(parent_category, X, IX)

    return IX
end


function induction_restriction(X::Object, simples::Vector = simples(parent(X)))
    @assert is_semisimple(parent(X)) "Requires semisimplicity"

    Z = direct_sum([s⊗X⊗dual(s) for s ∈ simples])[1]

    if !is_split_semisimple(parent(X))
        # factor out such that the left and right tensor product coincide.
        r = direct_sum([horizontal_direct_sum([image(f⊗id(X)⊗id(dual(b)) - id(b)⊗id(X)⊗dual(f))[2] for f in End(b)]) for b in simples])

        Z,r = cokernel(r)
    end

    Z
end

function induction_restriction(f::Morphism, simples::Vector = simples(parent(f)))

    ind_f = direct_sum([id(s)⊗f⊗id(dual(s)) for s in simples])

    if !is_split_semisimple(parent(f))
        # factor out such that the left and right tensor product coincide.
        Z,r,ir = induction_non_split_quotient_map(domain(f), simples)

        return r ∘ ind_f ∘ ir
    end

    ind_f
end


function end_of_induction(X::Object, IX = induction(X))
    @assert is_split_semisimple(parent(X))

    B = basis(Hom(X,object(IX)))

    if length(B) == 0 
        return HomSpace(IX, IX, morphism_type(parent(IX))[])
    end
    simpls = simples(parent(X))

    m = [zero_morphism(X,X) for _ in 1:length(simpls)]

    dims = right_dim.((simpls))

    @threads for i ∈ 1:length(simpls)
        xi = simpls[i]
        dxi = dual(xi)
        d = dims[i]
        if typeof(base_ring(X)) == CalciumField #|| base_ring(X) == QQBarField()
            m[i] = d*compose(
                associator(xi, object(IX), dxi),
                id(xi) ⊗ half_braiding(IX, dxi),
                inv_associator(xi, dxi, object(IX)),
                (ev(dxi) ∘ (pivotal(xi) ⊗ id(dxi))) ⊗ id(object(IX))
            ) 
        else
            m[i] = d*compose(
                inv(half_braiding(IX, xi)) ⊗ id(dxi),
                associator(object(IX), xi, dxi),
                id(object(IX)) ⊗ (ev(dxi) ∘ (pivotal(xi) ⊗ id(dxi)))
            ) 
        end
    end

    m = horizontal_direct_sum(m)

    ind_B = [m ∘ induction_restriction(f) for f ∈ B]

    ind_B = [morphism(IX,IX, f) for f ∈ ind_B]

    return HomSpace(IX,IX,ind_B)
end

function induction_adjunction(H::AbstractHomSpace, Y::CenterObject, IX = induction(domain(H), parent_category = parent(Y)))
    @assert is_split_semisimple(parent(H[1]))

    C = parent(domain(H))

    dims = right_dim.((simples(C)))
   
    simpls = simples(parent(H[1]))
    if typeof(base_ring(Y)) == CalciumField || base_ring(Y) == QQBarField()
        ind_f = [d * compose(
            associator(xi, object(Y), dual(xi)),
            id(xi) ⊗ half_braiding(Y, dual(xi)),
            inv_associator(xi, dual(xi), object(Y)),
            (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
        ) for (xi,d) ∈ zip(simpls, dims)]
    else
        ind_f = [d * compose(
            inv(half_braiding(Y, xi)) ⊗ id(dual(xi)),
            associator(object(Y), xi, dual(xi)),
            id(object(Y)) ⊗ (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi))))
        ) for (xi,d) ∈ zip(simpls,dims)]
    end
    # ind_f = [(dim(xi))*compose(
    #     (id(xi) ⊗ f) ⊗ id(dual(xi)),
    #     associator(xi, object(Y), dual(xi)),
    #     id(xi) ⊗ half_braiding(Y, dual(xi)),
    #     inv_associator(xi, dual(xi), object(Y)),
    #     (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
    # ) for xi ∈ simples(parent(f))]

    mors = [morphism(IX, Y, horizontal_direct_sum(ind_f) ∘ induction_restriction(f)) for f ∈ H]

    return HomSpace(IX, Y, mors)
    # morphism(IX,Y, horizontal_direct_sum([sqrt(dim(xi))*((ev(dual(xi)) ∘(pivotal(xi)⊗id(dual(xi))))⊗id(object(IX))) ∘ (id(xi)⊗half_braiding(Y,dual(xi))) ∘ associator(xi,object(Y),dual(xi)) ∘ ((id(xi)⊗f)⊗id(dual(xi))) for xi in simples(parent(X))]))
end

function induction_right_adjunction(H::AbstractHomSpace, Y::CenterObject, IX = induction(codomain(H[1]), parent_category = parent(Y)))

    length(basis(H)) == 0 && return HomSpace(Y, IX, CenterMorphism[])
    simpls = simples(parent(H[1]))

    C = parent(domain(H))


    dims = if is_spherical(C)
        dim.(simples(C))
    else
        sqrt.(squared_norm.(simples(C)))
    end


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
    #     (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
    # ) for xi ∈ simples(parent(f))]

    base = [morphism(Y, IX, induction_restriction(f) ∘ vertical_direct_sum(ind_f)) for f ∈ H]

    HomSpace(Y, IX, base)
end

function induction_right_adjunction(f::Morphism, Y::CenterObject, IX = induction(codomain(f), parent_category = parent(Y)))

    simpls = simples(parent(H[1]))

    dims = if is_spherical(C)
        dim.(simples(C))
    else
        sqrt.(squared_norm.(simples(C)))
    end


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
    #     (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
    # ) for xi ∈ simples(parent(f))]

    morphism(Y, IX, induction_restriction(f) ∘ vertical_direct_sum(ind_f))
end

function induction_adjunction(f::Morphism, Y::CenterObject, IX = induction(domain(f), parent_category = parent(Y)))
    @assert is_split_semisimple(parent(f))

    C = parent(domain(H))

    dims = right_dim.(simples(C))
   
    simpls = simples(parent(H[1]))
    if typeof(base_ring(Y)) == CalciumField || base_ring(Y) == QQBarField()
        ind_f = [d * compose(
            associator(xi, object(Y), dual(xi)),
            id(xi) ⊗ half_braiding(Y, dual(xi)),
            inv_associator(xi, dual(xi), object(Y)),
            (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
        ) for (xi,d) ∈ zip(simpls, dims)]
    else
        ind_f = [d * compose(
            inv(half_braiding(Y, xi)) ⊗ id(dual(xi)),
            associator(object(Y), xi, dual(xi)),
            id(object(Y)) ⊗ (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi))))
        ) for (xi,d) ∈ zip(simpls,dims)]
    end
    # ind_f = [(dim(xi))*compose(
    #     (id(xi) ⊗ f) ⊗ id(dual(xi)),
    #     associator(xi, object(Y), dual(xi)),
    #     id(xi) ⊗ half_braiding(Y, dual(xi)),
    #     inv_associator(xi, dual(xi), object(Y)),
    #     (ev(dual(xi)) ∘ (pivotal(xi) ⊗ id(dual(xi)))) ⊗ id(object(Y))
    # ) for xi ∈ simples(parent(f))]

    mors = morphism(IX, Y, horizontal_direct_sum(ind_f) ∘ induction_restriction(f))

end

function inverse_induction_adjunction(f::CenterMorphism, X::Object)
    @assert object(domain(f)) == induction_restriction(X)

    C = parent(X)
    simpls = simples(C)
    one_indices = [int_dim(Hom(one(C), s)) > 0 for s ∈ simpls]

    P₀ = direct_sum([i == true ? id(s)⊗id(X)⊗id(s) : zero_morphism(zero(C), s⊗X⊗dual(s)) for (i,s) ∈ zip(one_indices, simpls)])

    return morphism(f) ∘ P₀
end



function induction_non_split_quotient_map(X::Object,simples = simples(parent(X)))
    r = direct_sum([horizontal_direct_sum([image(f⊗id(X)⊗id(dual(b)) - id(b)⊗id(X)⊗dual(f))[2] for f in End(b)]) for b in simples])

    Z,r = cokernel(r)
    ir = right_inverse(r)

    Z,r,ir
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
    ϕ = (id(W)⊗ev(dual(T))) ∘ (id(W)⊗(pivotal(T)⊗id(dual(T)))) ∘ associator(W,T,dual(T)) ∘ (f⊗g) ∘ associator(S,dual(S),W) ∘ (coev(S)⊗id(W))

    int_dim(End(W)) == 1 && return dim(W) * ϕ
    
    return tr(ϕ)
end

function adjusted_sum_pairing(f::Morphism, g::Morphism, S::Object, W::Object, T::Object)
    # Compute the pairing for f: S → W⊗T, g:S*⊗W, T*
    # that is defined via the direct sum inclusions
    # and pojections of Hom(S, W⊗T)
    WT = W⊗T
    compose(
        f,
        coev(S) ⊗ id(WT),
        associator(S,dual(S),(WT)),
        id(S) ⊗ inv_associator(dual(S),W,T),  
        id(S) ⊗ compose(
            g ⊗ id(T),
            ev(T)
        )
    )
end


@doc raw""" 

    dual_basis(V::AbstractHomSpace, W::AbstractHomSpace)

Compute the dual basis for Hom(X,Y) and Hom(X̄,Ȳ)
"""
function dual_basis(V::AbstractHomSpace, W::AbstractHomSpace)
    dual_basis = []
    F = base_ring(V)
    n = length(basis(V))
    m = length(basis(W))
    for i ∈ 1:m
        M = matrix(F, n,m, [pairing(basis(V)[j], basis(W)[k]) for j ∈ 1:n, k ∈ 1:m])
        b = matrix(F,m,1, [i == j ? 1 : 0 for j ∈ 1:m])
        coeffs = collect(solve(M,b))
        push!(dual_basis, sum(coeffs .* basis(W)))
    end
    return basis(V), dual_basis
end


@doc raw""" 

    adjusted_dual_basis(V::AbstractHomSpace, U::AbstractHomSpace, S::Object, W::Object, T::Object)

Compute a dual basis for the spaces Hom(S, W⊗T) and Hom(S̄⊗W, T̄)
"""
function adjusted_dual_basis(V::AbstractHomSpace, U::AbstractHomSpace, S::Object, W::Object, T::Object)
    dual_basis = []
    F = base_ring(V)
    n = length(basis(V))
    m = length(basis(U))
    for i ∈ 1:m
        M = matrix(F, n,m, [adjusted_pairing(basis(V)[j], basis(U)[k], S, W, T) for j ∈ 1:n, k ∈ 1:m])
        b = matrix(F,m,1, [i == j ? 1 : 0 for j ∈ 1:m])
        coeffs = collect(solve(M,b))
        push!(dual_basis, sum(coeffs .* basis(U)))
    end
    return basis(V), dual_basis
end

function transform_dual_basis(U::Vector{M}, S::Object, W::Object, T::Object) where M <: Morphism

    length(U) == 0 && return M[] 
    
    dS = dual(S)
    dT = dual(T)

    [compose(
        id(dS) ⊗ (id(W) ⊗ coev(T)),
        id(dS) ⊗ inv_associator(W, T, dT),
        id(dS) ⊗ (g ⊗ id(dT)),
        inv_associator(dS,S,dT),
        ev(S) ⊗ id(dT)
    ) for g ∈ U]
end

function adjusted_dual_basis(V::Vector{<:Morphism}, U::Vector{<:Morphism}, S::Object, W::Object, T::Object)
    dual_basis = []
    F = base_ring(V[1])
    n = length(V)
    m = length(U)
    for i ∈ 1:m
        M = matrix(F, n,m, [adjusted_pairing(V[j], U[k], S, W, T) for j ∈ 1:n, k ∈ 1:m])
        b = matrix(F,m,1, [i == j ? 1 : 0 for j ∈ 1:m])
        coeffs = collect(solve(M,b))
        push!(dual_basis, sum(coeffs .* U))
    end
    return V, dual_basis
end


