#=----------------------------------------------------------
    Compute Monoidal Functors   
----------------------------------------------------------=#

function monoidal_structures(F::AbstractFunctor)
    C = domain(F)
    @assert is_fusion(C)

    S = simples(C)
    n = length(S)

    if n == 1 
        return [monoidal_functor(F, S, Dict(S[1] => id(F)))]
    end

    one_C = one(C)
    one_index = findfirst(==(one_C), S)
    
    #non_trivial_S = [s for s in S if s != one_C]

    bases = Dict()
    var_numbers = Dict()

    for i ∈ 1:n, j ∈ 1:n
        X,Y = S[[i,j]]
        if one_index ∈ [i,j]
            bases[(i,j)] = [id(F(X ⊗ Y))]
            var_numbers[(i,j)] = 0
        else
            bases[(i,j)] = basis(Hom(F(X)⊗F(Y), F(X⊗Y)))
            var_numbers[(i,j)] = length(bases[(i,j)])
        end
    end

    K = base_ring(C)
    R,x = polynomial_ring(K, sum(collect(values(var_numbers))))

    # split variables for bases
    y = copy(x)
    vars = Dict((s,t) => var_numbers[(s,t)] > 0 ? [popfirst!(y) for _ ∈ 1:var_numbers[(s,t)]] : [R(1)] for s ∈ 1:n, t ∈ 1:n)

    equations = []
    for (i,j,k) ∈ Base.product(1:n,1:n,1:n)
        if one_index ∈ [i,j,k] continue end

        X,Y,Z = S[[i,j,k]]

        eq_basis = basis(Hom(F((X) ⊗ F(Y)) ⊗ F(Z), F(X ⊗ (Y ⊗ Z))))
        eq = [zero(R) for _ ∈ length(eq_basis)]

        # Decompose X⊗Y and Y⊗Z
        _,_,ixy,pxy = direct_sum_decomposition(X ⊗ Y, S)
        _,_,iyz,pyz = direct_sum_decomposition(Y ⊗ Z, S)
        
        # Equations for J_XY,Z ∘ (J_X,Y ⊗ id_Z) = J_X,YZ ∘ (id_X ⊗ J_Y,Z)
        # First iterate over X,Y , then Y,Z
        
       
        for (a, J_XY) ∈ zip(vars[(i,j)], bases[(i,j)])
            # Iterate over factors of XY and YZ 

            for (ic,p) ∈ zip(ixy,pxy)
                V = findfirst(==(domain(ic)), S)

                for (c,t) ∈ zip(vars[(V,k)], bases[(V,k)])
                    J_XY_Z = (F(ic ⊗ id(Z))) ∘ t ∘ (F(p) ⊗ id(F(Z)))
                    coeffs = express_in_basis(compose(
                        J_XY ⊗ id(F(Z)),
                        J_XY_Z,
                        F(associator(X,Y,Z))),
                        eq_basis)
                    
                    
                    eq = eq .+ ((a*c) .* coeffs)
                end
               
            end
        end

        for (b, J_YZ) ∈ zip(vars[(j,k)], bases[(j,k)])

            for (ic,p) ∈ zip(iyz,pyz)
                V =  findfirst(==(domain(ic)), S)
                for (c,t) ∈ zip(vars[(i,V)], bases[(i,V)])
                    J_X_YZ = F(id(X) ⊗ ic) ∘ t ∘ (id(F(X)) ⊗ F(p))
                    coeffs = express_in_basis(compose(
                        associator(F(X),F(Y),F(Z)),
                        id(F(X)) ⊗ J_YZ,
                        J_X_YZ),
                        eq_basis)

                    eq = eq .- ((b*c) .* coeffs) 
                end
            end

        end
        equations = [equations; eq]
    end

    unique!(equations)
    filter!(!iszero, equations)

    # Equations for invertibility
    iso_mats = [sum([v .* matrix(f) for (v,f) ∈ zip(vars[(x,y)], bases[(x,y)])]) for (x,y) ∈ keys(vars) if !isempty(vars[(x,y)])]

    KR = fraction_field(R)
    inv_iso_mats = [inv(change_base_ring(KR, m)) for m ∈ iso_mats]

    # Require det = 1
    mats = [sum([a .* matrix(f) for (a,f) ∈ zip(vars[(x,y)], bases[(x,y)])]) for (x,y) ∈ keys(vars) if !isempty(vars[(x,y)])]

    filter!(m -> !is_constant(det(m)), mats) 
    dets = [det(m) for m ∈ mats]


    I = ideal(equations)
    @show dim(I)

    equations = [equations; dets .- 1]



    sols = real_solutions_over_base_field(ideal(equations))

    nats = [Dict((x,y) => sum([v(s...) * t for (v,t) ∈ zip(vars[(x,y)], bases[(x,y)])]) for (x,y) ∈ Base.product(1:n,1:n)) for s ∈ sols]

    [monoidal_functor(F, S, n) for n ∈ nats]
end

function monoidal_natural_transformations(F::AbstractMonoidalFunctor, G::AbstractMonoidalFunctor)

    S = F.indecomposables

    # Natural transformtions between F,G
    Nats = additive_natural_transformations(F.F, G.F, S)


    K = base_ring(domain(F))
    n = length(Nats)
    
    R,x = polynomial_ring(K,n)

    equations = []

    for X ∈ S, Y ∈ S
        eq_basis = basis(Hom(F(X)⊗F(Y), G(X ⊗ Y)))
        eqs = [zero(R) for _ ∈ eq_basis]

        for (a,η) ∈ zip(x, Nats), (b,ν) ∈ zip(x, Nats)

            left = monoidal_structure(G,X,Y) ∘ (η(X) ⊗ ν(Y))  

            coeffs = express_in_basis(left, eq_basis)

            eqs = eqs .+ ((a*b) .* coeffs)
        end

        for (a,η) ∈ zip(x, Nats)
            right = (η(X ⊗ Y)) ∘ monoidal_structure(F,X,Y) 
            coeffs = express_in_basis(right, eq_basis)
            eqs = eqs .+ ((a) .* coeffs)
        end

        equations = [equations; eqs]
    end

    sols = real_solutions_over_base_field(ideal(equations))
    filter!(s -> !all(iszero.(s)), sols)

    mon_nats = [sum(a .* Nats) for a ∈ sols]

    return [AdditiveNaturalTransformation(F,G,S,s.maps) for s ∈ mon_nats]
end

function monoidal_functor_axiom(F::AbstractMonoidalFunctor)
    S = indecomposables(domain(F))
    
    for X ∈ S, Y ∈ S, Z ∈ S 
        left = compose(
            monoidal_structure(F,X,Y) ⊗ id(F(Z)),
            monoidal_structure(F, X ⊗ Y, Z),
            F(associator(X,Y,Z))
        )

        right = compose(
            associator(F(X),F(Y),F(Z)),
            id(F(X)) ⊗ monoidal_structure(F,Y,Z),
            monoidal_structure(F, X, Y ⊗ Z)
        )

        if right != left 
            return false 
        end
    end

    true
end