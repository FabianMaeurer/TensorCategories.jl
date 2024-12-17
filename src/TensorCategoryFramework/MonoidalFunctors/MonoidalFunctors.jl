#=----------------------------------------------------------
    Compute Monoidal Functors   
----------------------------------------------------------=#

function monoidal_structures(F::AbstractFunctor)
    C = domain(F)
    @assert is_tensor(C)

    S = simples(C)
    one_C = one(C)
    non_trivial_S = [s for s in S if s != one_C]

    bases = Dict()
    var_numbers = Dict()

    for X ∈ S, Y ∈ S
        if int_dim(Hom(one_C, X)) > 0 || int_dim(Hom(one_C, Y)) > 0
            bases[(X,Y)] = [id(F(X ⊗ Y))]
            var_numbers[(X,Y)] = 0
        else
            bases[(X,Y)] = basis(Hom(F(X)⊗F(Y), F(X⊗Y)))
            var_numbers[(X,Y)] = length(bases[(X,Y)])
        end
    end

    K = base_ring(C)
    R,x = polynomial_ring(K, sum(collect(values(var_numbers))))

    # split variables for bases
    y = copy(x)
    vars = Dict((s,t) => [popfirst!(y) for _ ∈ 1:var_numbers[(s,t)]] for s ∈ S, t ∈ S)
    @show vars

    equations = []
    for (X,Y,Z) ∈ zip(non_trivial_S, non_trivial_S, non_trivial_S)
        eq_basis = basis(Hom(F((X) ⊗ F(Y)) ⊗ F(Z), F(X ⊗ (Y ⊗ Z))))

        # Decompose X⊗Y and Y⊗Z
        _,_,ixy,pxy = direct_sum_decomposition(X ⊗ Y, S)
        _,_,iyz,pyz = direct_sum_decomposition(Y ⊗ Z, S)
        
        # Equations for J_XY,Z ∘ (J_X,Y ⊗ id_Z) = J_X,YZ ∘ (id_X ⊗ J_Y,Z)
        # First iterate over X,Y and Y,Z
        
        for (a, J_XY) ∈ zip(vars[(X,Y)], bases[(X,Y)]), (b, J_YZ) ∈ zip(vars[(Y,Z)], bases[(Y,Z)])
            @show a,b
            # Iterate over factors of XY and YZ 
            eq = [zero(R) for _ ∈ length(eq_basis)]

            for (i,p) ∈ zip(ixy,pxy)
                V =  S[findfirst(==(domain(i)), S)]
                for (c,t) ∈ zip(vars[(V,Z)], bases[(V,Z)])
                    @show c
                    J_XY_Z = (F(i ⊗ id(Z))) ∘ t ∘ (F(p) ⊗ id(Z))
                    coeffs = express_in_basis(compose(
                        J_XY ⊗ id(F(Z)),
                        J_XY_Z,
                        F(associator(X,Y,Z))),
                        eq_basis)

                    eq = eq .+ (a*c) .* coeffs 
                end
            end

            for (i,p) ∈ zip(iyz,pyz)
                V =  S[findfirst(==(domain(i)), S)]
                for (c,t) ∈ zip(vars[(X,V)], bases[(X,V)])
                    J_X_YZ = F(id(X) ⊗ i) ∘ t ∘ (F(id(X))⊗ p)
                    coeffs = express_in_basis(compose(
                        associator(F(X),F(Y),F(Z)),
                        id(F(X)) ⊗ J_YZ,
                        J_X_YZ),
                        eq_basis)

                    eq = eq .- ((b*c) .* coeffs) 
                end
            end

            equations = [equations; eq]
        end
    end

    unique!(equations)
    filter!(!iszero, equations)

    # Equations for invertibility
    iso_mats = [sum([v .* matrix(f) for (v,f) ∈ zip(vars[(x,y)], bases[(x,y)])]) for (x,y) ∈ keys(vars) if !isempty(vars[(x,y)])]

    KR = fraction_field(R)
    inv_iso_mats = [inv(change_base_ring(KR, m)) for m ∈ iso_mats]

    @show [n for (m,n) in zip(iso_mats, inv_iso_mats)]

    # mats = [(dim(x)*dim(y), sum([a .* matrix(f) for (a,f) ∈ zip(vars[(x,y)], bases[(x,y)])])) for (x,y) ∈ keys(vars) if !isempty(vars[(x,y)])]
    # dets = [det(m) - inv(d) for (d,m) ∈ mats]
    # equations = [equations; dets]
    return equations
end