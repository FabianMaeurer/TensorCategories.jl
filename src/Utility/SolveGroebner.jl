function recover_solutions(p::Tuple, K::Field, var_order)
    p = p[2]
    f = p.elim
    g = p.denom
    v = p.param
    rs = roots(K,f)
    solutions = []
    

    if length(p.lf_cfs) == 0
        perm = indexin(var_order, p.vars)
        for r ∈ rs
            solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v];r])[perm]]
        end
    else
        perm = indexin(var_order, p.vars[1:end-1])
        for r ∈ rs
            solutions = [solutions; Tuple([vi(r)*inv(g(r)) for vi ∈ v])[perm]]
        end
    end

    if length(rs) < degree(f)
        @info "More solutions over splitting field of $f"
    end

    # for r ∈ rs
    #     solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v];r])[perm]]
    # end

    return solutions
end

function real_solutions_over_base_field(I::MPolyIdeal)
    K = base_ring(base_ring(I))

    if is_finite(K)
        return recover_solutions(real_solutions(I), K, symbols(base_ring(I)))
    end

    if K == QQ
        return recover_solutions(real_solutions(I), K, symbols(base_ring(I)))
    end

    G = gens(I)

    try 
        I = ideal([change_base_ring(QQ, f) for f ∈ gens(I)])
        S = recover_solutions(real_solutions(I), K, symbols(base_ring(I)))
        return [s for s ∈ S if all([g(s...) == 0 for g ∈ G])]
    catch
    end

    QI = rational_lift(I)

    S = recover_solutions(real_solutions(QI), K, symbols(base_ring(QI)))
    
    unique([s[1:end-1] for s ∈ S if all([g(s[1:end-1]...) == 0 for g ∈ G])])
end

function guess_real_solutions_over_base_field(I::MPolyIdeal)
    K = base_ring(base_ring(I))

    d = dim(I)
    G = gens(groebner_basis(I))
    x = gens(base_ring(I))

    i = findall(f -> is_univariate(f), G)
    
    y = [y for y ∈ x if sum([degree(G[j], y) for j ∈ i]) == 0]

    sort!(y, by = x -> sum([degree(g, x) for g ∈ G]))
    
    J = ideal(G)

    while d > 0 && length(y) > 0
        z = pop!(y)

        J2 = ideal([gens(J); z*(z^2 - 1)])

        d2 = dim(J2)
        if d2 ≥ 0 && d2 < d
            J = J2
            d = d2
        end
    end

    # y2 = collect(Base.product(y2,y2))[:]
    # filter!(e -> e[1]==e[2], y2)
    # groebner_basis(J)
    # while d > 0 && length(y2) > 0
    #     z1,z2 = popfirst!(y2)
    #     J2 = ideal([gens(J); [z1*(z1^2 - 1), z2*(z2^2 - 1)]])

    #     d2 = dim(J2)
    #     if d2 ≥ 0 && d2 < d
    #         J = J2
    #         d = d2
    #     end 
    # end

    real_solutions_over_base_field(J)
end
