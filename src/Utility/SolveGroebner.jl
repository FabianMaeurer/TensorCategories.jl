function recover_solutions(p::Tuple, K::Field)
    p = p[2]
    f = p.elim
    g = p.denom
    v = p.param
    rs = roots(K,f)
    solutions = []

    

    if length(p.lf_cfs) == 0
        perm = sortperm(p.vars)
        for r ∈ rs
            solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v];r])[perm]]
        end
    else
        perm = sortperm(p.vars[1:end-1])
        for r ∈ rs
            solutions = [solutions; Tuple([vi(r)*inv(g(r)) for vi ∈ v])[perm]]
        end
    end

    # for r ∈ rs
    #     solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v];r])[perm]]
    # end

    return solutions
end

function real_solutions_over_base_field(I::MPolyIdeal)
    K = base_ring(base_ring(I))

    if is_finite(K)
        return recover_solutions(real_solutions(I), K)
    end

    if K == QQ
        return recover_solutions(real_solutions(I), K)
    end

    G = gens(I)

    try 
        I = ideal([change_base_ring(QQ, f) for f ∈ gens(I)])
        S = recover_solutions(real_solutions(I), K)
        return [s for s ∈ S if all([g(s...) == 0 for g ∈ G])]
    catch
    end

    S = recover_solutions(real_solutions(rational_lift(I)), K)
    
    [s[1:end-1] for s ∈ S if all([g(s[1:end-1]...) == 0 for g ∈ G])]
end

function guess_real_solutions_over_base_field(I::MPolyIdeal)
    K = base_ring(base_ring(I))

    d = dim(I)
    G = gens(groebner_basis(I))
    x = gens(base_ring(I))

    i = findall(f -> is_univariate(f), G)
    
    y = [y for y ∈ x if sum([degree(G[j], y) for j ∈ i]) == 0]

    J = I
    while d > 0 && length(y) > 0
        z = pop!(y)
        J2 = ideal([gens(J); z*(z^2 - 1)])

        d2 = dim(J)
        if d2 ≥ 0
            J = J2
            d = d2
        end
    end
    real_solutions_over_base_field(J)
end
