function recover_solutions(p::Tuple, K::Field)
    p = p[2]
    f = p.elim
    g = p.denom
    v = length(p.lf_cfs) == length(p.param) ? p.lf_cfs .* p.param : p.param

    rs = roots(K,f)
    solutions = []

    for r ∈ rs
        solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v]; r])]
    end
    return solutions
end

function real_solutions_over_base_field(I::MPolyIdeal)
    K = base_ring(base_ring(I))

    if K == QQ
        return recover_solutions(real_solutions(I), K)
    end


    S = recover_solutions(real_solutions(rational_lift(I)), K)
    
    G = gens(I)

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
