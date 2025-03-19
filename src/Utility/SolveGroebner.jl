function recover_solutions(p::Tuple, K::Field, var_order, deg = -1; splitting_info = true, fix_solutions = [])
    p = p[2]
    f = p.elim
    g = p.denom
    v = p.param
    rs = roots(K,f)
    solutions = []
 
    if K == QQBarField()
        C = CalciumField()
        rs = C.(rs)
    end

    for (x, val) ∈ fix_solutions
        i = findfirst(==(x), p.vars)
        # find all solutions with xᵢ = v 
        if i > length(p.param) 
            all_i = findall(==((val)), K.(rs))
            rs = rs[all_i]
        else
            all_i = findall(==(val), [K(v[i](r)*inv(g(r))) for r ∈ rs])
            rs = rs[all_i]
        end
    end
  
    if length(p.lf_cfs) == 0
        perm = indexin(var_order, p.vars)
        for r ∈ rs
            solutions = [solutions; Tuple([[K(vi(r)*inv(g(r))) for vi ∈ v];K(r)])[perm]]
        end
    else
        perm = indexin(var_order, p.vars[1:end-1])
        for r ∈ rs
            solutions = [solutions; Tuple([K(vi(r)*inv(g(r))) for vi ∈ v])[perm]]
        end
    end

    if splitting_info && isempty(fix_solutions) && length(rs) < degree(f)
        @info "More solutions over splitting field of $f"
    end

    # for r ∈ rs
    #     solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v];r])[perm]]
    # end
    return solutions
end

function real_solutions_over_base_field(I::MPolyIdeal; splitting_info = true)
    K = base_ring(base_ring(I))

    if is_finite(K)
        return recover_solutions(real_solutions(I), K, symbols(base_ring(I)), splitting_info = splitting_info)
    end

    if K == QQ
        return recover_solutions(real_solutions(I), K, symbols(base_ring(I)),splitting_info = splitting_info)
    end

    G = gens(I)

    try 
        I = ideal([change_base_ring(QQ, f) for f ∈ gens(I)])
        S = recover_solutions(real_solutions(I), K, symbols(base_ring(I)),splitting_info = splitting_info)
        return [s for s ∈ S if all([g(s...) == 0 for g ∈ G])]
    catch
    end

    QI, fixed_sols = rational_lift(I)

    S = recover_solutions(real_solutions(QI), K, symbols(base_ring(QI)), fix_solutions = fixed_sols, splitting_info = splitting_info)
    
    n = dim(base_ring(I))

    unique([s[1:n] for s ∈ S if all([g(s[1:n]...) == 0 for g ∈ G])])
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

function witness_set(I::Ideal, bound = 100)
    d = dim(I)
    deg = degree(I)
    K = base_ring(base_ring(I))
    
    n = dim(base_ring(I))
    
    # a = gen(K)

    # rand_base = [[a^i for i ∈ 0:degree(K)-1]; [0]; [-a^i for i ∈ 0:degree(K)-1]]
    
    QI, fixed_sols = rational_lift(I)
    x = gens(base_ring(QI))[1:n]

    for _ ∈ 1:bound 
        c = matrix(QQ, rand(-1:1,n+1,d)) 

        J = QI + ideal([sum(e .* [x; 1]) for e ∈ eachcol(c)])

        rank(c) < d && continue 
        dim(J) > 0 && continue 
        
        S = recover_solutions(real_solutions(J), K, symbols(base_ring(QI)), fix_solutions = fixed_sols, splitting_info = false) 

        s = unique([s[1:n] for s ∈ S])

        if length(s) > 0
            return s
        end
    end

    error("Too many attempts")
end


