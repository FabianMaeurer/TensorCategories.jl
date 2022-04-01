function recover_solutions(p::Tuple, K::Field)
    p = p[1]
    f = p[1]
    g = p[2]
    v = p[3] .* p[4]
    F = splitting_field(f)

    if degree(F) > degree(K) @warn "would split over $F" end
    f = change_base_ring(K,f)
    rs = roots(f)
    solutions = []
    for r ∈ rs
        solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v]; [r]])]
    end
    return solutions
end
