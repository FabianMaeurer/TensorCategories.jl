function recover_solutions(p::Tuple, K::Field)
    p = p[2]
    f = p.elim
    g = p.denom
    v = length(p.lf_cfs) == length(p.param) ? p.lf_cfs .* p.param : p.param

    rs = roots(f)
    solutions = []
    for r ∈ rs
        solutions = [solutions; Tuple([[vi(r)*inv(g(r)) for vi ∈ v]; [r]])]
    end
    return solutions
end
