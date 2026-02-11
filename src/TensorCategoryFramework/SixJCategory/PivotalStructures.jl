#=----------------------------------------------------------
    Compute pivotal structures for SixJCategories 
----------------------------------------------------------=#

function pivotal_structures(C::SixJCategory)
    @assert is_multifusion(C) 

    n = rank(C) 

    K = base_ring(C) 

    Kx, x = polynomial_ring(K,n - sum(C.one))

    piv = [s ∈ one(C).components ? Kx(1) : popfirst!(x) for s ∈ 1:n]

    C_x = six_j_category(Kx, multiplication_table(C))
    set_associator!(C_x, [change_base_ring(Kx, m) for m ∈ C.ass])
    set_one!(C_x, C.one)
    
    set_pivotal!(C_x, piv)

    eqs = []
    for i ∈ 1:n, j ∈ 1:n

        eq = matrix(pivotal(C_x[i]) ⊗ pivotal(C_x[j])) - matrix(double_dual_monoidal_structure(C[i],C[j])) * matrix(pivotal(C_x[i] ⊗ C_x[j]))

        append!(eqs, filter!(!=(0), collect(matrix(eq))[:]))
    end
    
    sols = real_solutions_over_base_field(ideal(collect(groebner_basis(ideal(eqs)))))

    return [[p(s...) for p ∈ piv] for s ∈ sols]
end

spherical_structures(C::SixJCategory) = filter(p -> is_spherical(C,p), pivotal_structures(C))