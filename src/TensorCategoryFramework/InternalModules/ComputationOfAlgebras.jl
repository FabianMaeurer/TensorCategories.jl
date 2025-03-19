#=----------------------------------------------------------
    Find structures of Algebra objects in fiat categories 
----------------------------------------------------------=#

@doc raw""" 

    algebra_structures(X::Object)
    algebra_structures(X::Object, unit::Morphism)

Return a set of algebra objects over ``X``. An empty array is returned only if there are no algebra structures. If the algebr is not connected, i.e. ``Hom(𝟙,X) ≠ k``, then a unit should be provided.
"""
function algebra_structures(X::Object, unit = Hom(one(parent(X)), X)[1]; show_dimension = false)
    _algebra_structures(_algebra_structure_ideal, X, unit, show_dimension = show_dimension)
end

function algebra_extensions(A::AlgebraObject, X::BiModuleObject)
    AX, (inclusion,), (projection,) = direct_sum(object(A),object(X))

    u = inclusion ∘ unit(A)

    _algebra_structures(_algebra_structure_ideal, AX, u, subalgebra = A, extension = X, inclusion = inclusion, projection = projection, show_dimension = true)
end

@doc raw""" 

    separable_algebra_structures(X::Object)
    separable_algebra_structures(X::Object, unit::Morphism)

Return a set of separable algebra objects over ``X``. An empty array is returned only if there are no algebra structures. If the algebr is not connected, i.e. ``Hom(𝟙,X) ≠ k``, then a unit should be provided.
"""
function separable_algebra_structures(X::Object, unit = Hom(one(parent(X)), X)[1]; show_dimension = false)
    [A for A ∈ algebra_structures(X, unit, show_dimension = show_dimension) if is_separable(A)]

    # _algebra_structures(_separable_algebra_structure_ideal, X, unit, show_dimension = show_dimension)
end

@doc raw""" 

    commutative_algebra_structures(X::Object)
    commutative_algebra_structures(X::Object, unit::Morphism)

Return a set of commutative algebra objects over ``X``. An empty array is returned only if there are no algebra structures. If the algebr is not connected, i.e. ``Hom(𝟙,X) ≠ k``, then a unit should be provided.
"""
function commutative_algebra_structures(X::Object, unit = Hom(one(parent(X)), X)[1]; show_dimension = false)

    _algebra_structures(_commutative_algebra_structure_ideal, X, unit, show_dimension = show_dimension)
end 

@doc raw""" 

    etale_algebra_structures(X::Object)
    etale_algebra_structures(X::Object, unit::Morphism)

Return a set of separable algebra objects over ``X``. An empty array is returned only if there are no algebra structures. If the algebr is not connected, i.e. ``Hom(𝟙,X) ≠ k``, then a unit should be provided.
"""
function etale_algebra_structures(X::Object, unit = Hom(one(parent(X)), X)[1]; show_dimension = false)
    [A for A ∈ commutative_algebra_structures(X, unit, show_dimension = show_dimension) if is_separable(A)]
end

function _algebra_structures(structure_ideal::Function, X::Object, _unit = Hom(one(parent(X)), X)[1]; show_dimension = false, subalgebra = nothing, extension = nothing, inclusion = nothing, projection = nothing)

    mult_base = basis(Hom(X⊗X, X))

    if subalgebra !== nothing 
        @assert inclusion ∘ unit(subalgebra) == _unit 

        subalgebra_multiplication = compose(
            projection ⊗ projection,
            multiplication(subalgebra),
            inclusion
        )

        # Restrict to A-module morphisms 
        AX = bimodule(subalgebra) ⊕ extension

        # Get bimodule tensor prodct and projection
        AX_AX, c = bimodule_tensor_product(AX, AX)
        
        mult_base = [m ∘ c for m ∈ morphism.(basis(Hom(AX_AX, AX)))]

        coeffs = express_in_basis(subalgebra_multiplication, mult_base)

        free_basis = coeffs .== 0
        mult_base = [subalgebra_multiplication; mult_base[free_basis]]
    end

    m = length(mult_base)
    I = structure_ideal(X, mult_base, _unit)

    d = dim(I)

    if d > 0
        # get rid of isomorphic algebras

        # get a fixed solution for ϕ ∘ unit = unit together with
        # a basis for the linear part 
        if subalgebra === nothing 
            global iso_fixed_part, iso_var_basis = fix_subalgebra(basis(Hom(X,X)), _unit)
        else
            global iso_fixed_part, iso_var_basis = fix_subalgebra(basis(Hom(X,X)), inclusion)
        end

        n = length(iso_var_basis)
        K = base_ring(X)
        Kx, x = polynomial_ring(K, m + n)
        Q = fraction_field(Kx)

        # variables for the multiplication 
        m_vars = x[1:m]

        #variables for the isomorphism
        iso_vars = x[m+1:m+n]

        # Change base_ring of matrices
        iso_var_basis_mats = [change_base_ring(Kx, matrix(f)) for f ∈ iso_var_basis]

        # set up the matrix of ϕ 
        phi_mat = change_base_ring(Kx, matrix(iso_fixed_part)) + 
            sum([a * m for (a,m) ∈ zip(iso_vars, iso_var_basis_mats)])
        phi_mat = change_base_ring(Q, phi_mat)

        # set up matrix of ϕ⊗ϕ 
        phi_squared_mat = matrix(iso_fixed_part ⊗ iso_fixed_part) +
            sum([a .* matrix(iso_fixed_part ⊗ f) for (a,f) ∈ zip(iso_vars, iso_var_basis)]) + 
            sum([a .* matrix(f ⊗ iso_fixed_part) for (a,f) ∈ zip(iso_vars, iso_var_basis)]) + 
            sum([(a*b) .* matrix(f⊗g) for (a,f) ∈ zip(iso_vars, iso_var_basis), (b,g) ∈ zip(iso_vars, iso_var_basis)])
        phi_mat_squared = change_base_ring(Q, phi_squared_mat)

        # set up matrix for m 
        mult_mat = sum([a .* matrix(f) for (a,f) ∈ zip(m_vars, mult_base)])
        mult_mat = change_base_ring(Q, mult_mat)

        # get coefficients of the image multiplication
        image_mult = phi_squared_mat * mult_mat * inv(phi_mat)
        image_coeffs = express_in_basis(morphism(image_mult), morphism.(matrix.(mult_base)))

        # Find a coefficient that is linear in a for every a in iso_vars
        free_indices = []
        for a ∈ iso_vars 
            i = findall(c -> c ∉ free_indices && degree(numerator(image_coeffs[c]), a) ≥ 1, 1:length(image_coeffs)) 

            i = argmin(k -> degree(numerator(image_coeffs[k]), a), i)

            i !== nothing && push!(free_indices, (i,a)) 
        end

        # set free coefficients to 1
        y = gens(base_ring(I))
       
        free_coeffs = [y[i]*(y[i] - 1) for (i,_) ∈ free_indices[1:d]]
        
        I = ideal([gens(I); free_coeffs])
    end

    d = dim(I)

    show_dimension && @info "Dimension of solution set: $(d)"
    
    if d < 0 
        return AlgebraObject[]
    elseif d == 0
        sols = real_solutions_over_base_field(I)
    else
        sols = witness_set(I)
    end

    unique!(sols)

    length(sols) == 0 && error("Algebras exist but non found")

    ms = [sum(s .* mult_base) for s ∈ sols]

    [AlgebraObject(parent(X), X, m, _unit) for m ∈ ms]
end

function fix_unit(base::Vector{<:Morphism}, unit::Morphism)
    # Find a subset of basis such that the unit is fixed
    K = base_ring(parent(unit))
    n = length(base)
    Kx,vars = polynomial_ring(K, n)

    unit_basis = basis(Hom(domain(unit), codomain(unit)))

    eqs = [zero(Kx) for _ ∈ unit_basis]
    for (a,f) ∈ zip(vars, base) 
        eqs = eqs .+ (a .* express_in_basis(f ∘ unit, unit_basis))
    end

    # coefficients of unit
    m = length(unit_basis)
    unit_coeffs = matrix(K, m, 1, express_in_basis(unit, unit_basis))

    # extract coeffs as matrix 
    M = matrix(K, n, m, vcat([[coeff(e,a) for a ∈ vars] for e ∈ eqs]...))

    sol = solve(M,unit_coeffs)
    _,nullsp = nullspace(transpose(M))

    fixed_sol = sum(collect(sol)[:] .* base)
    
    var_sol = [sum(c .* base) for c ∈ eachcol(collect(nullsp))]

    return fixed_sol, var_sol
end

function fix_subalgebra(base::Vector{<:Morphism}, inclusion::Morphism)
    K = base_ring(parent(inclusion))
    n = length(base)
    Kx,vars = polynomial_ring(K, n)

    _basis = basis(Hom(domain(inclusion), codomain(inclusion)))

    eqs = [zero(Kx) for _ ∈ _basis]
    for (a,f) ∈ zip(vars, base) 
        eqs = eqs .+ (a .* express_in_basis(f ∘ inclusion, _basis))
    end

    # coefficients of inclusion
    m = length(_basis)
    _coeffs = matrix(K, m, 1, express_in_basis(inclusion, _basis))

    # extract coeffs as matrix 
    M = matrix(K, n, m, vcat([[coeff(e,a) for a ∈ vars] for e ∈ eqs]...))

    sol = solve(transpose(M),_coeffs, side = :right)
    _,nullsp = nullspace(transpose(M))

    fixed_sol = sum(collect(sol)[:] .* base)
    
    var_sol = [sum(c .* base) for c ∈ eachcol(collect(nullsp))]

    return fixed_sol, var_sol
end

function _algebra_structure_ideal(X::Object, mult_basis::Vector{<:Morphism},  unit::Morphism)
    C = parent(X)
    K = base_ring(C)

    mult_coeff_basis = basis(Hom((X⊗X)⊗X, X))
    unit_coeff_basis = basis(End(X))

    m = length(mult_basis)
    Kx,x_m = polynomial_ring(K, m)

    eqs_mult = [zero(Kx) for _ ∈ mult_coeff_basis]
    eqs_unit_l = [zero(Kx) for _ ∈ unit_coeff_basis]
    eqs_unit_r = [zero(Kx) for _ ∈ unit_coeff_basis]

    ass = associator(X,X,X)

    mult_Hom =  HomSpace(domain(mult_basis[1]), X, mult_basis)

    id_times_multbase = id(X) ⊗ mult_Hom 
    multbase_times_id = mult_Hom ⊗ id(X)

    u_id = (unit ⊗ id(X))
    id_u = (id(X) ⊗ unit)
    for (a, f, id_f, f_id) ∈ zip(x_m, mult_basis, id_times_multbase, multbase_times_id)
        
        for (a2, f2) ∈ zip(x_m, mult_basis)
            first = compose(ass, id_f, f2)
            second = compose(f_id, f2)
            coeffs = express_in_basis(first - second, mult_coeff_basis)

            eqs_mult = eqs_mult .+ ((a * a2) .* coeffs)
        end
    
        coeffs = express_in_basis(f ∘ u_id, unit_coeff_basis)
        eqs_unit_l = eqs_unit_l .+ (a .* coeffs) 

        coeffs = express_in_basis(f ∘ id_u, unit_coeff_basis)
        eqs_unit_r = eqs_unit_r .+ (a .* coeffs)
    
    end
    
    id_coeffs = express_in_basis(id(X), unit_coeff_basis)

    eqs_unit_l = eqs_unit_l .- id_coeffs
    eqs_unit_r = eqs_unit_r .- id_coeffs

    return ideal(filter(e -> e != 0, [eqs_mult; eqs_unit_l; eqs_unit_r]))
end

function _separable_algebra_structure_ideal(X::Object, mult_basis::Vector{<:Morphism},  unit::Morphism)

    J = _algebra_structure_ideal(X, mult_basis, unit)
    ideal([gens(J); non_degenerate_condition(X, mult_basis, gens(base_ring(J)), unit)])
end

function non_degenerate_condition(A::Object, mult_basis::Vector{<:Morphism}, vars::Vector, unit::Morphism)
    
    # An algebra is non-degenerate if there is a certain isomorphism
    # A → A∗. Reference https://doi.org/10.1016/S0550-3213(02)00744-7
    # That also means A is a frobenius Algebra

    C = parent(A)
    K = base_ring(C)

    m = length(mult_basis)
    Kx = parent(vars[1])

    dA = dual(A)
    mat = change_base_ring(Kx, matrix(zero_morphism(A,dA)))
    
    e = dim(A) * left_inverse(unit)

    before = compose(
        id(A) ⊗ coev(A),
        inv_associator(A,A,dA)
    )

    for (a,f) ∈ zip(vars, mult_basis)

        non_degenerate_condition = compose(
            before, 
            (e ∘ f) ⊗ id(dA)
        ) 
        mat = mat + a .* change_base_ring(Kx,matrix(non_degenerate_condition))
    end

    # Extract coeffificients for ϵ: A∗ → A 
    QKx = fraction_field(Kx)
    quo_mat = change_base_ring(QKx, mat)
    inv_quo_mat = inv(quo_mat)

    A_dual_to_A = basis(Hom(dA, A))

    coeffs = express_in_basis(morphism(inv_quo_mat), morphism.(matrix.(A_dual_to_A)))

    # set up comultiplication Δ: A → A⊗A 

    comult_mat = change_base_ring(QKx, matrix(zero_morphism(A, A ⊗ A)))
    ass = associator(A,dA,A)
    for (a, m) ∈ zip(vars, mult_basis)
        after = id(A) ⊗ m 
        for (b,ϕ) ∈ zip(coeffs, A_dual_to_A)
            eq = compose( 
                coev(A) ⊗ id(A),
                ass,
                id(A) ⊗ (ϕ ⊗ id(A)),
                after
            )

            comult_mat = comult_mat + a*b*change_base_ring(QKx,matrix(eq))
        end
    end

    comult_basis = basis(Hom(A, A ⊗ A))
    comult_coeffs = express_in_basis(morphism(comult_mat), morphism.(matrix.(comult_basis)))

    # Add equations for m ∘ Δ = id 
    eqs = sum(a*b*change_base_ring(QKx, matrix(m ∘ d)) for  (a,m) ∈ zip(vars, mult_basis), (b,d) ∈ zip(comult_coeffs, comult_basis)) .- matrix(id(A))

    number_of_new_variables = sum(!isconstant.(denominator.(eqs)))
    
    R,x = polymomial_ring(K, )
    eqs = [numerator(denominator(e) * e) for e ∈ eqs]
    return collect(eqs)[:]
end

function _commutative_algebra_structure_ideal(X::Object, mult_base::Vector{<:Morphism}, unit::Morphism)

    I = _algebra_structure_ideal(X,mult_base,unit)

    vars = gens(base_ring(I))
    
    eqs = [zero(base_ring(I)) for _ ∈ mult_base]

    braid = braiding(X,X)

    # Commutative algebras satisfy m ∘ c_{X,X} = m 
    for (x,f) ∈ zip(vars, mult_base)
        coeffs = express_in_basis(f ∘ braid - f, mult_base)
        eqs = eqs .+ (x .* coeffs)
    end

    ideal([gens(I); eqs])
end