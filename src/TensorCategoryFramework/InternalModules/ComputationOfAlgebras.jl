#=----------------------------------------------------------
    Find structures of Algebra objects in fiat categories 
----------------------------------------------------------=#

function algebra_structures(X::Object, unit = Hom(one(parent(X)), X)[1])

    mult_base = basis(Hom(X⊗X, X))
    I = _algebra_structure_ideal(X, mult_base, unit)

    d = dim(I)

    if d < 0 
        return AlgebraObject[]
    elseif d == 0
        sols = real_solutions_over_base_field(I)
    else
        sols = guess_real_solutions_over_base_field(I)
    end

    ms = [sum(s .* mult_base) for s ∈ sols]

    [AlgebraObject(parent(X), X, m, unit) for m ∈ ms]
end

function separable_algebra_structures(X::Object, unit = Hom(one(parent(X)), X)[1])
    [A for A ∈ algebra_structures(X, unit) if is_separable(A)]
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

    for (a, f) ∈ zip(x_m, mult_basis)
        for (a2, f2) ∈ zip(x_m, mult_basis)
            first = compose(ass, id(X) ⊗ f, f2)
            second = compose(f ⊗ id(X), f2)
            coeffs = express_in_basis(first - second, mult_coeff_basis)

            eqs_mult = eqs_mult .+ ((a * a2) .* coeffs)
        end
    
        coeffs = express_in_basis(f ∘ (unit ⊗ id(X)), unit_coeff_basis)
        eqs_unit_l = eqs_unit_l .+ (a .* coeffs) 

        coeffs = express_in_basis(f ∘ (id(X) ⊗ unit), unit_coeff_basis)
        eqs_unit_r = eqs_unit_r .+ (a .* coeffs)
    
    end
    
    id_coeffs = express_in_basis(id(X), unit_coeff_basis)

    eqs_unit_l = eqs_unit_l .- id_coeffs
    eqs_unit_r = eqs_unit_r .- id_coeffs

    return ideal(filter(e -> e != 0, [eqs_mult; eqs_unit_l; eqs_unit_r]))
end
