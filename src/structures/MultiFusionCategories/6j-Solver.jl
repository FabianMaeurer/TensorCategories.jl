#=----------------------------------------------------------
    Very crude attempt to solve pentagon equations 
----------------------------------------------------------=#

function six_j_symbols(mult::Array{Int,3}, one::Vector{Int}, K::Field = QQBar)
    C,eqs = pentagon_equations(mult, one)

    sols = recover_solutions(real_solutions(ideal(eqs)), K)
    return [C.ass.(s) for s ∈ sols]
end

function pentagon_equations(mult::Array{Int,3}, one::Vector{Int})
    #Dummy category with the corresponding multiplication to find number of indeterminates
    _C = RingCategory(QQ,mult)
    _C.one = one
    var_count = _number_of_variables_in_pentagon_equations(_C)

    duals = [findfirst( ==(dual(s)), simples(_C)) for s ∈ simples(_C)]
    
    R,x = PolynomialRing(QQ, var_count)

    y = deepcopy(x)
    
    poly_C = RingCategory(R, mult)
    poly_C.one = one

    m = poly_C.simples
    fpdims = fpdim.(simples(_C))

    eqs = elem_type(R)[]

    for i ∈ 1:m, j ∈ 1:m, k ∈ 1:m, o ∈ 1:m
        if sum(poly_C.one[[i,j,k]]) > 0
            continue
        end
        (r,t) = size(poly_C.ass[i,j,k,o])
        for p ∈ 1:r, q ∈ 1:t
            if i == duals[j] == k && p == q == 1
                eqs = [eqs; minpoly(inv(fpdims[i]))(y[end])]
            end
            
            poly_C.ass[i,j,k,o][p,q] = pop!(y)
        end
    end


    for X ∈ simples(poly_C), Y ∈ simples(poly_C), Z ∈ simples(poly_C), W ∈ simples(poly_C)
        f = (id(X)⊗associator(Y,Z,W)) ∘ associator(X,Y⊗Z,W) ∘ (associator(X,Y,Z)⊗id(W))
        g = associator(X,Y,Z⊗W) ∘ associator(X⊗Y,Z,W)

        eqs = [eqs; collect(matrix(f-g))[:]]
    end

    return poly_C, filter(e -> e != 0, unique(eqs))
end

function _number_of_variables_in_pentagon_equations(C::RingCategory)
    n = 0
    S = [s for s ∈ simples(C) if s != one(C)]
    sum(int_dim(End(x ⊗ y ⊗ z)) for x ∈ S, y ∈ S, z ∈ S)
end


function (f::QQMPolyRingElem)(x...)
    coeffs = coefficients(f)
    exps = exponents(f)

    return sum([Rational(c)* *((x.^e)...) for (c,e) ∈ zip(coeffs, exps)])
end