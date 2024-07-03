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
    _C = six_j_category(QQ,mult)
    _C.one = one
    var_count = _number_of_variables_in_pentagon_equations(_C)

    duals = [findfirst( ==(dual(s)), simples(_C)) for s ∈ simples(_C)]
    
    R,x = polynomial_ring(QQ, var_count)

    y = deepcopy(x)
    
    poly_C = six_j_category(R, mult)
    poly_C.one = one

    m = poly_C.simples
    fpdims = fpdim.(simples(_C))

    eqs = elem_type(R)[]

    for i ∈ 1:m, j ∈ 1:m, k ∈ 1:m, o ∈ 1:m
        if sum(poly_C.one[[i,j,k]]) > 0
            continue
        end
        (r,t) = size(poly_C.ass[i,j,k,o])
            
        poly_C.ass[i,j,k,o] = matrix(R,r,t,[pop!(y) for _ ∈ 1:r, _ ∈ 1:t])

    end


    for X ∈ simples(poly_C), Y ∈ simples(poly_C), Z ∈ simples(poly_C), W ∈ simples(poly_C)
        f = (id(X)⊗associator(Y,Z,W)) ∘ associator(X,Y⊗Z,W) ∘ (associator(X,Y,Z)⊗id(W))
        g = associator(X,Y,Z⊗W) ∘ associator(X⊗Y,Z,W)

        eqs = [eqs; collect(matrix(f-g))[:]]
    end

    return poly_C, filter(e -> e != 0, unique(eqs))
end

function _number_of_variables_in_pentagon_equations(C::SixJCategory)
    n = 0
    S = [s for s ∈ simples(C) if s != one(C)]
    sum(int_dim(End(x ⊗ y ⊗ z)) for x ∈ S, y ∈ S, z ∈ S)
end


function (f::QQMPolyRingElem)(x...)
    coeffs = coefficients(f)
    exps = exponents(f)

    return sum([Rational(c)* *((x.^e)...) for (c,e) ∈ zip(coeffs, exps)])
end

function split_ideal(I::MPolyIdeal)
    polys = filter(e -> e != 0, unique!(gens(I)))
    R = base_ring(I)
    vars = gens(R)
    n = length(vars)
    var_names = symbols(R)

   
    A = [polys[1]]
    not_taken = deepcopy(polys)[2:end]
    contained = [degree(polys[1],i) > 0 for i ∈ 1:n]
    i = 1
    while true
        f = not_taken[i]
        if sum(contained .* [degree(f,j) for j ∈ 1:n]) > 0
            A = [A;f]
            deleteat!(not_taken, i)
            i = 1
            contained[[degree(f,i) > 0 for i ∈ 1:n]] .= true
        else
            i = i+1
        end
        if i == length(not_taken) + 1
            break
        end
    end

    if length(not_taken) == 0
        return [ideal(A)]
    end
    
    var_names_in_A = var_names[contained]

    R2,y = polynomial_ring(base_ring(R), var_names_in_A)
    x = [c ? popfirst!(y) : 0 for c ∈ contained]

    var_names_in_rest = var_names[contained .⊻ true]

    R3,y = polynomial_ring(base_ring(R), var_names_in_rest)
    x2 = [!c ? popfirst!(y) : 0 for c ∈ contained]

    return [ideal([f(x...) for f ∈ A]); split_ideal(ideal([f(x2...) for f ∈ not_taken]))]
end


solve_fewnomial_system(I::Ideal) = _solve_fewnomial_system(I,1)

function _solve_fewnomial_system(I::Ideal, i::Int = 1)
    B = gens(I)
    n = nvars(base_ring(I))
    y = gens(base_ring(I))

    if n == i 
        if length(B) == 0 
            return [(1,)]
        end
        R,x = polynomial_ring(QQ,:x)
        C = [f([ones(Int,n-1);x]...) for f ∈ B]
        rs = roots(gcd(C))
        return [(r,) for r ∈ rs]
    end

    Bi = B[degree.(B,i) .> 0]
    Ci = filter(r -> r ∉ Bi, B)
    f = popfirst!(Bi)
    @show Bi = resultant.(f,Bi,i)
    unique!(Bi)
    filter!(r -> r ≠ 0, Bi)

    sols = _solve_fewnomial_system(ideal([Bi;Ci]), i+1)

    S = []
    R,x = QQBar[:x]
    B = [change_base_ring(QQBar, f) for f ∈ B]
    for s ∈ sols

        D = [f([R(1) for _ ∈ 1:i-1]...,x,R.(s)...) for f ∈ B]
        filter!(r -> r ≠ 0 ,D)
        unique!(D)
    
        if length(D) == 0 
            S = [S;(QQBar(1),s...)]
        else
            rs = roots(gcd(D))
            S = [S;[(r,s...) for r ∈ rs]]
        end
    end
    return S
end

