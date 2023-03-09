function eigenvalues(M::MatElem{qqbar})
    p = minpoly(M)

    rs = roots(p)
end

function factor(p::PolyElem{qqbar})

    u = leading_coefficient(p)

    rs = [r for r in roots(p)]

    y = gen(parent(p))
    return Fac(parent(p)(u), Dict(y - r => count(==(r), rs) for r ∈ rs))
end

# function roots(p::PolyElem{qqbar})
#     rs = [r for r in roots(rational_lift(p), base_ring(p)) if p(r) == 0]
# end

function roots(p::PolyElem{qqbar})
    coef_deg = maximum([degree(c) for c in coefficients(p)])
    deg = degree(p)
    max_deg = coef_deg * deg
    typeof(max_deg)
    max_height = maximum([height_bits(c) for c in coefficients(p)])
    
    # Set the nessecary prescision to find roots
    prec = 8*(max_height*max_deg)
    CC = ComplexField(prec)

    cp = change_base_ring(CC,p)

    rs = roots(cp, initial_prec = prec)
    return [guess(QQBar, r, max_deg) for r ∈ rs]
end

function rational_lift(p::T) where T <: Union{PolyElem, MPolyElem}
    K = base_ring(p)
    @assert characteristic(K) == 0

    if base_ring(p) == QQ
        return p
    end

    p_coeffs = collect(coefficients(p))
    filter!(e -> e != zero(K), p_coeffs)

    n = length(gens(parent(p)))

    if n == 1
        p_coeffs = reverse(p_coeffs)
    end

    non_rational = [!is_rational(c) for c ∈ p_coeffs]

    Qx,x = PolynomialRing(QQ, Symbol.([["x$i" for i ∈ 1:sum(non_rational)]; parent(p).S]))

    q = zero(Qx)
    l = 1
    y = x[end-n+1:end]
    polys = fmpq_mpoly[]
    for (i,m) in zip(1:lastindex(p_coeffs), monomials(p))
        if non_rational[i]
            w = minpoly(p_coeffs[i])(x[l])
            q = q + x[l]*m(y...)
            l = l+1
            polys = [polys; w]
        else
            typeof(m(y...))
            q = q + roots(minpoly(p_coeffs[i]))[1]*m(y...)
        end
    end

    polys = [polys; q]

    I = ideal(polys)
    
    GB = collect(groebner_basis(I, complete_reduction = true, ordering = lex(Qx)))

    _,x = QQ[String.(parent(p).S)...]

    return GB[1]([[1 for _ ∈ 1:sum(non_rational)]; x]...) 
end


function rational_lift(I::MPolyIdeal)
    ideal(rational_lift.(gens(I)))
end

function *(k::qqbar, p::T) where T <: Union{fmpq_poly, fmpq_mpoly}
    if is_rational(k) 
        return roots(minpoly(k))[1]*p
    else
        error("Cannot promote to common type")
    end
end

function +(k::qqbar, p::T) where T <: Union{fmpq_poly, fmpq_mpoly}
    if is_rational(k) 
        return roots(minpoly(k))[1] + p
    else
        error("Cannot promote to common type")
    end
end

function monomials(f::PolyElem{qqbar})
    x = gen(parent(f))
    return [x^i for i ∈ degree(f):-1:0 if coefficients(f)[i] != 0]
end