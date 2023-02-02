function eigenvalues(M::MatElem{qqbar})
    p = minpoly(M)

    rs = unique(roots(p))
end

function factor(p::Poly{qqbar})

    u = leading_coefficient(p)
    p = inv(u)*p

    p_coeffs = collect(coefficients(p))

    non_rational = [!is_rational(c) for c ∈ p_coeffs]

    Qx,x = PolynomialRing(QQ, sum(non_rational) + 1)

    q = zero(Qx)
    l = 1
    polys = fmpq_mpoly[]
    for i in 1:lastindex(p_coeffs)
        if non_rational[i]
            w = minpoly(p_coeffs[i])(x[l])
            q = q + x[l]*x[end]^(i-1)
            l = l+1
            polys = [polys; w]
        else
            q = q + roots(minpoly(p_coeffs[i]))[1]*x[end]^(i-1)
        end
    end

    polys = [polys; q]

    I = ideal(polys)

    _,y = QQ["x"]
    p2 = collect(groebner_basis(I, complete_reduction = true, ordering = lex(Qx)))[1]([[1 for i ∈ 1:length(x)-1]; y]...)

    rs = [r for r in roots(p2, parent(p_coeffs[1])) if p(r) == 0]

    y = gen(parent(p))
    return Fac(parent(p)(u), Dict(y - r => count(==(r), rs) for r ∈ rs))
end

function roots(p::Poly{qqbar})
    p_coeffs = collect(coefficients(p))

    non_rational = [!is_rational(c) for c ∈ p_coeffs]

    Qx,x = PolynomialRing(QQ, sum(non_rational) + 1)

    q = zero(Qx)
    l = 1
    polys = fmpq_mpoly[]
    for i in 1:lastindex(p_coeffs)
        if non_rational[i]
            w = minpoly(p_coeffs[i])(x[l])
            q = q + x[l]*x[end]^(i-1)
            l = l+1
            polys = [polys; w]
        else
            q = q + roots(minpoly(p_coeffs[i]))[1]*x[end]^(i-1)
        end
    end

    polys = [polys; q]

    I = ideal(polys)

    _,y = QQ["x"]
    p2 = collect(groebner_basis(I, complete_reduction = true, ordering = lex(Qx)))[1]([[1 for i ∈ 1:length(x)-1]; y]...)

    rs = [r for r in roots(p2, parent(p_coeffs[1])) if p(r) == 0]
end