function eigenvalues(M::MatElem{QQBarFieldElem}) 
    p = minpoly(M)

    rs = roots(p)
end

function eigenvalues(M::MatElem{CalciumFieldElem})
    m = change_base_ring(QQBarField(),M)
    base_ring(M).(eigenvalues(m))
end

function factor(p::PolyRingElem{T}) where T <: Union{CalciumFieldElem}

    u = leading_coefficient(p)

    rs = [r for r in roots(p)]

    y = gen(parent(p))
    return Fac(parent(p)(u), Dict(y - r => count(==(r), rs) for r ∈ rs))
end

function factor(p::Poly{<:FlintLocalFieldElem})
    rs = roots(p)
    u = leading_coefficient(p)

    if length(rs) == 0 
        return Fac(parent(p)(u), Dict(p => 1))
    end

    y = gen(parent(p))

    if length(rs) == degree(p)
        return Fac(parent(p)(u), Dict(y - r => count(==(r), rs) for r ∈ rs))
    end

    residue = divexact(p, *([y - r for r ∈ rs]...))

    return Fac(parent(p)(u), Dict([residue; [y - r => count(==(r), rs) for r ∈ rs]]...))
end



# function roots(p::PolyRingElem{QQBarFieldElem})
#     rs = [r for r in roots(rational_lift(p), base_ring(p)) if p(r) == 0]
# end

function roots(p::PolyRingElem{CalciumFieldElem})
    rs = roots(change_base_ring(QQBarField(), p))
    return base_ring(p).(rs)
end

# function roots(p::PolyRingElem{QQBarFieldElem})
#     coef_deg = maximum([degree(c) for c in coefficients(p)])
#     deg = degree(p)
#     max_deg = coef_deg * deg
#     typeof(max_deg)
#     max_height = maximum([height_bits(c) for c in coefficients(p)])
    
#     # Set the nessecary prescision to find roots
#     prec = 8*(max_height*max_deg)
#     prec == 0 ? 2^62 - 1 : prec
#     CC = AcbField(prec == 0 ? typemax(Int) >> 4 - 1 : prec)

#     cp = change_base_ring(CC,p)

#     rs = roots(cp, initial_prec = prec)

#     QQBar_roots = QQBarFieldElem[]

#     for r ∈ rs
#         try 
#             push!(QQBar_roots, guess(QQBarField(), r, max_deg))
#         catch
#             @warn "Maximal precision not high enough to recover root $r"
#         end
#     end
#     return QQBar_roots 
# end

function rational_lift(p::T) where T <: Union{PolyRingElem, MPolyRingElem}
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

    Qx,x = polynomial_ring(QQ, Symbol.([["x$i" for i ∈ 1:sum(non_rational)]; parent(p).S]))

    q = zero(Qx)

    l = 1
    y = x[end-n+1:end]
    polys = QQMPolyRingElem[]
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
    K = base_ring(base_ring((I)))

    #@assert typeof(K) <: Union{QQField, AbsSimpleNumField}
    if typeof(K) == CalciumField || K == QQBarField() 
        return qqbar_rational_lift(I)
    end

    if is_finite(K)
        return finite_field_rational_lift(I)
    end

    if K == QQ 
        return I,[]
    end

    G = [g for g ∈ groebner_basis(I) if g != 0]

    α = gen(K)
    μ = minpoly(α)

    Qx, x = polynomial_ring(QQ, length(gens(base_ring(I))) + 1) 
    base = [μ(x[end])]

    for f ∈ G
        coeffs = [polynomial(QQ, coefficients(a))(x[end]) for a ∈ coefficients(f)]
        monos = [change_base_ring(QQ,m)(x[1:end-1]...) for m ∈ monomials(f)]

        push!(base, sum(coeffs .* monos))
    end

    return ideal(base), [symbols(Qx)[end] => α]
end

function qqbar_rational_lift(I::MPolyIdeal)
    G = [g for g ∈ groebner_basis(I) if g != 0]
    n = length(gens(base_ring(I)))

    # Collect all non rational coefficients
    coeffs = vcat(collect.(coefficients.(G))...)
    filter!(! is_rational, coeffs)
    unique!(coeffs)

    # Add a variable for each root
    Qx, x = polynomial_ring(QQ, n + length(coeffs)) 

    coeffs_to_variable = Dict(a => z for (a,z) ∈ zip(coeffs, x[n+1:end]))
    minpolys = [minpoly(QQBarField()(c))(z) for (c,z) ∈ zip(coeffs, x[n+1:end])]
    
    base = []
    
    # Substitute coefficients by variables
    for f ∈ G
        local coeffs = [is_rational(a) ? QQ(a) : coeffs_to_variable[a] for a ∈ coefficients(f)]
        monos = [change_base_ring(QQ,m)(x[1:n]...) for m ∈ monomials(f)]

        push!(base, sum(coeffs .* monos))
    end

    return ideal([base; minpolys]), [symbols(Qx)[n+i] => c for (i,c) ∈ zip(1:length(coeffs), coeffs)]
end

function finite_field_rational_lift(I::MPolyIdeal)
    K = base_ring(base_ring(I))

    if dim(base_ring(I)) == 1 
        return ideal([change_base_ring(QQ, f) for f ∈ gens(I)]), []
    end

    α = gen(K)
    μ = minpoly(α)

    F = prime_field(K)
    Qx, x = polynomial_ring(QQ, length(gens(base_ring(I))) + 1) 
    base = [change_base_ring(QQ,μ)(x[end])]

    G = groebner_basis(I)

    for f ∈ G
        coeffs = [polynomial(QQ, QQ.([lift(ZZ,c) for c ∈ absolute_coordinates(a)]))(x[end]) for a ∈ coefficients(f)]

        monos = [change_base_ring(QQ,m)(x[1:end-1]...) for m ∈ monomials(f)]

        push!(base, sum(coeffs .* monos))
    end

    return ideal(base), [symbols(Qx)[end] => α]
end

function change_base_ring(K::QQField, p::FqMPolyRingElem)
    c = coefficients(p)

    _,vars = polynomial_ring(QQ, dim(parent(p)))

    return sum(QQ(lift(ZZ, c)) * sum(vars .* e) for (c,m,e) ∈ zip(coefficients(p), monomials(p), exponents(p)))
end

function change_base_ring(K::QQField, p::FqPolyRingElem)
    c = coefficients(p)

    _,x = QQ[:x]

    return sum(QQ(lift(ZZ, c)) * x^i for (c,i) ∈ zip(coefficients(p), 0:degree(p)))
end


function *(k::QQBarFieldElem, p::T) where T <: Union{QQPolyRingElem, QQMPolyRingElem}
    if is_rational(k) 
        return roots(minpoly(k))[1]*p
    else
        error("Cannot promote to common type")
    end
end

function +(k::QQBarFieldElem, p::T) where T <: Union{QQPolyRingElem, QQMPolyRingElem}
    if is_rational(k) 
        return roots(minpoly(k))[1] + p
    else
        error("Cannot promote to common type")
    end
end

function monomials(f::PolyRingElem{QQBarFieldElem})
    x = gen(parent(f))
    return [x^i for i ∈ degree(f):-1:0 if coefficients(f)[i] != 0]
end

function (::QQField)(x::QQBarFieldElem) 
    @assert is_rational(x)

    return roots(minpoly(x))[1]
end

function rand(::QQBarField, I::UnitRange)
    return QQBarField()(rand(ZZ, I))
end