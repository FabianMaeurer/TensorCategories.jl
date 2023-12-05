#=----------------------------------------------------------
    Limit constructions 
----------------------------------------------------------=#

function equilizer(f::T, g::T) where T <: Morphism 
    @assert domain(f) == domain(g) && codomain(f) == codomain(g)
    return kernel(f-g)
end

function coequilizer(f::T, g::T) where T <: Morphism 
    @assert domain(f) == domain(g) && codomain(f) == codomain(g)
    return cokernel(f-g)
end

function pullback(f::T,g::T) where T <: Morphism
    @assert codomain(f) == codomain(g)

    X = domain(f)
    Y = domain(g)
    Z = codomain(g)

    _, (p_X, p_Y) = product(X,Y)
    
    PullBack, e = equilizer(f ∘ p_X, g ∘ p_Y)
    
    return PullBack, [p_X ∘ e, p_Y ∘ e]
end

function pushout(f::T,g::T) where T <: Morphism
    @assert domain(f) == domain(g)

    X = domain(f)
    Y = codomain(f)
    Z = codomain(g)

    _, (i_Y, i_Z) = coproduct(Y,Z)
    
    PushOut, c = coequilizer(i_Y ∘ f, i_Z ∘ g)
    
    return PushOut, [c ∘ i_Y, c ∘ i_Z]
end

function pushout_product(f::T, g::T) where T <: Morphism
    dom, (ϕ₁, ϕ₂) = pushout(f ⊗ id(domain(g)), id(domain(f)) ⊗ g)

    A,B = domain(f), codomain(f)
    X,Y = domain(g), codomain(g)
    Z = B⊗Y

    base_g = basis(Hom(B⊗X, Z))
    base_f = basis(Hom(A⊗Y, Z))

    n = length(base_g) + length(base_f)

    if n == 0 
        return zero_morphism(dom, Z)
    end
    
    base_dom_cod = basis(Hom(dom, Z))

    # Set up equations to find the unique factoring morphism from the 
    # universal property of the pushout
    Rx, x = PolynomialRing(base_ring(f), length(base_dom_cod))

    eqs = [zero(Rx) for _ ∈ 1:n]

    m_g = id(B) ⊗ g
    m_f = f ⊗ id(Y)

    for (h,a) ∈ zip(base_dom_cod, x)
        e_g = express_in_basis(h ∘ ϕ₁, base_g)
        e_f = express_in_basis(h ∘ ϕ₂, base_f)

        eqs = eqs .+ (a.* [e_g; e_f])
    end
    
    M_arr = hcat([[coeff(e, a) for a ∈ x] for e ∈ eqs]...)
    b_arr = [express_in_basis(m_g, base_g); express_in_basis(m_f, base_f)]

    M = matrix(base_ring(f), length(base_dom_cod), length(eqs), M_arr)
    b = matrix(base_ring(f), 1, length(eqs), b_arr)

    s = solve_left(M,b)

    return sum(collect(s)[:] .* base_dom_cod)
end

⋆(f::T, g::T) where T <: Morphism = pushout_product(f,g)