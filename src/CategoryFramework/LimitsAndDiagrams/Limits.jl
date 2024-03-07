#=----------------------------------------------------------
    Abstract Cone and Limit types 
----------------------------------------------------------=#

abstract type Cone end
abstract type Limit <: Cone end

# Explicit Limit types

struct Kernel <: Limit 
    object::Object
    morphism::Morphism
end

struct Cokernel <: Limit
    object::Object
    morphism::Morphism
end

struct Equilizer <: Limit
    object::Object
    morphism::Morphism
end

struct Coequilizer <: Limit
    object::Object
    morphism::Morphism
end



#=----------------------------------------------------------
    Limit constructions 
----------------------------------------------------------=#

function _equilizer(f::T, g::T) where T <: Morphism 
    @assert domain(f) == domain(g) && codomain(f) == codomain(g)
    return kernel(f-g)
end

function equilizer(f::Morphism...)
    length(f) == 1 && return f[1]
    length(f) == 2 && return _equilizer(f...)

    g,h = f 

    E,e = _equilizer(g,h)

    f_new = [g ∘ e for g ∈ f[2:end]]

    return equilizer(f_new)
end


function _coequilizer(f::T, g::T) where T <: Morphism 
    @assert domain(f) == domain(g) && codomain(f) == codomain(g)
    return cokernel(f-g)
end

function coequilizer(f::Morphism...)
    length(f) == 1 && return f[1]
    length(f) == 2 && return _coequilizer(f...)

    g,h = f 

    C,c = _coequilizer(g,h)

    f_new = [c ∘ g for g ∈ f[2:end]]

    return coequilizer(f_new)
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

function universal_property_of_pushout(ϕ₁::Morphism, ϕ₂::Morphism, f::Morphism, g::Morphism)
    @assert codomain(f) == codomain(g)
    @assert domain(f) == domain(ϕ₁)
    @assert domain(g) == domain(ϕ₂)

    Z = codomain(f)
    base_f = basis(Hom(domain(ϕ₁), Z))
    base_g = basis(Hom(domain(ϕ₂), Z))

    n = length(base_g) + length(base_f)  
    
    if n == 0 
        return zero_morphism(codomain(ϕ₁), Z)
    end
    
    base_dom_cod = basis(Hom(codomain(ϕ₁), Z))

    # Set up equations to find the unique factoring morphism from the 
    # universal property of the pushout
    Rx, x = PolynomialRing(base_ring(f), length(base_dom_cod))

    eqs = [zero(Rx) for _ ∈ 1:n]

    for (h,a) ∈ zip(base_dom_cod, x)
        e_f = express_in_basis(h ∘ ϕ₁, base_f)
        e_g = express_in_basis(h ∘ ϕ₂, base_g)

        eqs = eqs .+ (a.* [e_f; e_g])
    end
    
    M_arr = hcat([[coeff(e, a) for a ∈ x] for e ∈ eqs]...)
    b_arr = [express_in_basis(f, base_f); express_in_basis(g, base_g)]

    M = matrix(base_ring(f), length(base_dom_cod), length(eqs), M_arr)
    b = matrix(base_ring(f), 1, length(eqs), b_arr)

    s = solve_left(M,b)

    return sum(collect(s)[:] .* base_dom_cod)
end


function pushout_product(f::T, g::T) where T <: Morphism 
    dom, (ϕ₁, ϕ₂) = pushout(f ⊗ id(domain(g)), id(domain(f)) ⊗ g)

    A,B = domain(f), codomain(f)
    X,Y = domain(g), codomain(g)
    Z = B⊗Y

    # @show base_g = basis(Hom(B⊗X, Z))
    #  base_f = basis(Hom(A⊗Y, Z))

    # n = length(base_g) + length(base_f)

    
    universal_property_of_pushout(ϕ₁,ϕ₂, id(B) ⊗ g, f ⊗ id(Y))

    #base_dom_cod = basis(Hom(dom, Z))

    # # Set up equations to find the unique factoring morphism from the 
    # # universal property of the pushout
    # Rx, x = PolynomialRing(base_ring(f), length(base_dom_cod))

    # eqs = [zero(Rx) for _ ∈ 1:n]

    # m_g = id(B) ⊗ g
    # m_f = f ⊗ id(Y)

    # for (h,a) ∈ zip(base_dom_cod, x)
    #     e_g = express_in_basis(h ∘ ϕ₁, base_g)
    #     e_f = express_in_basis(h ∘ ϕ₂, base_f)

    #     eqs = eqs .+ (a.* [e_g; e_f])
    # end
    
    # M_arr = hcat([[coeff(e, a) for a ∈ x] for e ∈ eqs]...)
    # b_arr = [express_in_basis(m_g, base_g); express_in_basis(m_f, base_f)]

    # M = matrix(base_ring(f), length(base_dom_cod), length(eqs), M_arr)
    # b = matrix(base_ring(f), 1, length(eqs), b_arr)

    # s = solve_left(M,b)

    # return sum(collect(s)[:] .* base_dom_cod)
end

⋆(f::T, g::T) where T <: Morphism = pushout_product(f,g)


#=----------------------------------------------------------
    Universal properties 
----------------------------------------------------------=#

function universal_morphism(X::Limit, Y::Cone)
end


#=----------------------------------------------------------
    Coend
----------------------------------------------------------=#

function is_bifunctor(F::AbstractFunctor)
    typeof(domain(F)) != ProductCategory && return false
    C = codomain(F)
    domain(F) == ProductCategory(OppositeCategory(C),C)
end

function coend(F::AbstractFunctor)
    @assert is_bifunctor(F)
    C = codomain(F)
    
    is_abelian(C) && return abelian_coend(F)

    error("not supported yet")
end

function abelian_coend(F::AbstractFunctor)
    C = codomain(F)
    @assert is_abelian(C)
    
    indecs = indecomposables(C)

    _,i,p = direct_sum([F(ProductObject(OppositeObject(X),X)) for X ∈ indecs])

end

