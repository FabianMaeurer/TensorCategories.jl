#------------------------------------------------------------------------
#   Structs for categories
#------------------------------------------------------------------------

abstract type Category end

abstract type Object end

abstract type Morphism end


"""
    VectorSpaceObject

An object in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceObject <: Object end

"""
    VectorSpaceMorphism

A morphism in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceMorphism <: Morphism end

abstract type HomSet end

abstract type HomSpace <: VectorSpaceObject end

domain(m::Morphism) = m.domain
codomain(m::Morphism) = m.codomain

"""
    parent(X::Object)

Return the parent category of the object X.
"""
parent(X::Object) = X.parent

"""
    function parent(f::Morphism)

Return the parent category of ``f``.
"""
parent(f::Morphism) = parent(domain(f))

"""
    base_ring(X::Object)

Return the base ring ```k``` of the ```k```-linear parent category of ```X```.
"""
base_ring(X::Object) = parent(X).base_ring

base_ring(X::Morphism) = parent(domain(X)).base_ring
"""
    base_ring(C::Category)

Return the base ring ```k```of the ```k```-linear category ```C```.
"""
base_ring(C::Category) = C.base_ring

base_group(C::Category) = C.base_group
base_group(X::Object) = parent(X).base_group

#---------------------------------------------------------
#   Direct Sums, Products, Coproducts
#---------------------------------------------------------

function ⊕(T::Tuple{S,Vector{R},Vector{R2}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism, R2 <: Morphism}
    Z,ix,px = dsum(T[1],X,true)
    incl = vcat([ix[1] ∘ t for t in T[2]], ix[2:2])
    proj = vcat([t ∘ px[1] for t in T[3]], px[2:2])
    return Z, incl, proj
end

⊕(X::S1,T::Tuple{S,Vector{R}, Vector{R2}}) where {S <: Object,S1 <: Object, R <: Morphism, R2 <: Morphism} = ⊕(T,X)

function dsum(X::Object...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y ∈ X[2:end]
        Z = dsum(Z,Y)
    end
    return Z
end

function dsum_with_morphisms(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1]),id(X[1])],[id(X[1]),id(X[1])]
    end
    Z,ix,px = dsum(X[1],X[2],true)
    for Y in X[3:end]
        Z,ix,px = ⊕((Z,ix,px),Y)
    end
    return Z,ix,px
end

function dsum(f::Morphism...)
    g = f[1]

    for h ∈ f[2:end]
        g = g ⊕ h
    end
    return g
end

function ×(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = product(T[1],X)
    m = vcat([t ∘ px[1] for t in T[2]], px[2])
    return Z, m
end

×(X::S1,T::Tuple{S,Vector{R}}) where {S <: Object,S1 <: Object, R <: Morphism} = ×(T,X)

function product(X::Object...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y ∈ X[2:end]
        Z = product(Z,Y)
    end
    return Z
end

function product_with_morphisms(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,px = product(X[1],X[2], true)
    for Y in X[3:end]
        Z,px = ×((Z,px),Y)
    end
    return Z,px
end

function ∐(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = coproduct(T[1],X)
    m = vcat([px[1] ∘ t for t in T[2]], px[2])
    return Z, m
end

∐(X::S1,T::Tuple{S,Vector{R}}) where {S <: Object,S1 <: Object, R <: Morphism} = ∐(T,X)

function coproduct(X::Object...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y in X[2:end]
        Z = coproduct(Z,Y)
    end
    return Z
end

function coproduct_with_morphisms(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,ix = coproduct(X[1],X[2])
    for Y in X[3:end]
        Z,ix = ∐((Z,ix),Y)
    end
    return Z,ix
end

"""
    ×(X::Object...)

Return the product Object and an array containing the projection morphisms.
"""
×(X::Object...) = product(X...)

"""
    ∐(X::Object...)

Return the coproduct Object and an array containing the injection morphisms.
"""
∐(X::Object...) = coproduct(X...)

"""
    ⊕(X::Object...)

Return the direct sum Object and arrays containing the injection and projection
morphisms.
"""

⊕(X::Object...) = dsum(X...)

⊕(X::Morphism...) = dsum(X...)

"""
    ⊗(X::Object...)

Return the tensor product object.
"""
⊗(X::Object...) = tensor_product(X...)

"""
    ^(X::Object, n::Integer)

Return the n-fold product object ```X^n```.
"""
^(X::Object,n::Integer) = n == 0 ? zero(parent(X)) : product([X for i in 1:n]...)

^(X::Morphism,n::Integer) = n == 0 ? zero_morphism(zero(parent(domain(X))), zero(parent(domain(X)))) : dsum([X for i in 1:n]...)
"""
    ⊗(f::Morphism, g::Morphism)

Return the tensor product morphism of ```f```and ```g```.
"""
⊗(f::Morphism, g::Morphism) where {T} = tensor_product(f,g)


dsum(X::T) where T <: Union{Vector,Tuple} = dsum(X...)
product(X::T) where T <: Union{Vector,Tuple} = product(X...)
coproduct(X::T) where T <: Union{Vector,Tuple} = coproduct(X...)

product(X::Object,Y::Object) = dsum(X,Y)
coproduct(X::Object, Y::Object) = dsum(X,Y)

#---------------------------------------------------------
#   Horizontal and Vertical direct sums
#---------------------------------------------------------

"""
    function horizontal_dsum(f::Morphism, g::Morphism)

Return the sum of ``f:X → Z``, ``g:Y → Z`` as ``f+g:X⊕Y → Z.
"""
function horizontal_dsum(f::Morphism, g::Morphism)
    #@assert codomain(f) == codomain(g) "Codomains do not coincide"

    sum = f ⊕ g
    _,_,(p1,p2) = dsum_with_morphisms(codomain(f),codomain(g))
    return p1∘sum + p2∘sum
end

"""
    function vertical_dsum(f::Morphism, g::Morphism)

Return the sum of ``f:X → Y``, ``g:X → Z`` as ``f+g: X → Y⊕Z.
"""
function vertical_dsum(f::Morphism, g::Morphism)
    #@assert domain(f) == domain(g) "Domains do not coincide"

    sum = f ⊕ g
    _,(i1,i2),_ = dsum_with_morphisms(domain(f), domain(g))
    return sum∘i1 + sum∘i2
end
#---------------------------------------------------------
#   tensor_product
#---------------------------------------------------------


function tensor_product(X::Object...)
    if length(X) == 1 return X end

    Z = X[1]
    for Y ∈ X[2:end]
        Z = Z⊗Y
    end
    return Z
end

tensor_product(X::T) where T <: Union{Vector,Tuple} = tensor_product(X...)
#------------------------------------------------------
#   Abstract Methods
#------------------------------------------------------
isfusion(C::Category) = false
ismultifusion(C::Category) = isfusion(C)

istensor(C::Category) = isfusion(C)
ismultitensor(C::Category) = ismultifusion(C) || istensor(C)

isring(C::Category) = istensor(C)
ismultiring(C::Category) = ismultitensor(C)

ismonoidal(C::Category) = ismultitensor(C)

isabelian(C::Category) = ismultitensor(C)

isadditive(C::Category) = isabelian(C)

islinear(C::Category) = isabelian(C)

issemisimple(C::Category) = ismultitensor(C)

function image(f::Morphism)
    C,c = cokernel(f)
    return kernel(c)
end

∘(f::Morphism...) = compose(reverse(f)...)

-(f::Morphism, g::Morphism) = f + (-1)*g
-(f::Morphism) = (-1)*f

#-------------------------------------------------------
# Hom Spaces
#-------------------------------------------------------

dim(V::HomSpace) = length(basis(V))

End(X::Object) = Hom(X,X)

zero_morphism(C::Category) = zero_morphism(zero(C), zero(C))

#-------------------------------------------------------
# Duals
#-------------------------------------------------------

left_dual(X::Object) = dual(X)
right_dual(X::Object) = dual(X)

dual(f::Morphism) = left_dual(f)

function left_dual(f::Morphism)
    X = domain(f)
    Y = codomain(f)
    a = ev(Y)⊗id(dual(X))
    b = id(dual(Y)⊗f)⊗id(dual(X))
    c = inv(associator(dual(Y),X,dual(X)))
    d = id(dual(Y)⊗coev(X))
    (a)∘(b)∘(c)∘(d)
end

tr(f::Morphism) = left_trace(f)

function left_trace(f::Morphism)
    V = domain(f)
    W = codomain(f)
    C = parent(V)
    if V == zero(C) || W == zero(C) return zero_morphism(one(C),one(C)) end

    if V == W
        return ev(left_dual(V)) ∘ ((spherical(V)∘f) ⊗ id(left_dual(V))) ∘ coev(V)
    end
    return ev(left_dual(V)) ∘ (f ⊗ id(left_dual(V))) ∘ coev(V)
end

function right_trace(f::Morphism)
    V = domain(f)
    W = codomain(f)
    dV = right_dual(V)
    _,i = isisomorphic(left_dual(dV),V)
    _,j = isisomorphic(right_dual(V), left_dual(right_dual(dV)))
    return (ev(right_dual(dV))) ∘ (j⊗(f∘i)) ∘ coev(right_dual(V))
end

#-------------------------------------------------------
# Spherical structure
#-------------------------------------------------------

function drinfeld_morphism(X::Object)
     (ev(X)⊗id(dual(dual(X)))) ∘ (braiding(X,dual(X))⊗id(dual(dual(X)))) ∘ (id(X)⊗coev(dual(X)))
 end

dim(X::Object) = base_ring(X)(tr(spherical(X)))

dim(C::Category) = sum(dim(s)^2 for s ∈ simples(C))
#-------------------------------------------------------
# S-Matrix
#-------------------------------------------------------

function smatrix(C::Category, simples = simples(C))
    @assert issemisimple(C) "Category has to be semisimple"
    F = base_ring(C)
    m = [tr(braiding(s,t)∘braiding(t,s)) for s ∈ simples, t ∈ simples]
    try
        return matrix(F,[F(n) for n ∈ m])
    catch
    end
    return matrix(F,m)
end

#-------------------------------------------------------
# decomposition morphism
#-------------------------------------------------------

function decompose(X::Object, S = simples(parent(X)))
    C = parent(X)
    @assert issemisimple(C) "Category not semisimple"
    dimensions = [dim(Hom(X,s)) for s ∈ S]
    return [(s,d) for (s,d) ∈ zip(S,dimensions) if d > 0]
end

function decompose_morphism(X::Object)
    C = parent(X)
    @assert issemisimple(C) "Semisimplicity required"

    S = simples(C)
    components = sort( decompose(X), by = e -> findfirst(s -> isisomorphic(s,e[1])[1], S))
    Z, incl, proj = dsum_with_morphisms([s^d for (s,d) ∈ components]...)

    # temporary solution!
    return isisomorphic(X,Z)[2]

    #----------------------------------
    f = zero_morphism(X,Z)

    for (p,i) ∈ zip(proj, incl)
        g = i∘p
        f = f + g
    end
    return f
end




#-------------------------------------------------------
# Semisimple: Subobjects
#-------------------------------------------------------

function eigenspaces(f::Morphism)
    @assert domain(f) == codomain(f) "Not an endomorphism"

    values = collect(keys(eigenspaces(matrix(f))))

    return Dict(λ => kernel(f-λ*id(domain(f)))[1] for λ ∈ values)
end

function irreducible_subobjects(X::Object)
    B = basis(End(X))

    if length(B) == 1 return [X] end

    for f ∈ B
        eig_spaces = eigenspaces(f)

        if length(eig_spaces) == 1 
            continue
        end

        simple_subs = vcat([irreducible_subobjects(K) for (_,K) ∈ eig_spaces]...)
        return unique_simples(simple_subs)
    end
    return [X]
end

function unique_simples(simples::Vector{<:Object})
    unique_simples = simples[1:1]
    for s ∈ simples[2:end]
        if sum([dim(Hom(s,u)) for u ∈ unique_simples]) == 0
            unique_simples = [unique_simples; s]
        end
    end
    return unique_simples
end


#-------------------------------------------------------
# Misc
#-------------------------------------------------------

*(f::Morphism, x) = x*f