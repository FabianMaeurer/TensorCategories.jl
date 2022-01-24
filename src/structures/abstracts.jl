#------------------------------------------------------------------------
#   Structs for categories
#------------------------------------------------------------------------


abstract type Category end

abstract type MultiTensorCategory{T} <: Category end

abstract type TensorCategory{T} <: MultiTensorCategory{T} end

abstract type MultiFusionCategory{T} <: MultiTensorCategory{T} end

abstract type FusionCategory{T} <: MultiFusionCategory{T} end

abstract type Object end

abstract type Morphism end


"""
    VectorSpaceObject{T}

An object in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceObject{T} <: Object end

"""
    VectorSpaceMorphism{T}()

A morphism in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceMorphism{T} <: Morphism end

abstract type HomSet end

abstract type HomSpace{T} <: VectorSpaceObject{T} end

domain(m::Morphism) = m.domain
codomain(m::Morphism) = m.codomain

"""
    parent(X::Object)

Return the parent category of the object X.
"""
parent(X::Object) = X.parent

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
#----------------------------------------------------------------------
#   Algebras
#----------------------------------------------------------------------

abstract type Algebra <: Ring end
abstract type FreeAlgebra{T} <: Algebra end
abstract type AlgebraQuo <: Algebra end


abstract type HopfAlgebra <: Ring end

abstract type AlgebraElem{T} <: RingElem end
abstract type HopfAlgebraElem{T} <: AlgebraElem{T} end


#---------------------------------------------------------
#   Direct Sums, Products, Coproducts
#---------------------------------------------------------

function ⊕(T::Tuple{S,Vector{R},Vector{R2}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism, R2 <: Morphism}
    Z,ix,px = dsum(T[1],X)
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

function dsum_morphisms(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1]),id(X[1])]
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

function product_morphisms(X::Object...)
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

function coproduct_morphisms(X::Object...)
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
^(X::Object,n::Integer) = product([X for i in 1:n]...)

"""
    ⊗(f::Morphism, g::Morphism)

Return the tensor product morphism of ```f```and ```g```.
"""
⊗(f::Morphism, g::Morphism) where {T} = tensor_product(f,g)


dsum(X::T) where T <: Union{Vector,Tuple} = dsum(X...)
product(X::T) where T <: Union{Vector,Tuple} = product(X...)
coproduct(X::T) where T <: Union{Vector,Tuple} = coproduct(X...)

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

issemisimple(C::Category) = :semisimple ∈ features(C)
isabelian(C::Category) = :abelian ∈ features(C)
ismonoidal(C::Category) = :monoidal ∈ features(C)


∘(f::Morphism...) = compose(reverse(f)...)

#-------------------------------------------------------
# Hom Spaces
#-------------------------------------------------------

dim(V::HomSpace) = length(basis(V))

End(X::Object) = Hom(X,X)


#-------------------------------------------------------
# Duals
#-------------------------------------------------------

left_dual(X::Object) = dual(X)
right_dual(X::Object) = dual(X)

tr(f::Morphism) = left_trace(f)

function left_trace(f::Morphism)
    V = domain(f)
    W = codomain(f)
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


 #-------------------------------------------------------
 # S-Matrix
 #-------------------------------------------------------

function smatrix(C::MultiTensorCategory, simples = simples(C))
    @assert issemisimple(C) "Category has to be semisimple"
    m = [tr(braiding(s,t)∘braiding(t,s)) for s ∈ simples, t ∈ simples]
    try
        F = base_ring(C)
        return [F(n) for n ∈ m]
    catch
        return m
    end
end
