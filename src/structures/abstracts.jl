#------------------------------------------------------------------------
#   Structs for categories
#------------------------------------------------------------------------


abstract type Category end

abstract type TensorCategory{T <: FieldElem} <: Category end

abstract type Object end

abstract type Morphism end

abstract type HomSet end

abstract type HomSpace <: HomSet end

domain(m::Morphism) = m.domain
codomain(m::Morphism) = m.codomain

#---------------------------------------------------------
#   Direct Sums, Products, Coproducts
#---------------------------------------------------------

function ⊕(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,ix = dsum(T[1],X)
    m = vcat([ix[1] ∘ t for t in T[2]], ix[2:2])
    return Z, m
end

⊕(X::S1,T::Tuple{S,Vector{R}}) where {S <: Object,S1 <: Object, R <: Morphism} = ×(T,X)

function dsum(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,ix = ⊕(X[1],X[2])
    for Y in X[3:end]
        Z,ix = ⊕((Z,ix),Y)
    end
    return Z,ix
end


function ×(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = product(T[1],X)
    m = vcat([t ∘ px[1] for t in T[2]], px[2])
    return Z, m
end

×(X::S1,T::Tuple{S,Vector{R}}) where {S <: Object,S1 <: Object, R <: Morphism} = ×(T,X)

function product(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,ix = product(X[1],X[2])
    for Y in X[3:end]
        Z,ix = ×((Z,ix),Y)
    end
    return Z,ix
end

function ∐(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = product(T[1],X)
    m = vcat([px[1] ∘ t for t in T[2]], px[2])
    return Z, m
end

∐(X::S1,T::Tuple{S,Vector{R}}) where {S <: Object,S1 <: Object, R <: Morphism} = ∐(T,X)

function coproduct(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,ix = coproduct(X[1],X[2])
    for Y in X[3:end]
        Z,ix = ∐((Z,ix),Y)
    end
    return Z,ix
end


×(X::Object...) = product(X...)
∐(X::Object...) = coproduct(X...)
⊕(X::Object...) = dsum(X...)

^(X::Object,n::Integer) = dsum([X for i in 1:n]...)

#------------------------------------------------------
#   Abstract Methods
#------------------------------------------------------

issemisimple(C::Category) = :semisimple ∈ features(C)
isabelian(C::Category) = :abelian ∈ features(C)
ismonoidal(C::Category) = :monoidal ∈ features(C)
