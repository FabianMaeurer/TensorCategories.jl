#------------------------------------------------------------------------
#   Structs for categories
#------------------------------------------------------------------------


abstract type Category end

abstract type TensorCategory{T <: FieldElem} <: Category end

abstract type Object end

abstract type Morphism end

#abstract type HomSet end

domain(m::Morphism) = m.domain
codomain(m::Morphism) = m.codomain

function dsum(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,ix = dsum(T[1],X)
    m = vcat([ix[1] ∘ t for t in T[2]], ix[2:2])
    return Z, m
end

function dsum(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,ix = dsum(X[1],X[2])
    for Y in X[3:end]
        Z,ix = dsum((Z,ix),Y)
    end
    return Z,ix
end

⊕(T::Tuple{S,Vector{R}},X::S1) where {S <: Object, S1 <: Object, R <: Morphism} = dsum(T,X)
^(X::Object,n::Integer) = dsum([X for i in 1:n]...)
