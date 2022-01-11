struct GroupAlgebra{T} <: FreeAlgebra{T}
    basis::GAPGroup
    base_ring::Ring
end

struct GroupAlgebraElem{T} <: AlgebraElem{T}
    parent::GroupAlgebra
    coeffs::Dict{GroupElem,T}
end

#------------------------------------------------------------------
#   Constructors
#------------------------------------------------------------------

function GroupAlgebra(K::Ring,G::GAPGroup)
    @assert iscommutative(K) "Ring is required to be commutative"
    return GroupAlgebra{elem_type(K)}(G,K)
end

getproperty(K::Ring,G::GAPGroup) = GroupAlgebra(K,G)

function GroupAlgebraElem(K::Ring,e::Dict{G,T}) where {G<:GroupElem, T<: RingElem}
    GroupAlgebraElem{T}(K[parent(keys(e)[1])], e)
end

(A::GroupAlgebra)(g::GroupElem) = GroupAlgebraElem(A.base_ring,Dict(g => one(A.base_ring)))

*(c::RingElem,g::GroupElem) = GroupAlgebraElem(parent(c),Dict(g => c))

#------------------------------------------------------------------
#   Arithmetic
#------------------------------------------------------------------

function +(x::GroupAlgebraElem{T},y::GroupAlgebraElem{T}) where T
    coeffs = Dict{GroupElem,T}()
    for g ∈ keys(x.coeffs) ∩ keys(y.coeffs)
        coeffs[g] = x.coeffs[g] + y.coeffs[g]
    end
    for g ∈ setminus(keys(x.coeffs), keys(y.coeffs))
        coeffs[g] = x.coeffs[g]
    end
    for g ∈ setminus(keys(y.coeffs),keys(x.coeffs))
        coeffs[g] = y.coeffs[g]
    end
    return GroupAlgebraElem(x.parent.base_ring,coeffs)
end

function -(x::GroupAlgebraElem{T},y::GroupAlgebraElem{T}) where T
    coeffs = Dict{GroupElem,T}()
    for g ∈ keys(x.coeffs) ∩ keys(y.coeffs)
        coeffs[g] = x.coeffs[g] - y.coeffs[g]
    end
    for g ∈ setminus(keys(x.coeffs), keys(y.coeffs))
        coeffs[g] = x.coeffs[g]
    end
    for g ∈ setminus(keys(y.coeffs),keys(x.coeffs))
        coeffs[g] = -y.coeffs[g]
    end
    return GroupAlgebraElem(x.parent.base_ring,coeffs)
end

zero(A::GroupAlgebra{T}) where T = GroupAlgebraElem(A.base_ring, Dict{elem_type(A.basis),T}())
