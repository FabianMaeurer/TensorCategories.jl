abstract type Functor end

domain(F::Functor) = F.domain
codomain(F::Functor) = F.codomain

#-------------------------------------------------------------------------------
#   Forgetful Functors
#-------------------------------------------------------------------------------

struct Forgetful <: Functor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

function Forgetful(C::GradedVectorSpaces{T,G}, D::VectorSpaces{T}) where {T,G}
    obj_map = x -> dsum([V for (g,V) in x.V]...)[1]
    mor_map = m -> dsum([f for (g,f) in m.m]...)
    return Forgetful(C,D,obj_map, mor_map)
end

(F::Functor)(x::T) where {T <: Object} = F.obj_map(x)
(F::Functor)(x::T) where {T <: Morphism} = F.mor_map(x)
