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

#-------------------------------------------------------------------------------
#   Hom Functors
#-------------------------------------------------------------------------------

struct HomFunctor <: Functor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

function Hom(X::Object,v)
    K = base_ring(parent(X))
    C = VectorSpaces(K)
    obj_map = Y -> Hom(X,Y)
    mor_map = f -> g -> g ∘ f
    return HomFunctor(parent(X),C,obj_map,mor_map)
end

function Hom(v,X::Object)
    K = base_ring(parent(X))
    C = VectorSpaces(K)
    obj_map = Y -> Hom(Y,X)
    mor_map = g -> f -> g ∘ f
    return HomFunctor(parent(X),C,obj_map,mor_map)
end

function Hom(X::SetObject,v)
    obj_map = Y -> Hom(X,Y)
    mor_map = f -> g -> g ∘ f
    return HomFunctor(parent(X),Sets(),obj_map,mor_map)
end

function Hom(v,X::SetObject)
    obj_map = Y -> Hom(Y,X)
    mor_map = g -> f -> g ∘ f
    return HomFunctor(parent(X),Sets(),obj_map,mor_map)
end
