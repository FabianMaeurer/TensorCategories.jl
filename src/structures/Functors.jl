abstract type Functor end


domain(F::Functor) = F.domain
codomain(F::Functor) = F.codomain

#-------------------------------------------------------------------------------
#   Functors Functionality
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#   Forgetful Functors
#-------------------------------------------------------------------------------

struct Forgetful <: Functor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

function Forgetful(C::GradedVectorSpaces, D::VectorSpaces)
    obj_map = x -> dsum([V for (g,V) in x.V]...)
    mor_map = m -> dsum([f for (g,f) in m.m]...)
    return Forgetful(C,D,obj_map, mor_map)
end

(F::Functor)(x::T) where {T <: Object} = F.obj_map(x)
(F::Functor)(x::T) where {T <: Morphism} = F.mor_map(x)

function show(io::IO, F::Forgetful)
    print(io, "Forgetful functor from $(domain(F)) to $(codomain(F))")
end
#-------------------------------------------------------------------------------
#   Hom Functors
#-------------------------------------------------------------------------------

struct HomFunctor <: Functor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

function Hom(X::Object,::Colon)
    K = base_ring(parent(X))
    C = VectorSpaces(K)
    obj_map = Y -> Hom(X,Y)
    mor_map = f -> g -> g ∘ f
    return HomFunctor(parent(X),C,obj_map,mor_map)
end

function Hom(::Colon,X::Object)
    K = base_ring(parent(X))
    C = VectorSpaces(K)
    obj_map = Y -> Hom(Y,X)
    mor_map = g -> f -> g ∘ f
    return HomFunctor(OppositeCategory(parent(X)),C,obj_map,mor_map)
end

function Hom(X::SetObject,::Colon)
    obj_map = Y -> Hom(X,Y)
    mor_map = f -> g -> g ∘ f
    return HomFunctor(parent(X),Sets(),obj_map,mor_map)
end

function Hom(::Colon,X::SetObject)
    obj_map = Y -> Hom(Y,X)
    mor_map = g -> f -> g ∘ f
    return HomFunctor(parent(X),Sets(),obj_map,mor_map)
end

function show(io::IO, H::HomFunctor)
    print(io,"Hom-functor in $(H.domain)")
end
#-------------------------------------------------------------------------------
#   Tensor Product Functors
#-------------------------------------------------------------------------------

struct TensorFunctor <: Functor
    domain::ProductCategory{2}
    codomain::Category
    obj_map
    mor_map
end

function TensorFunctor(C::Category)
    domain = ProductCategory(C,C)
    obj_map = X -> X[1]⊗X[2]
    mor_map = f -> f[1]⊗f[2]
    return TensorFunctor(domain, C, obj_map, mor_map)
end

⊗(C::Category) = TensorFunctor(C)


#-------------------------------------------------------------------------------
#   Restriction and Induction
#-------------------------------------------------------------------------------

struct GRepRestriction <: Functor
    domain::GroupRepresentationCategory
    codomain::GroupRepresentationCategory
    obj_map
    mor_map
end

function Restriction(C::GroupRepresentationCategory, D::GroupRepresentationCategory)
    @assert base_ring(C) == base_ring(D) "Not compatible"
    #@assert issubgroup(base_group(D), base_group(C))[1] "Not compatible"
    obj_map = X -> restriction(X, base_group(D))
    mor_map = f -> restriction(f, base_group(D))
    return GRepRestriction(C,D,obj_map,mor_map)
end

struct GRepInduction <: Functor
    domain::GroupRepresentationCategory
    codomain::GroupRepresentationCategory
    obj_map
    mor_map
end

function Induction(C::GroupRepresentationCategory, D::GroupRepresentationCategory)
    @assert base_ring(C) == base_ring(D) "Not compatible"
    #@assert issubgroup(base_group(C), base_group(D))[1] "Not compatible"
    obj_map = X -> induction(X, base_group(D))
    mor_map = f -> induction(f, base_group(D))
    return GRepInduction(C,D, obj_map, mor_map)
end

function show(io::IO, F::GRepRestriction)
    print(io,"Restriction functor from $(domain(F)) to $(codomain(F)).")
end

function show(io::IO, F::GRepInduction)
    print(io,"Induction functor from $(domain(F)) to $(codomain(F)).")
end
