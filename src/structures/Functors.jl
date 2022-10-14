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
    obj_map = x -> X.V
    mor_map = f -> Morphism(domain(f).V, codomain(f).V, f.m)
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
    print(io,"$(typeof(H.domain) == OppositeCategory ? "Contravariant" : "Covariant") Hom-functor in $(H.domain)")
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


#=------------------------------------------------
    Dual
------------------------------------------------=#

# struct DualFunctor <: Functor
#     domain::Category
#     codomain::Category
#     obj_map
#     mor_map
# end

# function DualFunctor(C::Category)
#     @assert ismonoidal(C) "Category has to be monoidal"
#     obj_map = X -> dual(X)
#     mor_map = f -> dual(f)
#     return DualFunctor(C,C,obj_map,mor_map)
# end

function dual_monoidal_structure(X::Object, Y::Object)
    (ev(X⊗Y)⊗id(dual(Y)⊗dual(X))) ∘ inv(associator(dual(X⊗Y),X⊗Y,dual(Y)⊗dual(X))) ∘ (id(dual(X⊗Y))⊗product_coev(X,Y))
end

function product_coev(X::Object, Y::Object)
    associator(X⊗Y,dual(Y),dual(X)) ∘ (inv(associator(X,Y,dual(Y)))⊗id(dual(X))) ∘ (id(X)⊗coev(Y)⊗id(dual(X))) ∘ coev(X)
end

function product_ev(X::Object, Y::Object)
    ev(Y) ∘ (id(dual(Y))⊗ev(X)⊗id(Y)) ∘ (associator(dual(Y),dual(X),X)⊗id(Y)) ∘ inv(associator(dual(Y)⊗dual(X),X,Y))
end


# function show(io::IO, F::DualFunctor)
#     show(io, "Duality Functor in $(F.domain)")
# end
