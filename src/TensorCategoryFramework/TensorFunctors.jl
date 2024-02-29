
#-------------------------------------------------------------------------------
#   Tensor Product Functors
#-------------------------------------------------------------------------------

struct TensorProductFunctor <: AbstractFunctor
    domain::ProductCategory{2}
    codomain::Category
    obj_map
    mor_map
end

function TensorProductFunctor(C::Category)
    domain = ProductCategory(C,C)
    obj_map = X -> X[1]⊗X[2]
    mor_map = f -> f[1]⊗f[2]
    return TensorProductFunctor(domain, C, obj_map, mor_map)
end

⊗(C::Category) = TensorProductFunctor(C)

#=----------------------------------------------------------
    Functor  X ⊗ - : C → C 
----------------------------------------------------------=#

struct LeftTensorProductFunctor <: AbstractFunctor
    domain::Category
    codomain::Category 
    object::Object
end

(F::LeftTensorProductFunctor)(X::Object) = F.object ⊗ X
(F::LeftTensorProductFunctor)(f::Morphism) = id(F.object) ⊗ f

function LeftTensorProductFunctor(X::Object) 
    LeftTensorProductFunctor(parent(X),parent(X), X)
end

⊗(X::Object, ::typeof(-)) = LeftTensorProductFunctor(X)

is_additive(T::LeftTensorProductFunctor) = is_additive(domain(T))

#=----------------------------------------------------------
    Functor  - ⊗ X : C → C 
----------------------------------------------------------=#

struct RightTensorProductFunctor <: AbstractFunctor
    domain::Category
    codomain::Category
    object::Object 
end

(F::RightTensorProductFunctor)(X::Object) = X ⊗ F.object 
(F::RightTensorProductFunctor)(f::Morphism) = f ⊗ id(F.object)

function RightTensorProductFunctor(X::Object) 
    RightTensorProductFunctor(parent(X),parent(X), X)
end

⊗(::typeof(-), X::Object) = RightTensorProductFunctor(X)

is_additive(T::RightTensorProductFunctor) = is_additive(domain(T))

#-------------------------------------------------------------------------------
#   Restriction and Induction
#-------------------------------------------------------------------------------

struct GRepRestriction <: AbstractFunctor
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

struct GRepInduction <: AbstractFunctor
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
#     @assert is_monoidal(C) "Category has to be monoidal"
#     obj_map = X -> dual(X)
#     mor_map = f -> dual(f)
#     return DualFunctor(C,C,obj_map,mor_map)
# end

#= @memoize Dict =# function dual_monoidal_structure(X::Object, Y::Object)
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


#=----------------------------------------------------------
    Induction Functor 
----------------------------------------------------------=#

struct InductionFunctor <: AbstractFunctor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

function Induction(C::Category)
    @assert is_multifusion(C)
    obj_map = x -> induction(x)
    mor_map = induction_mor_map
    return InductionFunctor(C,Center(C),obj_map,mor_map)
end

function induction_mor_map(f::Morphism)
    S = simples(parent(domain(f)))

    return direct_sum([id(s)⊗f⊗id(dual(s)) for s ∈ S])
end