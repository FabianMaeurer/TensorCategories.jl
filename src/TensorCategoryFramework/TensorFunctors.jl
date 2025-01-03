
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
    identity 
----------------------------------------------------------=#

struct IdentityFunctor <: AbstractFunctor 
    C::Category
end

domain(F::IdentityFunctor) = F.C
codomain(F::IdentityFunctor) = F.C 

(::IdentityFunctor)(X::Object) = X
(::IdentityFunctor)(f::Morphism) = f

id(C::Category) = IdentityFunctor(C)

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
    return InductionFunctor(C,center(C),obj_map,mor_map)
end

function induction_mor_map(f::Morphism)
    S = simples(parent(domain(f)))

    return direct_sum([id(s)⊗f⊗id(dual(s)) for s ∈ S])
end


#=----------------------------------------------------------
    Inner Endofunctors
----------------------------------------------------------=#

abstract type MonoidalFunctor <: AbstractFunctor end 

struct InnerAutoequivalence <: MonoidalFunctor 
    domain::Category
    object::Object 
    dual_object::Object

    monoidal_structure::Dict
end 



function inner_autoequivalence(C::Category, X::Object) 
    InnerAutoequivalence(C,X,dual(X), Dict())
end

inner_autoequivalence(X::Object) = new(parent(X), X)

function inner_autoequivalence(X::Object) 
    @req is_invertible(X) "Object has to be invertible"
    InnerAutoequivalence(X)
end

(F::InnerAutoequivalence)(X::Object) = (F.object ⊗ X) ⊗ F.dual_object 
(F::InnerAutoequivalence)(f::Morphism) = (id(F.object) ⊗ f) ⊗ id(F.dual_object)

function monoidal_structure(F::InnerAutoequivalence, S::Object, T::Object) 
    get(F.monoidal_structure, (S,T)) do 
        X = F.object 
        dX = F.dual_object

        F.monoidal_structure[(S,T)] = compose(
            associator(X⊗S, dX, (X⊗T)⊗dX),
            id(X⊗S) ⊗ compose(
                (id(dX) ⊗ associator(X,T,dX)),
                inv_associator(dX,X,T⊗dX),
                ev(X) ⊗ id(T⊗dX)
            ),
            inv_associator(X⊗S,T,dX),
            associator(X,S,T) ⊗ id(dX)
        )
    end
end
    



    
