
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
    domain = product_category(C,C)
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

compose(F::IdentityFunctor...) = IdentityFunctor(F[1].C)

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

#=  =# function dual_monoidal_structure(X::Object, Y::Object)
    #(ev(X⊗Y)⊗id(dual(Y)⊗dual(X))) ∘ inv_associator(dual(X⊗Y),X⊗Y,dual(Y)⊗dual(X)) ∘ (id(dual(X⊗Y))⊗product_coev(X,Y))
    dXY = dual(X⊗Y)
    dX = dual(X)
    dY = dual(Y)
    XY = (X ⊗ Y)
    compose(
        id(dXY) ⊗ compose(
            coev(X),
            (id(X) ⊗ coev(Y)) ⊗ id(dX),
            inv_associator(X, Y, dual(Y)) ⊗ id(dX),
        ),
        inv_associator(dXY, XY ⊗ dY, dual(X)),
        inv_associator(dXY, XY, dY) ⊗ id(dX),
        ev(dXY) ⊗ id(dY) ⊗ id(dX)
    )   
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

abstract type AbstractMonoidalFunctor <: AbstractFunctor end 

struct InnerAutoequivalence <: AbstractMonoidalFunctor 
    domain::Category
    object::Object 
    dual_object::Object
    indecomposables::Vector{<:Object}
    monoidal_structure::Dict
end 



function inner_autoequivalence(C::Category, X::Object) 
    dX = dual(X)

    indecs = indecomposables(C) 
    n = length(indecs) 

    mon_structure = Dict()
    for i ∈ 1:n, j ∈ 1:n
        S,T = indecs[[i,j]]
        mon_structure[(i,j)] = compose(
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
    InnerAutoequivalence(C,X, dX, indecs, mon_structure)
end

indecomposables(F::InnerAutoequivalence) = indecomposables(domain(F))


function inner_autoequivalence(X::Object) 
    @req is_invertible(X) "Object has to be invertible"
    inner_autoequivalence(parent(X),X)
end

codomain(F::InnerAutoequivalence) = domain(F)

function functor(F::InnerAutoequivalence) 
    functor(
        domain(F),
        codomain(F),
        X -> F(X),
        f -> F(f)
    )
end

(F::InnerAutoequivalence)(X::Object) = (F.object ⊗ X) ⊗  F.dual_object
(F::InnerAutoequivalence)(f::Morphism) = (id(F.object) ⊗ f) ⊗ id( F.dual_object)

# function monoidal_structure(F::InnerAutoequivalence, S::Object, T::Object)  
#     X = F.object 
#     dX = F.dual_object

#     compose(
#         associator(X⊗S, dX, (X⊗T)⊗dX),
#         id(X⊗S) ⊗ compose(
#             (id(dX) ⊗ associator(X,T,dX)),
#             inv_associator(dX,X,T⊗dX),
#             ev(X) ⊗ id(T⊗dX)
#         ),
#         inv_associator(X⊗S,T,dX),
#         associator(X,S,T) ⊗ id(dX)
#     )
# end

function compose(F::InnerAutoequivalence, G::InnerAutoequivalence)
    X,dX = F.object, F.dual_object
    Y,dY = G.object, G.dual_object
    FG = Functor(
        domain(F),
        codomain(G),
        V -> (Y ⊗ ((X ⊗ V) ⊗ dX)) ⊗ dY,
        f -> (id(Y) ⊗ ((id(X) ⊗ f) ⊗ id(dX))) ⊗ id(dY)
    )

    S = indecomposables(domain(F))

    MonoidalFunctor(
        FG,
        S,
        Dict((i,j) => compose(
                monoidal_structure(G, F(S[i]), F(S[j])),
                G(monoidal_structure(F, S[i], S[j]))
            ) for i ∈ 1:length(S), j ∈ 1:length(S)
        )   
    )
end

function show(io::IO, F::InnerAutoequivalence)
    print(io, "Inner autoequivalence defined by $(F.object)")
end

#=----------------------------------------------------------
    TensorFunctors 
----------------------------------------------------------=#

mutable struct MonoidalFunctor <: AbstractMonoidalFunctor
    F::AbstractFunctor 
    indecomposables::Vector{<:Object}
    monoidal_structure::Dict
end

(F::MonoidalFunctor)(X::Object) = F.F(X)
(F::MonoidalFunctor)(f::Morphism) = F.F(f)

functor(F::MonoidalFunctor) = F.F

domain(F::MonoidalFunctor) = domain(F.F)
codomain(F::MonoidalFunctor) = codomain(F.F)
indecomposables(F) = F.indecomposables

function monoidal_functor(F::AbstractFunctor, S::Vector{<:Object}, nats::Dict)
    MonoidalFunctor(F,S,nats)
end

function compose(F::AbstractMonoidalFunctor, G::AbstractMonoidalFunctor)
    S = indecomposables(F)
    MonoidalFunctor(
        compose(functor(F), functor(G)),
        S,
        Dict((i,j) => compose(
                monoidal_structure(G, F(S[i]), F(S[j])),
                G(monoidal_structure(F, S[i], S[j]))
            ) for i ∈ 1:length(S), j ∈ 1:length(S)
        )            
    )
end


function show(io::IO, F::MonoidalFunctor)
    print(io, """Tensor functor with 
    domain:   $(domain(F))
    codomain: $(codomain(F))""")
end

function monoidal_structure(F::AbstractMonoidalFunctor, X::Object, Y::Object) 
    S = indecomposables(F)

    if is_indecomposable(X) && is_indecomposable(Y) 
        local iso_X
        local iso_Y
        local i,j

        if X ∈ S 
            iso_X = id(X)
            i = findfirst(==(X), S)
        else
            i = findfirst(x -> is_isomorphic(x,X)[1], S)
            _,iso_X = is_isomorphic(X, S[i])
        end

        if Y ∈ S
            iso_Y = id(Y)
            j = findfirst(==(Y), S)
        else
            j = findfirst(x -> is_isomorphic(x,Y)[1], S)
            _,iso_Y = is_isomorphic(Y, S[j])
        end

        return compose(
            F(iso_X) ⊗ F(iso_Y),
            F.monoidal_structure[(i,j)],
            F(inv(iso_X) ⊗ inv(iso_Y))
        )
    end

    _,_,ix,px = direct_sum_decomposition(X,S)
    _,_,iy,py = direct_sum_decomposition(Y,S)

    sum([compose(
        F(p)⊗F(q),
        monoidal_structure(F,codomain(p), codomain(q)),
        F(i ⊗ j)
    ) for (p,i) ∈ zip(px,ix), (q,j) ∈ zip(py,iy)])
end