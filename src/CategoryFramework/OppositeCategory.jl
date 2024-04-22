struct OppositeCategory <: Category
    C::Category
end

struct OppositeObject <: Object
    parent::OppositeCategory
    X::Object
end


struct OppositeMorphism <: Morphism 
    domain::OppositeObject
    codomain::OppositeObject
    m::Morphism
end

op(C::Category) = OppositeCategory(C)
op(X::Object) = OppositeObject(op(parent(X)),X)
op(f::Morphism) = OppositeMorphism(op(codomain(f)),op(domain(f)), f)

base_ring(C::OppositeCategory) = base_ring(C.C)
base_ring(X::OppositeObject) = base_ring(X.X)
parent(X::OppositeObject) = OppositeCategory(parent(X.X))

compose(f::OppositeMorphism, g::OppositeMorphism) = OppositeMorphism(compose(g.m,f.m))

function product(X::OppositeObject, Y::OppositeObject)
    Z,px = product(X.X,Y.X)
    return OppositeObject(Z), OppositeMorphism.(px)
end

function coproduct(X::OppositeObject, Y::OppositeObject)
    Z,ix = coproduct(X.X,Y.X)
    return OppositeObject(Z), OppositeMorphism.(ix)
end

function direct_sum(X::OppositeObject, Y::OppositeObject)
    Z,ix,px = coproduct(X.X,Y.X)
    return OppositeObject(Z), OppositeMorphism.(ix), OppositeMorphism.(px)
end

tensor_product(X::OppositeObject, Y::OppositeObject) = OppositeObject(tensor_product(X.X,Y.X))

one(C::OppositeCategory) = OppositeObject(one(C.C))

simples(C::OppositeCategory) = op.(simples(C.C))

zero(C::OppositeCategory) = op(zero(C.C))

#=----------------------------------------------------------
    Constructors 
----------------------------------------------------------=#

@doc raw""" 

    opposite_category(C::Category)

Construct the category ``Cᵒᵖ``.
"""
opposite_category(C::Category) = OppositeCategory(C)

@doc raw""" 

    opposite_object(X::Object)

Regard the object ``X ∈ C`` as an object in ``Cᵒᵖ``.
"""
opposite_object(X::Object) = OppositeObject(opposite_category(parent(X)),X)

@doc raw""" 

    opposite_morphism(f::Morphism)

Regard the morphism ``f ∈ C`` as a morphism in ``Cᵒᵖ``.
"""
opposite_morphism(f::Morphism) = OppositeMorphism(OppositeObject(codomain(f)), OppositeObject(domain(f)), f)

(C::OppositeCategory)(X::Object) = OppositeObject(C,X)
(C::OppositeCategory)(f::Morphism) = OppositeMorphism(C(codomain(f)), C(domain(f)), f)
#-------------------------------------------------------------------------------
#   Inversion
#-------------------------------------------------------------------------------

OppositeCategory(C::OppositeCategory) = C.C
OppositeObject(X::OppositeObject) = X.X
OppositeMorphism(f::OppositeMorphism) = f.m

#=------------------------------------------------
    Checks
------------------------------------------------=#

is_semisimple(C::OppositeCategory) = is_semisimple(C.C)
is_abelian(C::OppositeCategory) = is_abelian(C.C)
is_linear(C::OppositeCategory) = is_linear(C.C)
is_monoidal(C::OppositeCategory) = is_monoidal(C.C)
isfinite(C::OppositeCategory) = isfinte(C.C)
is_tensor(C::OppositeCategory) = is_tensor(C.C)
is_multitensor(C::OppositeCategory) = is_multitensor(C.C)
is_fusion(C::OppositeCategory) = is_fusion(C.C)
is_multifusion(C::OppositeCategory) = is_multifusion(C.C)

#=------------------------------------------------
    Hom-Spaces
------------------------------------------------=#

function Hom(X::OppositeObject, Y::OppositeObject)
    opposite_basis = op.(basis(Hom(X.X,Y.X)))
    HomSpace(X,Y,opposite_basis)
end