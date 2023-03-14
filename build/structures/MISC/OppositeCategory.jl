struct OppositeCategory <: Category
    C::Category
end

struct OppositeCategoryObject <: CategoryObject
    parent::OppositeCategory
    X::CategoryObject
end


struct OppositeMorphism <: CategoryMorphism 
    domain::OppositeCategoryObject
    codomain::OppositeCategoryObject
    m::CategoryMorphism
end

op(C::Category) = OppositeCategory(C)
op(X::CategoryObject) = OppositeCategoryObject(op(parent(X)),X)
op(f::CategoryMorphism) = OppositeMorphism(op(codomain(f)),op(domain(f)), f)

base_ring(C::OppositeCategory) = base_ring(C.C)
base_ring(X::OppositeCategoryObject) = base_ring(X.X)
parent(X::OppositeCategoryObject) = OppositeCategory(parent(X.X))

compose(f::OppositeMorphism, g::OppositeMorphism) = OppositeMorphism(compose(g.m,f.m))

function product(X::OppositeCategoryObject, Y::OppositeCategoryObject)
    Z,px = product(X.X,Y.X)
    return OppositeCategoryObject(Z), OppositeMorphism.(px)
end

function coproduct(X::OppositeCategoryObject, Y::OppositeCategoryObject)
    Z,ix = coproduct(X.X,Y.X)
    return OppositeCategoryObject(Z), OppositeMorphism.(ix)
end

function direct_sum(X::OppositeCategoryObject, Y::OppositeCategoryObject)
    Z,ix,px = coproduct(X.X,Y.X)
    return OppositeCategoryObject(Z), OppositeMorphism.(ix), OppositeMorphism.(px)
end

tensor_product(X::OppositeCategoryObject, Y::OppositeCategoryObject) = OppositeCategoryObject(tensor_product(X.X,Y.X))

one(C::OppositeCategory) = OppositeCategoryObject(one(C.C))

simples(C::OppositeCategory) = op.(simples(C.C))

zero(C::OppositeCategory) = op(zero(C.C))

#-------------------------------------------------------------------------------
#   Inversion
#-------------------------------------------------------------------------------

OppositeCategory(C::OppositeCategory) = C.C
OppositeCategoryObject(X::OppositeCategoryObject) = X.X
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

function Hom(X::OppositeCategoryObject, Y::OppositeCategoryObject)
    opposite_basis = op.(basis(Hom(X.X,Y.X)))
    CategoryHomSpace(X,Y,opposite_basis, VectorSpaces(base_ring(X)))
end