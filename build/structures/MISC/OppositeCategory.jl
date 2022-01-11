struct OppositeCategory <: Category
    C::Category
end

struct OppositeMorphism <: Morphism
    m::Morphism
end

struct OppositeObject <: Object
    X::Object
end

domain(f::OppositeMorphism) = codomain(f.m)
codomain(f::OppositeMorphism) = domain(f.m)

base_ring(C::OppositeCategory) = base_ring(C.C)
base_ring(X::OppositeObject) = base_ring(X.X)
parent(X::OppositeObject) = OppositeCategory(parent(X.X))

compose(f::OppositeMorphism, g::OppositeMorphism) = OppositeMorphism(compose(f.m,g.m))

function product(X::OppositeObject, Y::OppositeObject)
    Z,px = product(X.X,Y.X)
    return OppositeObject(Z), OppositeMorphism.(px)
end

function coproduct(X::OppositeObject, Y::OppositeObject)
    Z,ix = coproduct(X.X,Y.X)
    return OppositeObject(Z), OppositeMorphism.(ix)
end

function dsum(X::OppositeObject, Y::OppositeObject)
    Z,ix,px = coproduct(X.X,Y.X)
    return OppositeObject(Z), OppositeMorphism.(ix), OppositeMorphism.(px)
end

tensor_product(X::OppositeObject, Y::OppositeObject) = OppositeObject(tensor_product(X.X,Y.X))

#-------------------------------------------------------------------------------
#   Inversion
#-------------------------------------------------------------------------------

OppositeCategory(C::OppositeCategory) = C.C
OppositeObject(X::OppositeObject) = X.X
OppositeMorphism(f::OppositeMorphism) = f.m
