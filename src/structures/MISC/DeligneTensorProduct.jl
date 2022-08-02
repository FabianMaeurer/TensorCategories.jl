struct DeligneProduct <: Category
    C::Category 
    D::Category 
    base_ring::Field
end

abstract type DeligneProdObject <: Object end

struct SimpleDeligneProdObject <: DeligneProdObject
    X::Object 
    Y::Object 
    parent::DeligneProduct 
end

struct CompositeDeligneProdObject <: DeligneProdObject
    summands::Vector{SimpleDeligneProdObject}
    parent::DeligneProduct
end

abstract type DeligneProdMorphism <: Morphism end

struct SimpleDeligneProdMorphism <: DeligneProdMorphism 
    domain::SimpleDeligneProdObject 
    codomain::SimpleDeligneProdObject 
    m::Tuple{Morphism,Morphism}
end

struct CompositeDeligneProdMorphism <: DeligneProdMorphism
    domain::CompositeDeligneProdMorphism
    codomain::CompositeDeligneProdMorphism
    m::Matrix{simpleDeligneProdMorphism}
end

#=------------------------------------------------
    Constructors
------------------------------------------------=#
function DeligneProduct(C::Category, D::Category)
    @assert base_ring(C) == base_ring(D) "Missmatching base fields"
    @assert issemisimple(C) "Categories have to be semisimple"
    @assert issemisimple(D) "Categories have to be semisimple"
    DeligneProduct(C,D,base_ring(C))
end    

⊠(C::Category, D::Category) = DeligneProduct(C,D)

function DeligneProdMorphism(f::Morphism, g::Morphism)
    dom = domain(f) ⊠ domain(g)
    cod = codomain(f) ⊠ codomain(g)
    SimpleDeligneProdMorphism(dom, cod, (f,g))
end

⊠(f::Morphism, g::Morphism) = DeligneProdMorphism(f,g)

function DeligneProdObject(X::Object, Y::Object)
    SimpleDeligneProdObject(X,Y,parent(X) ⊠ parent(Y))
end

⊠(X::Object, Y::Object) = DeligneProdObject(X,Y)

#=------------------------------------------------
    simple objects
------------------------------------------------=#

function simples(C::DeligneProduct) 
    C_simples = simples(C.C)
    D_simples = simples(C.D)
    [DeligneProdObject(x,y) for x ∈ C_simples, y ∈ D_simples][:]
end

#=------------------------------------------------
    Implied Structures
------------------------------------------------=#

one(C::DeligneProduct) = one(C.C) ⊠ one(C.D)

zero(C::DeligneProduct) = zero(C.C) ⊠ zero(C.D)

dsum(X::DeligneProdObject, Y::DeligneProdObject) = ()

#=------------------------------------------------
    Checks
------------------------------------------------=#

issemisimple(C::DeligneProduct) = issemisimple(C.C) && issemisimple(C.D)
isabelian(C::DeligneProduct) = isabelian(C.C) && isabelian(C.D)
islinear(C::DeligneProduct) = islinear(C.C) && islinear(C.D)
ismonoidal(C::DeligneProduct) = ismonoidal(C.C) && ismonoidal(C.D)
isfinite(C::DeligneProduct) = isfinte(C.C) && isfinite(C.D)
istensor(C::DeligneProduct) = istensor(C.C) && istensor(C.D)
ismultitensor(C::DeligneProduct) = ismultitensor(C.C) && ismultitensor(C.D)
isfusion(C::DeligneProduct) = isfusion(C.C) && isfusion(C.D)
ismultifusion(C::DeligneProduct) = ismultifusion(C.C) && ismultifusion(C.D)

#=------------------------------------------------
    Homspaces
------------------------------------------------=#

function Hom(X::SimpleDeligneProdObject, Y::SimpleDeligneProdObject)
    hom_X1_Y1 = basis(Hom(X.X,Y.X))
    hom_X2_Y2 = basis(Hom(X.Y,Y.Y))
    
    hom_basis = [f⊠g for f ∈ hom_X1_Y1, g ∈ hom_X2_Y2][:]

    return HomSpace(X,Y,hom_basis, VectorSpaces(base_ring(X)))
end