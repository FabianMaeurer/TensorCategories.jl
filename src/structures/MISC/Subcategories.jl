abstract type AbstractSubcategory <: Category end

struct RingSubcategory <: AbstractSubcategory
    category::Category
    simples::Vector{<:Object}
    projector::Object
end

struct SubcategoryObject <: Object
    parent::AbstractSubcategory
    object::Object
end

struct SubcategoryMorphism <: Morphism
    domain::SubcategoryObject
    codomain::SubcategoryObject
    m::Morphism
end


object(X::SubcategoryObject) = X.object
morphism(f::SubcategoryMorphism) = f.m

#=-------------------------------------------------
    Constructors 
-------------------------------------------------=#

function RingSubcategory(C::Category,i::Int)
    @assert ismultitensor(C)
    𝟙ᵢ = decompose(one(C))[i][1] 
    projection = [𝟙ᵢ⊗S⊗𝟙ᵢ for S ∈ simples(C)]
    filter!(e -> e != zero(C), projection)
    return RingSubcategory(C,projection,𝟙ᵢ)
end

#=-------------------------------------------------
    Functionality 
-------------------------------------------------=#
function dsum(X::SubcategoryObject, Y::SubcategoryObject)
    @assert parent(X) == parent(Y)
    obj = dsum(object(X), object(Y))
    return SubcategoryObject(parent(X), obj)
end

function dsum(f::SubcategoryMorphism, g::SubcategoryMorphism)
    @assert parent(f) == parent(g)
    mor = dsum(morphism(f), morphism(g))
    return SubcategoryMorphism(domain(f)⊕domain(g), codomain(f)⊕codomain(g), mor)
end

function dsum_with_morphisms(X::SubcategoryObject, Y::SubcategoryObject)
    @assert parent(X) == parent(Y)
    obj,ix,px = dsum_with_morphisms(object(X), object(Y))
    sub_obj = SubcategoryObject(parent(X), obj)
    sub_ix = [SubcategoryMorphism(x,sub_obj,i) for (i,x) ∈ zip(ix,[X,Y])]
    sub_px = [SubcategoryMorphism(sub_obj,y,p) for (p,y) ∈ zip(px,[X,Y])]
    return sub_obj, sub_ix, sub_px
end

function tensor_product(X::SubcategoryObject, Y::SubcategoryObject)
    @assert parent(X) == parent(Y)
    obj = tensor_product(object(X),object(Y))
    return SubcategoryObject(parent(X), obj)
end

function tensor_product(f::SubcategoryMorphism, g::SubcategoryMorphism)
    @assert parent(f) == parent(g)
    mor = tensor_product(morphism(f),morphism(g))
    return SubcategoryMorphism(domain(f)⊗domain(g), codomain(f)⊗codomain(g), mor)
end

function compose(f::SubcategoryMorphism, g::SubcategoryMorphism)
    @assert parent(f) == parent(g)
    return SubcategoryMorphism(domain(f),codomain(g), compose(morphism(f),morphism(g)))
end

function dual(X::SubcategoryObject)
    return SubcategoryObject(parent(X), dual(object(X)))
end

function ev(X::SubcategoryObject)
    dom = dual(X)⊗X
    cod = one(parent(X))
    return SubcategoryMorphism(dom,cod, ev(object(X)))
end

function coev(X::SubcategoryObject)
    return SubcategoryMorphism(X⊗dual(X),one(parent(X)), coev(object(X)))
end
    
function spherical(X::SubcategoryObject)
    return SubcategoryMorphism(X,dual(dual(X)), spherical(object(X)))
end

issemisimple(C::AbstractSubcategory) = issemisimple(C.category)
ismultifusion(C::AbstractSubcategory) = ismultifusion(C.category)
#=-------------------------------------------------
    Functionality for RingSubcategory 
-------------------------------------------------=#

one(C::RingSubcategory) = SubcategoryObject(C,C.projector)
simples(C::RingSubcategory) = [SubcategoryObject(C,s) for s in C.simples]

function associator(X::SubcategoryObject, Y::SubcategoryObject, Z::SubcategoryObject)
    dom = (X ⊗ Y) ⊗ Z
    cod = X ⊗ (Y ⊗ Z)
    return SubcategoryMorphism(dom,cod, associator(object(X), obejct(Y), object(Z)))
end