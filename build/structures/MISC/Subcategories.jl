abstract type AbstractSubcategory <: Category end

struct RingSubcategory <: AbstractSubcategory
    category::Category
    simples::Vector{<:CategoryObject}
    projector::CategoryObject
end

struct SubcategoryCategoryObject <: CategoryObject
    parent::AbstractSubcategory
    object::CategoryObject
end

struct SubcategoryMorphism <: CategoryMorphism
    domain::SubcategoryCategoryObject
    codomain::SubcategoryCategoryObject
    m::CategoryMorphism
end


object(X::SubcategoryCategoryObject) = X.object
morphism(f::SubcategoryMorphism) = f.m
base_ring(C::RingSubcategory) = base_ring(C.category)
#=-------------------------------------------------
    Constructors 
-------------------------------------------------=#


function RingSubcategory(C::Category,i::Int...)
    @assert is_multitensor(C)
    𝟙ᵢ = direct_sum([decompose(one(C))[iₖ][1] for iₖ ∈ i]...)
    projection = [𝟙ᵢ⊗S⊗𝟙ᵢ for S ∈ simples(C)]
    filter!(e -> e != zero(C), projection)
    return RingSubcategory(C,projection,𝟙ᵢ)
end

RingSubcategory(C::Category,i::Vector{Int}) = RingSubcategory(C,i...)
#=-------------------------------------------------
    Functionality 
-------------------------------------------------=#
function direct_sum(X::SubcategoryCategoryObject, Y::SubcategoryCategoryObject)
    @assert parent(X) == parent(Y)
    obj = direct_sum(object(X), object(Y))
    return SubcategoryCategoryObject(parent(X), obj)
end

function direct_sum(f::SubcategoryMorphism, g::SubcategoryMorphism)
    @assert parent(f) == parent(g)
    mor = direct_sum(morphism(f), morphism(g))
    return SubcategoryMorphism(domain(f)⊕domain(g), codomain(f)⊕codomain(g), mor)
end

function direct_sum(X::SubcategoryCategoryObject, Y::SubcategoryCategoryObject)
    @assert parent(X) == parent(Y)
    obj,ix,px = direct_sum(object(X), object(Y))
    sub_obj = SubcategoryCategoryObject(parent(X), obj)
    sub_ix = [SubcategoryMorphism(x,sub_obj,i) for (i,x) ∈ zip(ix,[X,Y])]
    sub_px = [SubcategoryMorphism(sub_obj,y,p) for (p,y) ∈ zip(px,[X,Y])]
    return sub_obj, sub_ix, sub_px
end

function tensor_product(X::SubcategoryCategoryObject, Y::SubcategoryCategoryObject)
    @assert parent(X) == parent(Y)
    obj = tensor_product(object(X),object(Y))
    return SubcategoryCategoryObject(parent(X), obj)
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

function dual(X::SubcategoryCategoryObject)
    return SubcategoryCategoryObject(parent(X), dual(object(X)))
end

function ev(X::SubcategoryCategoryObject)
    dom = dual(X)⊗X
    cod = one(parent(X))
    proj = basis(Hom(one(parent(X).category), parent(X).projector))[1]
    return SubcategoryMorphism(dom,cod, proj ∘ ev(object(X)))
end

function coev(X::SubcategoryCategoryObject)
    incl = basis(Hom(parent(X).projector, one(parent(X).category)))[1]
    return SubcategoryMorphism(X⊗dual(X),one(parent(X)), coev(object(X)) ∘ incl)
end
    
function spherical(X::SubcategoryCategoryObject)
    return SubcategoryMorphism(X,dual(dual(X)), spherical(object(X)))
end

function id(X::SubcategoryCategoryObject) 
    return SubcategoryMorphism(X,X, id(object(X)))
end

function zero_morphism(X::SubcategoryCategoryObject, Y::SubcategoryCategoryObject)
    return SubcategoryMorphism(X,Y, zero_morphism(object(X),object(Y)))
end

function Hom(X::SubcategoryCategoryObject, Y::SubcategoryCategoryObject)
    sub_basis = [SubcategoryMorphism(X,Y,f) for f ∈ Hom(object(X),object(Y))]
    return CategoryHomSpace(X,Y,sub_basis,VectorSpaces(base_ring(X)))
end

function is_isomorphic(X::SubcategoryCategoryObject, Y::SubcategoryCategoryObject)
    b, iso = is_isomorphic(object(X),object(Y))
    if !b 
        return false, nothing
    end
    return true, SubcategoryMorphism(X,Y, iso)
end

function kernel(f::SubcategoryMorphism)
    @assert is_abelian(parent(f))
    k,i = kernel(morphism(f))
    sub_k = SubcategoryCategoryObject(parent(f), k)
    return sub_k, SubcategoryMorphism(sub_k, domain(f), i)
end

function cokernel(f::SubcategoryMorphism)
    @assert is_abelian(parent(f))
    c,i = cokernel(morphism(f))
    sub_c = SubcategoryCategoryObject(parent(f), c)
    return sub_c, SubcategoryMorphism(codomain(f),sub_c, i)
end

left_inverse(f::SubcategoryMorphism) = SubcategoryMorphism(codomain(f),domain(f), left_inverse(morphism(f)))
right_inverse(f::SubcategoryMorphism) = SubcategoryMorphism(codomain(f),domain(f), right_inverse(morphism(f)))

is_simple(X::SubcategoryCategoryObject) = is_simple(object(X))

matrix(f::SubcategoryMorphism) = matrix(morphism(f))

is_semisimple(C::AbstractSubcategory) = is_semisimple(C.category)
is_multifusion(C::AbstractSubcategory) = is_multifusion(C.category)

*(x, f::SubcategoryMorphism) = SubcategoryMorphism(domain(f),codomain(f), x*morphism(f))
+(f::SubcategoryMorphism, g::SubcategoryMorphism) = SubcategoryMorphism(domain(f),codomain(f), morphism(f) + morphism(g))
inv(f::SubcategoryMorphism) = SubcategoryMorphism(codomain(f),domain(f), inv(morphism(f)))

function express_in_basis(f::SubcategoryMorphism, B::Vector{SubcategoryMorphism})
    express_in_basis(morphism(f), morphism.(B))
end
#=-------------------------------------------------
    Functionality for RingSubcategory 
-------------------------------------------------=#

one(C::RingSubcategory) = SubcategoryCategoryObject(C,C.projector)
zero(C::AbstractSubcategory) = SubcategoryCategoryObject(C,zero(C.category))
simples(C::RingSubcategory) = [SubcategoryCategoryObject(C,s) for s in C.simples]

function associator(X::SubcategoryCategoryObject, Y::SubcategoryCategoryObject, Z::SubcategoryCategoryObject)
    dom = (X ⊗ Y) ⊗ Z
    cod = X ⊗ (Y ⊗ Z)
    return SubcategoryMorphism(dom,cod, associator(object(X), object(Y), object(Z)))
end

is_fusion(C::RingSubcategory) = is_multifusion(C.category) && length(decompose(C.projector)) == 1
#=-------------------------------------------------
    Pretty Printing 
-------------------------------------------------=#

function show(io::IO, X::SubcategoryCategoryObject)
    print(io, """(Subcategory) $(object(X))""")
end

function show(io::IO, f::SubcategoryMorphism)
    print(io, """(Subcategory) $(morphism(f))""")
end

function show(io::IO, C::RingSubcategory)
   # i = findfirst(e -> e == C.projector, [c for (c,k) ∈ decompose(one(C.category))])
    print(io, """Ring subcategory of $(C.category)""")
end