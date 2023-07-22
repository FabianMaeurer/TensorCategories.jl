mutable struct TensorPowerCategory <: Category
    generator::CategoryObject
    simples::Vector
    complete::Bool
    max_exponent::Int

    function TensorPowerCategory(X::CategoryObject) 
        C = new()
        C.generator = X
        C.simples = typeof(X)[]

        C.complete = X == zero(parent(X)) ? true : false
        C.max_exponent = 0
        return C
    end
end

function ==(C::TensorPowerCategory, D::TensorPowerCategory)
    C.generator == D.generator
end

struct TensorPowerCategoryObject <: CategoryObject 
    parent::TensorPowerCategory
    object::CategoryObject
end

struct TensorPowerCategoryMorphism <: CategoryMorphism
    domain::TensorPowerCategoryObject
    codomain::TensorPowerCategoryObject
    morphism::CategoryMorphism
end

object(X::TensorPowerCategoryObject) = X.object
morphism(f::TensorPowerCategoryMorphism) = f.m
category(C::TensorPowerCategory) = parent(C.generator)
Morphism(X::TensorPowerCategoryObject, Y::TensorPowerCategoryObject, f::TensorPowerCategoryMorphism) = TensorPowerCategoryMorphism(X,Y,f)

base_ring(C::TensorPowerCategory) = base_ring(category(C))

""" 

    tensor_power(X::CategoryObject, k::Int) -> CategoryObject

Return the ``k``-th tensor power ``X^{\\otimes k}``.
"""
function tensor_power(X::CategoryObject, k::Int)
    if k < 0 
        error("Negative exponent")
    elseif k == 0
        return one(parent(X))
    elseif k == 1
        return X
    end

    if isodd(k)
        return X ⊗ tensor_power(X,k-1)
    else 
        Y = tensor_power(X, div(k,2))
        return Y ⊗ Y
    end
    return Y
end

⊗(X::CategoryObject,k::Int) = tensor_power(X,k)

function direct_sum(X::TensorPowerCategoryObject, Y::TensorPowerCategoryObject)
    direct_sum(object(X), object(Y))
end

function tensor_product(X::TensorPowerCategoryObject, Y::TensorPowerCategoryObject)
    tensor_product(object(X), object(Y))
end

function direct_sum(X::TensorPowerCategoryMorphism, Y::TensorPowerCategoryMorphism)
    direct_sum(morphism(X), morphism(Y))
end

function tensor_product(X::TensorPowerCategoryMorphism, Y::TensorPowerCategoryMorphism)
    tensor_product(morphism(X), morphism(Y))
end


one(C::TensorPowerCategory) = TensorPowerCategoryObject(C, one(category(C)))

#=----------------------------------------------------------
    Simples/Indecompodables 
----------------------------------------------------------=#

function indecomposable_subobjects(X::TensorPowerCategoryObject)
    subs = indecomposable_subobjects(object(X))
    return [TensorPowerCategoryObject(parent(X), s) for s ∈ subs]
end

function simples(C::TensorPowerCategory, k = Inf)
    simpls = object_type(category(C))[]
    n1 = 0
    j = 0
    X = C.generator
    Y = one(category(C))
    while j ≤ k+1
        simpls = unique_simples([simpls; simple_subobjects(Y)])
        if length(simpls) == n1
            simpls = [TensorPowerCategoryObject(C,s) for s ∈ simpls]
            C.simples = simpls
            C.complete = true
            C.max_exponent = j-1
            return simpls
        end
        n1 = length(simpls)
        Y = Y ⊗ X
        j = j+1
    end
    simpls = [TensorPowerCategoryObject(C,s) for s ∈ simpls]
    C.simpls = simpls
    C.complete = false
    C.max_exponent = k
    return simpls
end


function indecomposables(C::TensorPowerCategory, k = Inf)
    X = C.generator
    Y = one(category(C))
    indecs_in_X = indecomposable_subobjects(X)
    simpls = indecs_in_X
    n1 = 0
    j = 2
    
    while j ≤ k
        for V ∈ indecs_in_X, W ∈ simpls
            simpls = unique_indecomposables([simpls; indecomposable_subobjects(W ⊗ V)])
        end
        if length(simpls) == n1
            simpls = [TensorPowerCategoryObject(C,s) for s ∈ simpls]
            C.simples = simpls
            C.complete = true
            C.max_exponent = j-1
            return simpls
        end
        n1 = length(simpls)
        Y = Y ⊗ X
        j = j+1
    end
    simpls = [TensorPowerCategoryObject(C,s) for s ∈ simpls]
    C.simpls = simpls
    C.complete = false
    C.max_exponent = k
    return simpls
end

function show(io::IO, C::TensorPowerCategory)
    print(io, "Tensor power category with genrator $(C.generator)")
end

function show(io::IO, X::TensorPowerCategoryObject)
    print(io, "Tensor power category object: $(object(X))")
end

function show(io::IO, f::TensorPowerCategoryMorphism)
    print(io, "Tensor power category morphism:  $(morphism(f))")
end


#=----------------------------------------------------------
    Functionality 
----------------------------------------------------------=#

compose(f::TensorPowerCategoryMorphism, g::TensorPowerCategoryMorphism) = TensorPowerCategoryMorphism(domain(f), codomain(g), compose(morphism(f),morphism(g)))


*(k,f::TensorPowerCategoryMorphism) = TensorPowerCategoryMorphism(domain(f),codomain(f), k*morphism(f))

+(f::TensorPowerCategoryMorphism, g::TensorPowerCategoryMorphism) = TensorPowerCategoryMorphism(domain(f),codomain(f), morphism(f) + morphism(g))

matrix(f::TensorPowerCategoryMorphism) = matrix(morphism(f))

function Hom(X::TensorPowerCategoryObject, Y::TensorPowerCategoryObject)
    H = Hom(object(X), object(Y))
    B = [TensorPowerCategoryMorphism(X,Y,f) for f ∈ basis(H)]
    CategoryHomSpace(X,Y,B,VectorSpaces(base_ring(X)))
end