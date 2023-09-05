mutable struct TensorPowerCategory <: Category
    generator::Object
    simples::Vector
    complete::Bool
    max_exponent::Int

    function TensorPowerCategory(X::Object) 
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

struct TensorPowerObject <: Object 
    parent::TensorPowerCategory
    object::Object
end

struct TensorPowerMorphism <: Morphism
    domain::TensorPowerObject
    codomain::TensorPowerObject
    morphism::Morphism
end

object(X::TensorPowerObject) = X.object
morphism(f::TensorPowerMorphism) = f.m
category(C::TensorPowerCategory) = parent(C.generator)
Morphism(X::TensorPowerObject, Y::TensorPowerObject, f::TensorPowerMorphism) = TensorPowerMorphism(X,Y,f)

base_ring(C::TensorPowerCategory) = base_ring(category(C))

""" 

    tensor_power(X::Object, k::Int) -> Object

Return the ``k``-th tensor power ``X^{\\otimes k}``.
"""
function tensor_power(X::Object, k::Int)
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

⊗(X::Object,k::Int) = tensor_power(X,k)

function direct_sum(X::TensorPowerObject, Y::TensorPowerObject)
    direct_sum(object(X), object(Y))
end

function tensor_product(X::TensorPowerObject, Y::TensorPowerObject)
    tensor_product(object(X), object(Y))
end

function direct_sum(X::TensorPowerMorphism, Y::TensorPowerMorphism)
    direct_sum(morphism(X), morphism(Y))
end

function tensor_product(X::TensorPowerMorphism, Y::TensorPowerMorphism)
    tensor_product(morphism(X), morphism(Y))
end


one(C::TensorPowerCategory) = TensorPowerObject(C, one(category(C)))

#=----------------------------------------------------------
    Simples/Indecompodables 
----------------------------------------------------------=#

function indecomposable_subobjects(X::TensorPowerObject)
    subs = indecomposable_subobjects(object(X))
    return [TensorPowerObject(parent(X), s) for s ∈ subs]
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
            simpls = [TensorPowerObject(C,s) for s ∈ simpls]
            C.simples = simpls
            C.complete = true
            C.max_exponent = j-1
            return simpls
        end
        n1 = length(simpls)
        Y = Y ⊗ X
        j = j+1
    end
    simpls = [TensorPowerObject(C,s) for s ∈ simpls]
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
            simpls = [TensorPowerObject(C,s) for s ∈ simpls]
            C.simples = simpls
            C.complete = true
            C.max_exponent = j-1
            return simpls
        end
        n1 = length(simpls)
        Y = Y ⊗ X
        j = j+1
    end
    simpls = [TensorPowerObject(C,s) for s ∈ simpls]
    C.simpls = simpls
    C.complete = false
    C.max_exponent = k
    return simpls
end

function show(io::IO, C::TensorPowerCategory)
    print(io, "Tensor power category with genrator $(C.generator)")
end

function show(io::IO, X::TensorPowerObject)
    print(io, "Tensor power category object: $(object(X))")
end

function show(io::IO, f::TensorPowerMorphism)
    print(io, "Tensor power category morphism:  $(morphism(f))")
end


#=----------------------------------------------------------
    Functionality 
----------------------------------------------------------=#

compose(f::TensorPowerMorphism, g::TensorPowerMorphism) = TensorPowerMorphism(domain(f), codomain(g), compose(morphism(f),morphism(g)))


*(k,f::TensorPowerMorphism) = TensorPowerMorphism(domain(f),codomain(f), k*morphism(f))

+(f::TensorPowerMorphism, g::TensorPowerMorphism) = TensorPowerMorphism(domain(f),codomain(f), morphism(f) + morphism(g))

matrix(f::TensorPowerMorphism) = matrix(morphism(f))

function Hom(X::TensorPowerObject, Y::TensorPowerObject)
    H = Hom(object(X), object(Y))
    B = [TensorPowerMorphism(X,Y,f) for f ∈ basis(H)]
    HomSpace(X,Y,B,VectorSpaces(base_ring(X)))
end