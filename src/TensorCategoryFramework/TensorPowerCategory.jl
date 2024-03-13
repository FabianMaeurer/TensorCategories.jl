mutable struct TensorPowerCategory <: Category
    generator::Vector{Object}
    simples::Vector
    complete::Bool
    max_exponent::Int
    #multiplication_table::Dict

    function TensorPowerCategory(X::Object...) 
        C = new()
        C.generator = unique_simples(vcat([[k for (k,_) ∈ decompose(x)] for x ∈ X]...))
        C.simples = typeof(X)[]
        
        C.complete = X == zero(parent(X[1])) ? true : false
        C.max_exponent = 0
        return C
    end

    function TensorPowerCategory(X::Vector{<:Object})
        C = new()
        C.generator = unique_simples(vcat([[k for (k,_) ∈ decompose(x)] for x ∈ X]...))

        C.simples = typeof(X)[]
        
        C.complete = X == zero(parent(X[1])) ? true : false
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

is_additive(::TensorPowerCategory) = true

object(X::TensorPowerObject) = X.object
morphism(f::TensorPowerMorphism) = f.morphism
category(C::TensorPowerCategory) = parent(C.generator[1])
Morphism(X::TensorPowerObject, Y::TensorPowerObject, f::Morphism) = TensorPowerMorphism(X,Y,f)
matrix(f::TensorPowerObject) = matrix(f.morphism)

base_ring(C::TensorPowerCategory) = base_ring(category(C))

dim(X::TensorPowerObject) = dim(object(X))


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

function direct_sum(X::TensorPowerObject...)
    S, incl, proj= direct_sum(object.(X))
    S = TensorPowerObject(parent(X[1]), S)
    incl = [Morphism(x, S, i) for (i,x) ∈ zip(incl, X)]
    proj = [Morphism(S, x, p) for (p,x) ∈ zip(proj, X)]
    S, incl, proj
end

function tensor_product(X::TensorPowerObject, Y::TensorPowerObject)
    TensorPowerObject(parent(X), tensor_product(object(X), object(Y)))
end

function direct_sum(f::TensorPowerMorphism, g::TensorPowerMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    Morphism(dom,cod, direct_sum(morphism(f), morphism(g)))
end

function tensor_product(f::TensorPowerMorphism, g::TensorPowerMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)
    Morphism(dom, cod, tensor_product(morphism(f), morphism(g)))
end

function associator(X::TensorPowerObject, Y::TensorPowerObject, Z::TensorPowerObject)
    ass = associator(object.([X,Y,Z])...)
    dom = TensorPowerObject(parent(X), domain(ass))
    cod = TensorPowerObject(parent(X), codomain(ass))
    Morphism(dom,cod, ass)
end

inv(f::TensorPowerMorphism) = Morphism(codomain(f), domain(f), inv(morphism(f)))

one(C::TensorPowerCategory) = TensorPowerObject(C, one(category(C)))

function id(X::TensorPowerObject) 
    Morphism(X,X, id(object(X)))
end

#=----------------------------------------------------------
    Simples/Indecompodables 
----------------------------------------------------------=#

function indecomposable_subobjects(X::TensorPowerObject)
    subs = indecomposable_subobjects(object(X))
    return [TensorPowerObject(parent(X), s) for s ∈ subs]
end

# function simples(C::TensorPowerCategory, k = Inf)
#     simpls = object_type(category(C))[]
#     n1 = 0
#     j = 0
#     X = C.generator
#     Y = one(category(C))
#     while j ≤ k+1
#         simpls = unique_simples([simpls; simple_subobjects(Y)])
#         if length(simpls) == n1
#             simpls = [TensorPowerObject(C,s) for s ∈ simpls]
#             C.simples = simpls
#             C.complete = true
#             C.max_exponent = j-1
#             return simpls
#         end
#         n1 = length(simpls)
#         Y = Y ⊗ X
#         j = j+1
#     end
#     simpls = [TensorPowerObject(C,s) for s ∈ simpls]
#     C.simples = simpls
#     C.complete = false
#     C.max_exponent = k
#     return simpls
# end

function braiding(X::TensorPowerObject, Y::TensorPowerObject)
    b = braiding(object(X), object(Y))
    Morphism(X⊗Y,Y⊗X,b)
end

dual(X::TensorPowerObject) = TensorPowerObject(parent(X), dual(object(X)))

function ev(X::TensorPowerObject) 
    evaluation = ev(object(X))
    dom = TensorPowerObject(parent(X), domain(evaluation))
    cod = TensorPowerObject(parent(X), codomain(evaluation))
    Morphism(dom, cod, evaluation)
end

function coev(X::TensorPowerObject) 
    coevaluation = coev(object(X))
    dom = TensorPowerObject(parent(X), domain(coevaluation))
    cod = TensorPowerObject(parent(X), codomain(coevaluation))
    Morphism(dom, cod, coevaluation)
end

function spherical(X::TensorPowerObject)
    sp = spherical(object(X))
    dom = TensorPowerObject(parent(X), domain(sp))
    cod = TensorPowerObject(parent(X), codomain(sp))
    Morphism(dom, cod, sp)
end

zero(T::TensorPowerCategory) = TensorPowerObject(T, zero(parent(T.generator[1])))

function zero_morphism(X::TensorPowerObject, Y::TensorPowerObject)
    Morphism(X,Y, zero_morphism(object(X), object(Y)))
end

function indecomposables(C::TensorPowerCategory, k = Inf)
    if C.complete
        return C.simples
    end

   
    n1 = 0
    j = 2

    indecs_in_X = C.generator
    new_indecs = [one(category(C))]
    new_indecs = unique_indecomposables([new_indecs; indecs_in_X])
    simpls = new_indecs

    while j ≤ k
        new_indecs_temp = []
        for V ∈ indecs_in_X, W ∈ new_indecs
            summands_of_VW = [x for (x,k) ∈ decompose(W ⊗ V)]
            new_indecs_temp = [new_indecs_temp; [x for x ∈ summands_of_VW]]
            simpls = unique_indecomposables(Object[simpls; new_indecs_temp])
        end
        new_indecs = new_indecs_temp
        if length(simpls) == n1
            simpls = TensorPowerObject[TensorPowerObject(C,s) for s ∈ simpls]
            C.simples = simpls
            C.complete = true
            C.max_exponent = j-1
            return simpls
        end
        n1 = length(simpls)
        #Y = Y ⊗ X
        j = j+1
    end
    simpls = [TensorPowerObject(C,s) for s ∈ simpls]
    C.simples = simpls
    C.complete = false
    C.max_exponent = k

    return simpls
 
end

function decompose(X::TensorPowerObject)
    dec = decompose(object(X))
    [(TensorPowerObject(parent(X), x), k) for (x,k) ∈ dec]
end

function decompose(X::TensorPowerObject, S::Vector{TensorPowerObject})
    dec = decompose(object(X), object.(S))
    [(TensorPowerObject(parent(X), x), k) for (x,k) ∈ dec]
end

function kernel(f::TensorPowerMorphism)
    K,incl = kernel(morphism(f))
    K = TensorPowerObject(parent(f), K)
    return K, Morphism(K,domain(f), incl)
end

function cokernel(f::TensorPowerMorphism)
    C,proj = cokernel(morphism(f))
    C = TensorPowerObject(parent(f), C)
    return C, Morphism(codomain(f), C, proj)
end

function is_isomorphic(X::TensorPowerObject, Y::TensorPowerObject)
    is_iso, iso = is_isomorphic(object(X), object(Y))
    if is_iso
        return true, Morphism(X,Y,iso)
    else
        return false, nothing
    end
end

function show(io::IO, C::TensorPowerCategory)
    print(io, "Tensor power category with generator $(C.generator)")
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
    HomSpace(X,Y,B)
end