#=----------------------------------------------------------
    Construct a semisimplification of any
    tensor category.
    Reference: https://doi.org/10.48550/arXiv.1801.04409 
----------------------------------------------------------=#

mutable struct Semisimplification <: Category
    category::Category
    simples::Vector{Object}

    Semisimplification(C::Category) = new(C)
end

struct SemisimplifiedObject <: Object
    parent::Semisimplification
    object::Object
end

struct SemisimplifiedMorphism <: Morphism 
    domain::SemisimplifiedObject
    codomain::SemisimplifiedObject
    morphism::Morphism
end

function Morphism(dom::SemisimplifiedObject, cod::SemisimplifiedObject, m::Morphism) 
    SemisimplifiedMorphism(dom,cod, m)
end


function ==(X::Semisimplification, Y::Semisimplification)
    category(X) == category(Y)
end
    
morphism(f::SemisimplifiedMorphism) = f.morphism
object(X::SemisimplifiedObject) = X.object
category(C::Semisimplification) = C.category
base_ring(C::Semisimplification) = base_ring(category(C))

function Semisimplification(X::Object, C::Semisimplification) 
    @assert parent(X) == parent(C.category)
    SemisimplifiedObject(C,X)
end

function Semisimplification(X::Object)
    C = Semisimplification(parent(X))
    SemisimplifiedObject(C,X)
end

function Semisimplification(f::Morphism, C::Semisimplification)
    dom = SemisimplifiedObject(C, domain(f))
    cod = SemisimplifiedObject(C, codomain(f))
    Morphism(dom, cod, f)
end

function Semisimplification(f::Morphism)
    C = Semisimplification(parent(domain(f)))
    Semisimplification(f,C)
end

is_semisimple(C::Semisimplification) = true
is_multiring(C::Semisimplification) = is_multiring(category(C))
is_multifusion(C::Semisimplification) = is_multiring(category(C)) && is_rigid(category(C))

semisimplify(C::Category) = Semisimplification(C)
semisimplify(X::Object) = SemisimplifiedObject(Semisimplification(parent(X)), X)
semisimplify(X::Object, C::Category) = SemisimplifiedObject(C,X)

function semisimplify(f::Morphism, C::Semisimplification)   
    dom = semisimplify(domain(f),C)
    cod = semisimplify(codomain(f),C)
    SemisimplifiedMorphism(dom,cod,f)
end

semisimplify(f::Morphism) = SemisimplifiedMorphism(f,semisimplify(parent(f)))

#=----------------------------------------------------------
    Morphism functionality 
----------------------------------------------------------=#

function compose(f::SemisimplifiedMorphism...)
    dom = domain(f[1])
    codom = codomain(f[end])
    Morphism(dom, codom, compose(morphism.(f)...))
end

function direct_sum(f::SemisimplifiedMorphism...)
    dom = direct_sum(domain.(f)...)
    codom = direct_sum(codomain.(f)...)
    map = direct_sum(morphism.(f)...)
    Morphism(dom, codom, map)
end

function tensor_product(f::SemisimplifiedMorphism...)
    dom = tensor_product(domain.(f)...)
    codom = tensor_product(codomain.(f)...)
    map = tensor_product(morphism.(f)...)
    Morphism(dom, codom, map)
end

function associator(X::SemisimplifiedObject, Y::SemisimplifiedObject, Z::SemisimplifiedObject)
    C = parent(X)
    a = associator(object.([X,Y,Z])...)
    dom = SemisimplifiedObject(C, domain(a))
    cod = SemisimplifiedObject(C, codomain(a))
    SemisimplifiedMorphism(dom,cod, a)
end

function is_negligible(f::Morphism)
    X = domain(f)
    Y = codomain(f)
    iszero([tr(f ∘ g) == zero_morphism(Y,Y) for g ∈ Hom(Y,X)])
end

function ==(f::SemisimplifiedMorphism, g::SemisimplifiedMorphism)
    is_negligible(f-g)
end

function Hom(X::SemisimplifiedObject, Y::SemisimplifiedObject)
    base_XY = basis(Hom(object(X), object(Y)))
    base_YX = basis(Hom(object(Y), object(X)))


    F = base_ring(X)

    if length(base_XY) == 0
        return HomSpace(X,Y,SemisimplifiedMorphism[])
    end

    if length(base_YX) == 0
        base = base_XY
    else
        # linear system to find negligible morphisms
        # (tr(f ∘ g) = 0 for all g)
        M = zero_matrix(F, length(base_YX), length(base_XY))
        

        m = [F(tr(f ∘ g)) for f ∈ base_XY, g ∈ base_YX]
        
        M = matrix(F, length(base_XY), length(base_YX), m)

        # The image of M except 0 are all non-negligible morphisms
        base_coeffs = hnf(M)
        #base_coeffs = base_coeffs[1:rank(base_coeffs), 1]
        
        # Basis 
        base =  [sum(collect(base_coeffs[i,:])[:] .* base_XY) for i ∈ 1:length(base_coeffs[:,1])]
    end
    filter!(e -> e != zero_morphism(object(X), object(Y)), base)

    base = [SemisimplifiedMorphism(X,Y, f) for f ∈ base]
    HomSpace(X,Y,base)
end

function id(X::SemisimplifiedObject)
    Morphism(X,X, id(object(X)))
end

function matrix(f::SemisimplifiedMorphism)
    matrix(morphism(f))
end

function zero_morphism(X::SemisimplifiedObject, Y::SemisimplifiedObject)
    Morphism(X,Y, zero_morphism(object(X), object(Y)))
end

function kernel(f::SemisimplifiedMorphism)
    SC = parent(f)
    K,k = kernel(morphism(f))
    SK = SemisimplifiedObject(SC, K)
    SK, SemisimplifiedMorphism(SK, domain(f), k)
end

function cokernel(f::SemisimplifiedMorphism)
    SC = parent(f)
    C,c = cokernel(morphism(f))
    S = SemisimplifiedObject(SC, C)
    S, SemisimplifiedMorphism(codomain(f), S, c)
end
 


#=----------------------------------------------------------
    Object Functionality 
----------------------------------------------------------=#

function direct_sum(X::SemisimplifiedObject...)
    C = parent(X[1])
    S,i,p = direct_sum(object.(X)...)
    S = SemisimplifiedObject(C, S)
    i = [semisimplify(f, C) for f ∈ i]
    p = [semisimplify(f, C) for f ∈ p]
    S,i,p
end

function tensor_product(X::SemisimplifiedObject...)
    T = tensor_product(object.(X)...)
    SemisimplifiedObject(parent(X[1]), T)
end

function zero(C::Semisimplification)
    SemisimplifiedObject(C, zero(category(C)))
end

function one(C::Semisimplification)
    SemisimplifiedObject(C, one(category(C)))
end

function simples(C::Semisimplification)
    indecs = filter(e -> dim(e) != 0, indecomposables(category(C)))
    simpls = [semisimplify(x,C) for x ∈ indecs]
    C.simples = unique(simpls)
end

function inv(f::SemisimplifiedMorphism)
    Morphism(codomain(f), domain(f), inv(morphism(f)))
end

spherical(X::SemisimplifiedObject) = Morphism(X,dual(dual(X)), spherical(object(X)))

dual(X::SemisimplifiedObject) = SemisimplifiedObject(parent(X), dual(object(X)))

function ev(X::SemisimplifiedObject)
    e = ev(object(X))
    C = parent(X)
    dom = SemisimplifiedObject(C, domain(e))
    cod = SemisimplifiedObject(C, codomain(e))
    Morphism(dom, cod, e)
end

function coev(X::SemisimplifiedObject)
    c = coev(object(X))
    C = parent(X)
    dom = SemisimplifiedObject(C, domain(c))
    cod = SemisimplifiedObject(C, codomain(c))
    Morphism(dom, cod, c)
end

#=----------------------------------------------------------
    Printing 
----------------------------------------------------------=#

function show(io::IO, X::SemisimplifiedObject)
    if int_dim(End(X)) == 0 
        print(io, "Semisimplified: $(zero(parent(object(X))))")
    else
        print(io, "Semisimplified: $(object(X))")
    end
end

function show(io::IO, C::Semisimplification)
    print(io, "Semisimplification of $(category(C))")
end