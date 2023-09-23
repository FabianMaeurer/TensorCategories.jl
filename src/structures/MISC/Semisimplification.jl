#=----------------------------------------------------------
    Construct a semisimplification of any
    tensor category.
    Reference: https://doi.org/10.48550/arXiv.1801.04409 
----------------------------------------------------------=#

struct Semisimplification <: Category
    category::Category
end

struct SemisimplifiedObject <: Object
    parent::Semisimplification
    object::Object
end

struct SemisimplifiedMorphism <: Morphism 
    domain::SemisimplifiedObject
    codomain::SemisimplifiedObject
    morphism::morphism
end

function Morphism(dom::SemisimplifiedObject, cod::SemisimplifiedObject, m::Morphism) 
    SemisimplifiedMorphism(dom,cod, m)
end

morphism(f::SemisimplifiedMorphism) = f.morphism
object(X::SemisimplifiedObject) = X.object
category(C::Semisimplification) = C.category

#=----------------------------------------------------------
    Morphism functionality 
----------------------------------------------------------=#

function compose(f::SemisimplifiedMorphism...)
    domain = domain(f[1])
    codomain = codomain(f[end])
    Morphism(domain, codomain, compose(morphism.(f)...))
end

function direct_sum(f::SemisimplifiedMorphism...)
    domain = direct_sum(domain.(f)...)
    codomain = direct_sum(codomain.(f)...)
    map = direct_sum(morphism.(f)...)
    Morphism(domain, codomain, map)
end

function tensor_product(f::SemisimplifiedMorphism...)
    domain = tensor_product(domain.(f)...)
    codomain = tensor_product(codomain.(f)...)
    map = tensor_product(morphism.(f)...)
    Morphism(domain, codomain, map)
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
    base = basis(Hom(object(X), object(Y)))
