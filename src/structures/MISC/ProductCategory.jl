struct ProductCategory{N} <: Category
    factors::Tuple
end

struct ProductObject{N} <: Object
    parent::ProductCategory{N}
    factors::Tuple
end

struct ProductMorphism{N} <: Morphism
    domain::ProductObject{N}
    codomain::ProductObject{N}
    factors::Tuple
end


ProductCategory(C::Category...) = ProductCategory{length(C)}(C)
ProductObject(X::Object...) = ProductObject{length(X)}(ProductCategory(parent.(X)...), X)

×(C::Category, D::Category) = ProductCategory(C,D)

function Morphism(f::Morphism...)
    dom = ProductObject(domain.(f)...)
    cod = ProductObject(codomain.(f)...)
    ProductMorphism{length(X)}(dom,cod,f)
end

getindex(C::ProductCategory,x) = C.factors[x]
getindex(X::ProductObject,x) = X.factors[x]
getindex(f::ProductMorphism,x) = f.factors[x]

#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------

function dsum(X::ProductObject, Y::ProductObject)
    return ProductObject([dsum(x,y) for x ∈ X.factors, y ∈ Y.factors]...)
end

function tensor_product(X::ProductObject, Y::ProductObject)
    return ProductObject([tensor_product(x,y) for x ∈ X.factors, y ∈ Y.factors]...)
end

function dsum(f::ProductMorphism, g::ProductMorphism)
    ProductMorphism([dsum(fi,gi) for fi ∈ f.factors, gi ∈ g.factors])
end

function tensor_product(f::ProductMorphism, g::ProductMorphism)
    ProductMorphism([tensor_product(fi,gi) for fi ∈ f.factors, gi ∈ g.factors])
end

function simples(C::ProductCategory{N}) where N
    zeros = [zero(Ci) for Ci ∈ C.factors]
    simpls = ProductObject{N}[]
    for i ∈ 1:N
        for s ∈ simples C.factors[i]
            so = zeros
            so[i] = s
            push!(simpls, ProductObject(so...))
        end
    end
    return simpls
end
