struct ProductCategory{N} <: Category
    factors::Tuple
end

struct ProductObject{N} <: Object
    factors::Tuple
end

struct ProductMorphism{N} <: Morphism
    factors::Tuple
end


ProductCategory(C::Category...) = ProductCategory{length(C)}(C)
ProductObject(X::Object...) = ProductObject{length(X)}(X)
ProductMorphism(f::Morphism...) = ProductMorphism{length(X)}(f)

domain(f::ProductMorphism) = ProductObject(domain(f.factors))
codomain(f::ProductMorphism) = ProductObject(codomain(f.factors))

getindex(C::ProductCategory,x) = C.factors[x]
getindex(X::ProductObject,x) = X.factors[x]
getindex(f::ProductMorphism,x) = f.factors[x] 
