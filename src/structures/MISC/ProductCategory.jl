struct ProductCategory{N} <: Category
    factors::Tuple
end

struct ProductCategoryObject{N} <: CategoryObject
    parent::ProductCategory{N}
    factors::Tuple
end

struct ProductMorphism{N} <: CategoryMorphism
    domain::ProductCategoryObject{N}
    codomain::ProductCategoryObject{N}
    factors::Tuple
end


ProductCategory(C::Category...) = ProductCategory{length(C)}(C)
ProductCategoryObject(X::CategoryObject...) = ProductCategoryObject{length(X)}(ProductCategory(parent.(X)...), X)

×(C::Category, D::Category) = ProductCategory(C,D)

function Morphism(f::CategoryMorphism...)
    dom = ProductCategoryObject(domain.(f)...)
    cod = ProductCategoryObject(codomain.(f)...)
    ProductMorphism{length(X)}(dom,cod,f)
end

function Morphism(X::ProductCategoryObject{N}, Y::ProductCategoryObject{N}, f::NTuple{N,CategoryMorphism}) where N
    ProductMorphism{N}(X,Y,f)
end

function compose(f::ProductMorphism, g::ProductMorphism)
    return ProductMorphism(domain(f), codomain(g), Tuple(compose(fi,gi) for (fi,gi) ∈ zip(f.factors, g.factors)))
end

function +(f::ProductMorphism, g::ProductMorphism)
    return ProductMorphism(domain(f), codomain(g), Tuple(fi+gi for (fi,gi) ∈ zip(f.factors, g.factors)))
end

function *(λ, f::ProductMorphism)
    return ProductMorphism(domain(f), codomain(f), Tuple(λ*fi for (fi) ∈ f.factors))
end

getindex(X::ProductCategoryObject,x) = X.factors[x]
getindex(f::ProductMorphism,x) = f.factors[x]

base_ring(C::ProductCategory) = parent(*(gen.(base_ring.([c for c ∈ C.factors]))...))

function matrix(f::ProductMorphism) 
    diagonal_matrix([change_base_ring(base_ring(f), matrix(fi)) for fi ∈ f.factors])
end

matrices(f::ProductMorphism) = vcat([matrices(fi) for fi ∈ f.factors]...)

inv(f::ProductMorphism) = ProductMorphism(codomain(f),domain(f), Tuple(inv(fi) for fi ∈ f.factors))

#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------


function direct_sum(X::ProductCategoryObject, Y::ProductCategoryObject)
    sums = [direct_sum(x,y) for (x,y) ∈ zip(X.factors,Y.factors)]
    Z = ProductCategoryObject(parent(X), Tuple(s[1] for s ∈ sums))
    ix = ProductMorphism(X,Z, Tuple(s[2][1] for s ∈ sums))
    iy = ProductMorphism(Y,Z, Tuple(s[2][2] for s ∈ sums))
    px = ProductMorphism(Z,X, Tuple(s[3][1] for s ∈ sums))
    py = ProductMorphism(Z,Y, Tuple(s[3][2] for s ∈ sums))
    return Z,[ix,iy],[px,py]
end


function tensor_product(X::ProductCategoryObject, Y::ProductCategoryObject)
    return ProductCategoryObject(parent(X), Tuple([tensor_product(x,y) for (x,y) ∈ zip(X.factors, Y.factors)]))
end

function direct_sum(f::ProductMorphism, g::ProductMorphism)
    ProductMorphism(domain(f)⊕domain(g), codomain(f)⊕codomain(g),Tuple([direct_sum(fi,gi) for (fi,gi) ∈ zip(f.factors, g.factors)]))
end

function tensor_product(f::ProductMorphism, g::ProductMorphism)
    ProductMorphism(domain(f)⊗domain(g), codomain(f)⊗codomain(g), Tuple([tensor_product(fi,gi) for (fi,gi) ∈ zip(f.factors, g.factors)]))
end

function simples(C::ProductCategory{N}) where N
    simpls = ProductCategoryObject{N}[]
    for i ∈ 1:N
        for s ∈ simples(C.factors[i])
            so = [zero(Ci) for Ci ∈ C.factors]
            so[i] = s
            push!(simpls, ProductCategoryObject(C,Tuple(so)))
        end
    end
    return simpls
end


dual(X::ProductCategoryObject) = ProductCategoryObject(parent(X), Tuple(dual(x) for x ∈ X.factors))

ev(X::ProductCategoryObject) = ProductMorphism(dual(X)⊗X, one(parent(X)), Tuple(ev(x) for x ∈ X.factors))
coev(X::ProductCategoryObject) = ProductMorphism(one(parent(X)), X⊗dual(X), Tuple(coev(x) for x ∈ X.factors))

spherical(X::ProductCategoryObject) =  ProductMorphism(X,dual(dual(X)), Tuple(spherical(x) for x ∈ X.factors))

associator(X::ProductCategoryObject, Y::ProductCategoryObject, Z::ProductCategoryObject) = ProductMorphism((X⊗Y)⊗Z, X⊗(Y⊗Z), Tuple(associator(x,y,z) for (x,y,z) ∈ zip(X.factors,Y.factors,Z.factors)))

zero(C::ProductCategory) = ProductCategoryObject(C, zero.(C.factors))
one(C::ProductCategory) = ProductCategoryObject(C, one.(C.factors))

zero_morphism(X::ProductCategoryObject, Y::ProductCategoryObject) = ProductMorphism(X,Y, Tuple(zero_morphism(x,y) for (x,y) ∈ zip(X.factors, Y.factors)))
id(X::ProductCategoryObject) = ProductMorphism(X,X, Tuple(id(x) for x ∈ X.factors))

is_multitensor(C::ProductCategory) = *(is_multitensor.(C.factors)...)

function Hom(X::ProductCategoryObject{N}, Y::ProductCategoryObject{N}) where N
    basis = ProductMorphism[]
    for i ∈ 1:N
        for f ∈ Hom(X[i], Y[i])
            m_tuple = [zero_morphism(x,y) for (x,y) ∈ zip(X.factors, Y.factors)]
            m_tuple[i] = f
            basis = [basis; Morphism(X,Y,Tuple(m_tuple))]
        end
    end
    return CategoryHomSpace(X,Y,basis, VectorSpaces(base_ring(X)))
end
        
