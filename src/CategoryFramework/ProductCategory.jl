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

function Morphism(X::ProductObject{N}, Y::ProductObject{N}, f::NTuple{N,Morphism}) where N
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

getindex(X::ProductObject,x) = X.factors[x]
getindex(f::ProductMorphism,x) = f.factors[x]

function base_ring(C::ProductCategory) 
    all(e -> e == base_ring(C.factors[1]), base_ring.(C.factors)) && return base_ring(C.factors[1])
    parent(*(gen.(base_ring.([c for c ∈ C.factors]))...))
end

function matrix(f::ProductMorphism) 
    diagonal_matrix([change_base_ring(base_ring(f), matrix(fi)) for fi ∈ f.factors])
end

matrices(f::ProductMorphism) = vcat([matrices(fi) for fi ∈ f.factors]...)

inv(f::ProductMorphism) = ProductMorphism(codomain(f),domain(f), Tuple(inv(fi) for fi ∈ f.factors))

#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------


function direct_sum(X::ProductObject, Y::ProductObject)
    sums = [direct_sum(x,y) for (x,y) ∈ zip(X.factors,Y.factors)]
    Z = ProductObject(parent(X), Tuple(s[1] for s ∈ sums))
    ix = ProductMorphism(X,Z, Tuple(s[2][1] for s ∈ sums))
    iy = ProductMorphism(Y,Z, Tuple(s[2][2] for s ∈ sums))
    px = ProductMorphism(Z,X, Tuple(s[3][1] for s ∈ sums))
    py = ProductMorphism(Z,Y, Tuple(s[3][2] for s ∈ sums))
    return Z,[ix,iy],[px,py]
end


function tensor_product(X::ProductObject, Y::ProductObject)
    return ProductObject(parent(X), Tuple([tensor_product(x,y) for (x,y) ∈ zip(X.factors, Y.factors)]))
end

function direct_sum(f::ProductMorphism, g::ProductMorphism)
    ProductMorphism(domain(f)⊕domain(g), codomain(f)⊕codomain(g),Tuple([direct_sum(fi,gi) for (fi,gi) ∈ zip(f.factors, g.factors)]))
end

function tensor_product(f::ProductMorphism, g::ProductMorphism)
    ProductMorphism(domain(f)⊗domain(g), codomain(f)⊗codomain(g), Tuple([tensor_product(fi,gi) for (fi,gi) ∈ zip(f.factors, g.factors)]))
end

function simples(C::ProductCategory{N}) where N
    simpls = ProductObject{N}[]
    for i ∈ 1:N
        for s ∈ simples(C.factors[i])
            so = [zero(Ci) for Ci ∈ C.factors]
            so[i] = s
            push!(simpls, ProductObject(C,Tuple(so)))
        end
    end
    return simpls
end


function indecomposables(C::ProductCategory{N}) where N
    indecs = ProductObject{N}[]
    for i ∈ 1:N
        for s ∈ indecomposables(C.factors[i])
            so = [zero(Ci) for Ci ∈ C.factors]
            so[i] = s
            push!(indecs, ProductObject(C,Tuple(so)))
        end
    end
    return indecs
end



dual(X::ProductObject) = ProductObject(parent(X), Tuple(dual(x) for x ∈ X.factors))

ev(X::ProductObject) = ProductMorphism(dual(X)⊗X, one(parent(X)), Tuple(ev(x) for x ∈ X.factors))
coev(X::ProductObject) = ProductMorphism(one(parent(X)), X⊗dual(X), Tuple(coev(x) for x ∈ X.factors))

spherical(X::ProductObject) =  ProductMorphism(X,dual(dual(X)), Tuple(spherical(x) for x ∈ X.factors))

associator(X::ProductObject, Y::ProductObject, Z::ProductObject) = ProductMorphism((X⊗Y)⊗Z, X⊗(Y⊗Z), Tuple(associator(x,y,z) for (x,y,z) ∈ zip(X.factors,Y.factors,Z.factors)))

zero(C::ProductCategory) = ProductObject(C, zero.(C.factors))
one(C::ProductCategory) = ProductObject(C, one.(C.factors))

zero_morphism(X::ProductObject, Y::ProductObject) = ProductMorphism(X,Y, Tuple(zero_morphism(x,y) for (x,y) ∈ zip(X.factors, Y.factors)))
id(X::ProductObject) = ProductMorphism(X,X, Tuple(id(x) for x ∈ X.factors))

is_multitensor(C::ProductCategory) = *(is_multitensor.(C.factors)...)

function Hom(X::ProductObject{N}, Y::ProductObject{N}) where N
    basis = ProductMorphism[]
    for i ∈ 1:N
        for f ∈ Hom(X[i], Y[i])
            m_tuple = [zero_morphism(x,y) for (x,y) ∈ zip(X.factors, Y.factors)]
            m_tuple[i] = f
            basis = [basis; Morphism(X,Y,Tuple(m_tuple))]
        end
    end
    return HomSpace(X,Y,basis)
end
        
