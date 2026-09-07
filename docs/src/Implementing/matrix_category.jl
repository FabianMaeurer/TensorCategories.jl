module MatrixCategoryTutorial

using TensorCategories, Oscar

export MatCategory, MatObject

struct MatCategory <: Category
    base_ring::Field
end

struct MatObject <: Object
    parent::MatCategory
    n::Int
    function MatObject(C::MatCategory, n::Int)
        n >= 0 || throw(ArgumentError("dimension must be nonnegative"))
        new(C, n)
    end
end

struct MatMorphism <: Morphism
    domain::MatObject
    codomain::MatObject
    m::MatElem
end

Base.:(==)(C::MatCategory, D::MatCategory) = base_ring(C) === base_ring(D)
Base.:(==)(X::MatObject, Y::MatObject) = parent(X) == parent(Y) && X.n == Y.n
Base.:(==)(f::MatMorphism, g::MatMorphism) =
    domain(f) == domain(g) && codomain(f) == codomain(g) && f.m == g.m

function TensorCategories.morphism(X::MatObject, Y::MatObject, M::MatElem)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    base_ring(M) === base_ring(X) || throw(ArgumentError("different fields"))
    size(M) == (X.n, Y.n) || throw(ArgumentError("wrong matrix dimensions"))
    MatMorphism(X, Y, M)
end

TensorCategories.matrix(f::MatMorphism) = f.m
TensorCategories.int_dim(X::MatObject) = X.n
TensorCategories.id(X::MatObject) = morphism(X, X, identity_matrix(base_ring(X), X.n))

function TensorCategories.compose(f::MatMorphism, g::MatMorphism)
    codomain(f) == domain(g) || throw(ArgumentError("incompatible endpoints"))
    morphism(domain(f), codomain(g), matrix(f)*matrix(g))
end

function Base.:+(f::MatMorphism, g::MatMorphism)
    domain(f) == domain(g) && codomain(f) == codomain(g) ||
        throw(ArgumentError("maps must be parallel"))
    morphism(domain(f), codomain(f), matrix(f) + matrix(g))
end

Base.:*(a, f::MatMorphism) =
    morphism(domain(f), codomain(f), base_ring(f)(a)*matrix(f))
TensorCategories.zero_morphism(X::MatObject, Y::MatObject) =
    morphism(X, Y, zero_matrix(base_ring(X), X.n, Y.n))
Base.zero(C::MatCategory) = MatObject(C, 0)

function TensorCategories.Hom(X::MatObject, Y::MatObject)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    B = MatMorphism[]
    for j in 1:Y.n, i in 1:X.n
        M = zero_matrix(base_ring(X), X.n, Y.n)
        M[i,j] = 1
        push!(B, morphism(X, Y, M))
    end
    HomSpace(X, Y, B)
end

function TensorCategories.direct_sum(X::MatObject, Y::MatObject)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    D = MatObject(parent(X), X.n + Y.n)
    K = base_ring(X)
    ix = zero_matrix(K, X.n, D.n)
    iy = zero_matrix(K, Y.n, D.n)
    for j in 1:X.n
        ix[j,j] = 1
    end
    for j in 1:Y.n
        iy[j,X.n+j] = 1
    end
    i = [morphism(X,D,ix), morphism(Y,D,iy)]
    p = [morphism(D,X,transpose(ix)), morphism(D,Y,transpose(iy))]
    D, i, p
end

# Reuse the package's vector-space linear algebra, then reconstruct our objects.
function TensorCategories.kernel(f::MatMorphism)
    V, i = kernel(morphism(matrix(f)))
    K = MatObject(parent(f), int_dim(V))
    K, morphism(K, domain(f), matrix(i))
end

function TensorCategories.cokernel(f::MatMorphism)
    V, p = cokernel(morphism(matrix(f)))
    Q = MatObject(parent(f), int_dim(V))
    Q, morphism(codomain(f), Q, matrix(p))
end

function TensorCategories.tensor_product(X::MatObject, Y::MatObject)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    MatObject(parent(X), X.n*Y.n)
end
TensorCategories.tensor_product(f::MatMorphism, g::MatMorphism) =
    morphism(domain(f)⊗domain(g), codomain(f)⊗codomain(g),
             kronecker_product(matrix(f), matrix(g)))
Base.one(C::MatCategory) = MatObject(C, 1)
TensorCategories.associator(X::MatObject, Y::MatObject, Z::MatObject) =
    id((X⊗Y)⊗Z)

# These declarations describe the model just implemented; they do not test it.
TensorCategories.is_linear(::MatCategory) = true
TensorCategories.is_abelian(::MatCategory) = true
TensorCategories.is_monoidal(::MatCategory) = true

end
