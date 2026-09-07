# [Linear and abelian categories](@id linear-categories)

Fix a coefficient field $k$. More information about constructing exact and
numerical coefficient fields is given in the [next section](@ref base-fields).
Here we first describe the mathematical structures and the corresponding Julia
interfaces. Only afterwards do we implement a complete example.

## Linear categories

A $k$-linear structure on a category $\mathcal C$ consists of a
$k$-vector-space structure on every $\operatorname{Hom}_{\mathcal C}(X,Y)$
such that composition

```math
\label{eq:linear-category-composition}
\operatorname{Hom}_{\mathcal C}(Y,Z)\times
\operatorname{Hom}_{\mathcal C}(X,Y)
\longrightarrow \operatorname{Hom}_{\mathcal C}(X,Z)
```

is $k$-bilinear. We call a category equipped with such a structure
$k$-linear. In particular, its Hom sets are abelian groups, so a $k$-linear
category is preadditive.

!!! note "Terminology"
    [EGNO; Definition 1.2.2](@citet) begins with an additive category and then
    calls it $k$-linear when its Hom groups carry compatible $k$-vector-space
    structures. Here, as in the more general enrichment convention, *linear*
    refers only to the Hom-space structure and bilinear composition; finite
    biproducts are recorded separately as additivity. Thus a category that is
    $k$-linear in the sense of EGNO is both linear and additive in the
    terminology used by the interface.

### What an implementation must provide

A finite-dimensional $k$-linear model should provide the following methods:

| Required method | Meaning |
|:---|:---|
| `base_ring(C)` | the coefficient field $k$ |
| `Hom(X,Y)` | a representation of $\operatorname{Hom}_{\mathcal C}(X,Y)$ |
| `basis(H)` | a finite ordered basis of the represented Hom space |
| `f + g` | addition of parallel morphisms |
| `a*f` | scalar multiplication for $a\in k$ |
| `zero_morphism(X,Y)` | the zero morphism $X\to Y$ |
| `express_in_basis(f,H)` | the coordinates of $f$ in the ordered basis of $H$ |

The last method can be obtained from a faithful matrix realization, but a
matrix realization is not required: a model may implement coordinates by any
valid method. After providing these operations and checking bilinearity, the
model reports the structure through `is_linear(C)`.

### Generic consequences

Several functions are then available from the generic interface:

| Function | Generic meaning or construction |
|:---|:---|
| `End(X)` | `Hom(X,X)` |
| `int_dim(H)` | the dimension of a standard `HomSpace` from its stored basis |
| `f - g` and `-f` | subtraction using addition and scalar multiplication |
| `is_zero(f)` | comparison with `zero_morphism(domain(f),codomain(f))` |
| `endomorphism_ring(X)` | the $k$-algebra $\operatorname{End}_{\mathcal C}(X)$, using a finite basis and coordinates |

The generic functions have hypotheses. For example, `endomorphism_ring(X)`
needs an effective finite Hom basis and a way to express composites in that
basis. Merely declaring `is_linear(C) = true` does not create those operations
or verify their axioms.

## Additive categories

A category is *semiadditive* if it has finite biproducts. It is *additive* if
it is both semiadditive and preadditive; equivalently, it has zero morphisms,
finite biproducts, and compatible abelian-group structures on its Hom sets. A
$k$-linear category is already preadditive, so it becomes additive once finite
biproducts are supplied. This is equivalent to the axioms in
[EGNO; Definition 1.2.1](@citet).

### What an implementation must provide

For an additive model, implement `zero(C)` and binary `direct_sum(X,Y)`. The
latter returns

```julia
D, i, p = direct_sum(X, Y)
```

where $D=X\oplus Y$, `i` contains the inclusions, and `p` contains the
projections. For $D=\bigoplus_{s=1}^nX_s$, these maps must satisfy

```math
\label{eq:biproduct-identities}
p_r\circ i_s=
\begin{cases}
\operatorname{id}_{X_s},&r=s,\\
0_{X_s,X_r},&r\ne s,
\end{cases}
\qquad
\sum_r i_r\circ p_r=\operatorname{id}_D.
```

An implementation not already carrying a linear structure must also provide
the preadditive Hom-group operations. Once the axioms are established, it
reports the property with `is_additive(C)`.

### Generic consequences

From binary direct sums and a zero object, the interface derives:

| Function | Generic construction |
|:---|:---|
| `direct_sum(X₁,...,Xₙ)` | iterated binary biproduct, with all inclusions and projections |
| `X ⊕ Y` | the biproduct object without its structural maps |
| `product(X,Y)` and `coproduct(X,Y)` | the biproduct with projections or inclusions |
| `initial_object(C)` and `terminal_object(C)` | the zero object |

The empty direct sum is `zero(C)`. An empty collection of objects does not by
itself determine its parent category.

## [Abelian categories](@id abelian-operations)

An abelian category is an additive category with kernels and cokernels in which
every morphism has the usual image--coimage factorization; equivalently, every
monomorphism is a kernel and every epimorphism is a cokernel. We use the
conventions of [EGNO; Definition 1.3.1](@citet).

### What an implementation must provide

For every $f:X\to Y$, an abelian model must implement:

| Required method | Return value | Universal property begins with |
|:---|:---|:---|
| `kernel(f)` | $(K,i)$ with $i:K\to X$ | $f\circ i=0$ |
| `cokernel(f)` | $(Q,p)$ with $p:Y\to Q$ | $p\circ f=0$ |

The displayed zero composites do not suffice: $i$ and $p$ must satisfy the
kernel and cokernel universal properties. The implementation is also
responsible for the abelian normality conditions, which are not verified by
dispatch. Once these conditions are known, it reports the property through
`is_abelian(C)`.

### Generic consequences

The generic interface then provides:

| Function | Generic construction |
|:---|:---|
| `image(f)` | the kernel of the cokernel of $f$ |
| `is_monomorphism(f)` | tests whether the kernel object is zero |
| `is_epimorphism(f)` | tests whether the cokernel object is zero |

These functions rely on the abelian-category contract. A matrix nullspace is
not yet a categorical kernel: the implementation must reconstruct an object of
the category and the corresponding universal morphism.

## [Example: Matrix category](@id implementing-matrices)

Let $\operatorname{Mat}_k$ be the category whose objects are the nonnegative
integers. The object $n$ represents the standard coordinate space $k^n$, and

```math
\label{eq:matrix-category-hom}
\operatorname{Hom}_{\operatorname{Mat}_k}(n,m)
=\operatorname{Mat}_{n\times m}(k).
```

We use row coordinates, so a matrix in equation
\eqref{eq:matrix-category-hom} acts from the left on row vectors. Identities
are identity matrices and

```math
\label{eq:matrix-category-composition}
M_{g\circ f}=M_fM_g.
```

This is a skeletal, or ``discrete'', coordinate model of the category
$\operatorname{Vec}_k$ of finite-dimensional vector spaces: every vector space
is replaced by the standard space of its dimension. It is not a discrete
category in the categorical sense, since it has many nonidentity morphisms.

The following blocks form one continuous Julia session. We begin with the
represented data and the basic category interface:

```@example matrix_category_tutorial
using TensorCategories, Oscar

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
    matrix::MatElem
end

Base.:(==)(C::MatCategory, D::MatCategory) = base_ring(C) === base_ring(D)
Base.:(==)(X::MatObject, Y::MatObject) =
    parent(X) == parent(Y) && X.n == Y.n
Base.:(==)(f::MatMorphism, g::MatMorphism) =
    domain(f) == domain(g) && codomain(f) == codomain(g) && matrix(f) == matrix(g)

function TensorCategories.morphism(X::MatObject, Y::MatObject, M::MatElem)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    base_ring(M) === base_ring(X) || throw(ArgumentError("different fields"))
    size(M) == (X.n, Y.n) || throw(ArgumentError("wrong matrix dimensions"))
    MatMorphism(X, Y, M)
end

TensorCategories.matrix(f::MatMorphism) = f.matrix
TensorCategories.int_dim(X::MatObject) = X.n
TensorCategories.id(X::MatObject) =
    morphism(X, X, identity_matrix(base_ring(X), X.n))

function TensorCategories.compose(f::MatMorphism, g::MatMorphism)
    codomain(f) == domain(g) || throw(ArgumentError("incompatible endpoints"))
    morphism(domain(f), codomain(g), matrix(f)*matrix(g))
end
```

The field name `base_ring` and the endpoint field names activate the generic
accessors. We now add the linear interface. The basis of a Hom space consists
of the elementary matrices.

```@example matrix_category_tutorial
function Base.:+(f::MatMorphism, g::MatMorphism)
    domain(f) == domain(g) && codomain(f) == codomain(g) ||
        throw(ArgumentError("maps must be parallel"))
    morphism(domain(f), codomain(f), matrix(f) + matrix(g))
end

Base.:*(a, f::MatMorphism) =
    morphism(domain(f), codomain(f), base_ring(f)(a)*matrix(f))

TensorCategories.zero_morphism(X::MatObject, Y::MatObject) =
    morphism(X, Y, zero_matrix(base_ring(X), X.n, Y.n))

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

TensorCategories.is_linear(::MatCategory) = true
```

The additive structure uses the zero-dimensional object and the standard block
inclusions and projections:

```@example matrix_category_tutorial
Base.zero(C::MatCategory) = MatObject(C, 0)

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

TensorCategories.is_additive(::MatCategory) = true
```

Finally, matrix nullspaces give kernels and cokernels. The methods turn the
resulting matrices back into objects and morphisms of our category:

```@example matrix_category_tutorial
function TensorCategories.kernel(f::MatMorphism)
    M = kernel(matrix(f))
    K = MatObject(parent(domain(f)), number_of_rows(M))
    K, morphism(K, domain(f), M)
end

function TensorCategories.cokernel(f::MatMorphism)
    M = kernel(matrix(f), side=:right)
    Q = MatObject(parent(f), number_of_columns(M))
    Q, morphism(codomain(f), Q, M)
end

TensorCategories.is_abelian(::MatCategory) = true
```

We can now use generic linear, additive, and abelian functions on this model:

```@example matrix_category_tutorial
C = MatCategory(QQ)
X, Y = MatObject(C, 2), MatObject(C, 3)
f = morphism(X, Y, matrix(QQ, [1 0 0; 0 0 0]))

H = Hom(X,Y)
@assert int_dim(H) == 6
@assert int_dim(End(X)) == 4

D, i, p = direct_sum(X,Y)
@assert p[1] ∘ i[1] == id(X) && p[2] ∘ i[2] == id(Y)
@assert is_zero(p[1] ∘ i[2]) && is_zero(p[2] ∘ i[1])
@assert i[1] ∘ p[1] + i[2] ∘ p[2] == id(D)

K, inclusion = kernel(f)
Q, projection = cokernel(f)
@assert int_dim(K) == 1 && int_dim(Q) == 2
@assert is_zero(f ∘ inclusion) && is_zero(projection ∘ f)

I, image_inclusion = image(f)
@assert int_dim(I) == 1
(int_dim(H), int_dim(D), int_dim(K), int_dim(Q), int_dim(I))
```

## [Matrix realizations and built-in models](@id matrix-realizations)

The method `matrix(f)` is meaningful only when a category supplies compatible
coordinates. Frequently these coordinates come from a faithful $k$-linear
functor

```math
\label{eq:matrix-realization-functor}
U:\mathcal C\longrightarrow\operatorname{Vec}_k
```

that is implicit in the stored data. A model need not store $U$ as a separate
Julia functor. It must nevertheless document the bases and the direction in
which its matrices act.

`vector_spaces(k)` is the built-in version of the matrix category above:

```@example built_in_linear_models
using TensorCategories, Oscar
V = vector_spaces(QQ)
X = VectorSpaceObject(V, 2)
Y = VectorSpaceObject(V, 3)
f = morphism(X, Y, matrix(QQ, [1 0 2; 0 1 3]))
@assert size(matrix(f)) == (2,3)
@assert int_dim(Hom(X,Y)) == 6
matrix(f)
show(stdout, MIME"text/plain"(), matrix(f)); println() # hide
```

## [Example: Group representations](@id concrete-models)

Let $G$ be a finite group. The category $\operatorname{Rep}_k(G)$ is another
$k$-linear abelian category, but not every matrix between underlying vector
spaces is a morphism. With the package's row-vector convention, a matrix
$M:X\to Y$ is a morphism precisely when

```math
\label{eq:representation-intertwiner-foundations}
\rho_X(g)M=M\rho_Y(g)
```

for every chosen generator $g$ of $G$. The constructor
`representation_category(k,G)` stores this model, and `Hom(X,Y)` solves these
intertwining equations. Forgetting the action gives the faithful functor in
equation \eqref{eq:matrix-realization-functor}, but it is generally not full.

```@example built_in_linear_models
G = cyclic_group(3)
R = representation_category(QQ, G)
A = matrix(QQ, [0 1; -1 -1])
X = Representation(R, gens(G), [A]; check=true)
@assert A^3 == identity_matrix(QQ, 2)
@assert int_dim(X) == 2
@assert int_dim(End(X)) == 2
int_dim(End(X))
```

The coefficient field affects the solutions of the intertwining equations and
the resulting abelian category. We therefore discuss coefficient fields before
turning to simple objects and composition factors.

Continue with [coefficient fields and numeric computations](@ref base-fields).
