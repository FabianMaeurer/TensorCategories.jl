# [Linear and abelian categories](@id linear-categories)

Fix a coefficient field $k$. More information about constructing exact and
numerical coefficient fields is given in the [next section](@ref base-fields).
Here we first describe the mathematical structures and the corresponding Julia
interfaces. We then illustrate them with the built-in categories of vector
spaces and finite-group representations.

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

### The interface

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

The last method can be obtained from a matrix realization, but matrices are not
part of the definition of a linear category: a model may implement coordinates
by any valid method. After providing these operations and checking
bilinearity, the model reports the structure through `is_linear(C)`.

#### [Matrix realizations](@id matrix-realizations)

A **matrix realization** of a $k$-linear category is a faithful $k$-linear
functor

```math
\label{eq:matrix-realization-functor}
U:\mathcal C\longrightarrow\operatorname{Vec}_k
```

together with an ordered basis of every vector space $U(X)$. The implementation
need not store $U$ as a Julia functor; it may be implicit in the data used to
represent objects and morphisms. If the model supplies `matrix(f)` for
$f:X\to Y$, it must return the matrix $M_f$ of $U(f)$ in these chosen bases.

TensorCategories.jl uses row coordinates in its concrete matrix models. Thus
`matrix(f)` has entries in $k$ and size
$\dim_k U(X)\times\dim_k U(Y)$, a row vector $v\in U(X)$ is sent to
$vM_f$, and the matrices must satisfy

```math
\label{eq:matrix-realization-composition}
M_{\operatorname{id}_X}=I,
\qquad
M_{g\circ f}=M_fM_g,
\qquad
M_{af+bg}=aM_f+bM_g.
```

Faithfulness means that two parallel morphisms are equal whenever their
matrices are equal. A model providing `matrix(f)` must document the underlying
realization, the ordered bases, and the direction in which its matrices act.
No compatibility with a monoidal structure is assumed here; that additional
condition belongs to the later notion of a fiber functor.

#### Generic functions

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

### The interface

An additive model should provide the following methods in addition to the
preadditive Hom-group operations:

| Required method | Meaning |
|:---|:---|
| `zero(C)` | a chosen zero object of $\mathcal C$ |
| `direct_sum(X,Y)` | the binary biproduct $X\oplus Y$ with its inclusions and projections |

The second method returns

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

#### Generic functions

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

### The interface

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

#### Generic functions

The generic interface then provides:

| Function | Generic construction |
|:---|:---|
| `image(f)` | the kernel of the cokernel of $f$ |
| `is_monomorphism(f)` | tests whether the kernel object is zero |
| `is_epimorphism(f)` | tests whether the cokernel object is zero |

These functions rely on the abelian-category contract. A matrix nullspace is
not yet a categorical kernel: the implementation must reconstruct an object of
the category and the corresponding universal morphism.

## [Working with built-in abelian categories](@id built-in-abelian-categories)

### Vector spaces

The category $\operatorname{Vec}_k$ has finite-dimensional $k$-vector spaces
as objects and linear maps as morphisms; see
[EGNO; Example 2.3.3, p. 26](@citet). TensorCategories.jl uses the coordinate
model constructed by `vector_spaces(k)`. An object created by
`VectorSpaceObject(V,n)` represents $k^n$ with its standard ordered basis, and
`morphism(X,Y,M)` constructs the linear map whose row-coordinate matrix is
$M$. Consequently, $M$ must have size
$\dim_k(X)\times\dim_k(Y)$.

The category implements the linear, additive, and abelian operations described
above. For example, we can compute Hom and endomorphism spaces, direct
sums, kernels, cokernels, and images through the common interface:

```@example linear_abelian_tour
using TensorCategories, Oscar

V = vector_spaces(QQ)
X = VectorSpaceObject(V, 2)
Y = VectorSpaceObject(V, 3)
f = morphism(X, Y, matrix(QQ, [1 0 0; 0 0 0]))

H = Hom(X, Y)
E = End(X)
@assert int_dim(H) == 6
@assert int_dim(E) == 4

D, inclusions, projections = direct_sum(X, Y)
@assert projections[1] ∘ inclusions[1] == id(X)
@assert projections[2] ∘ inclusions[2] == id(Y)
@assert is_zero(projections[1] ∘ inclusions[2])
@assert is_zero(projections[2] ∘ inclusions[1])

K, kernel_inclusion = kernel(f)
Q, cokernel_projection = cokernel(f)
I, image_inclusion = image(f)
@assert is_zero(f ∘ kernel_inclusion)
@assert is_zero(cokernel_projection ∘ f)
(int_dim(H), int_dim(E), int_dim(D), int_dim(K), int_dim(Q), int_dim(I))
```

The result is `(6,4,5,1,2,1)`: the dimensions of the Hom and endomorphism
spaces are the expected matrix dimensions, while the last three entries record
rank–nullity for the rank-one map $f$. The function `endomorphism_ring(X)`
turns `End(X)` into an explicit $k$-algebra when that representation is needed.

### [Finite-group representations](@id concrete-models)

Let $G$ be a finite group. An object of $\operatorname{Rep}_k(G)$ is a
finite-dimensional $k$-vector space $X$ together with a representation
$\rho_X:G\to\operatorname{GL}(X)$, and a morphism $f:X\to Y$ is a
$G$-equivariant linear map; see
[EGNO; Examples 2.3.4 and 2.10.13, pp. 26 and 43](@citet). This is a
$k$-linear abelian category over any field $k$. Its forgetful functor

```math
\label{eq:representation-forgetful-functor}
U:\operatorname{Rep}_k(G)\longrightarrow\operatorname{Vec}_k
```

is a matrix realization: it is faithful, though generally not full.
The implementation records dimensions, action matrices, and intertwiner
matrices; it does not store this forgetful functor as a separate Julia value.

TensorCategories.jl constructs this category with
`representation_category(k,G)`. An explicit representation can be given by
the images of group generators:

```julia
Representation(C, generators, matrices; check=true)
```

The keyword `check=true` verifies that the matrices satisfy the relations of
$G$. The implementation uses row coordinates. Thus an action matrix
$\rho_X(g)$ acts on the right of a row vector, and a matrix $M:X\to Y$ is an
intertwiner precisely when

```math
\label{eq:representation-intertwiner-foundations}
\rho_X(g)M=M\rho_Y(g)
```

for every generator $g$. The function `Hom(X,Y)` solves these simultaneous
linear equations, and `matrix(f)` returns the matrix $M$ of a represented
intertwiner.

For $G=C_3$ over $\mathbb Q$, the following matrix has order three and defines
a two-dimensional representation. We compare it with the trivial
one-dimensional representation and then use abelian operations on their direct
sum:

```@example linear_abelian_tour
G = cyclic_group(3)
C = representation_category(QQ, G)
A = matrix(QQ, [0 1; -1 -1])
X = Representation(C, gens(G), [A]; check=true)
T = Representation(C, gens(G), [identity_matrix(QQ, 1)]; check=true)

@assert A^3 == identity_matrix(QQ, 2)
@assert int_dim(Hom(T, X)) == 0
@assert int_dim(Hom(X, T)) == 0
@assert int_dim(End(T)) == 1
@assert int_dim(End(X)) == 2
@assert length(basis(End(X))) == 2
@assert all(size(matrix(f)) == (2, 2) for f in basis(End(X)))

D, inclusions, projections = direct_sum(T, X)
K, kernel_inclusion = kernel(projections[1])
Q, cokernel_projection = cokernel(inclusions[1])
@assert int_dim(K) == 2 && is_zero(projections[1] ∘ kernel_inclusion)
@assert int_dim(Q) == 2 && is_zero(cokernel_projection ∘ inclusions[1])
(int_dim(D), int_dim(K), int_dim(Q))
```

The two-dimensional endomorphism space of $X$ consists of matrices commuting
with $A$. The kernel and cokernel computations return representations, not
merely vector-space nullspaces: the implementation restricts or descends the
$G$-action to the computed subspace or quotient.

The coefficient field affects the intertwining equations and the resulting
abelian category. We therefore discuss coefficient fields before turning to
simple objects and composition factors.

Continue with [coefficient fields and numeric computations](@ref base-fields).
