# [Idempotents and Krull–Schmidt categories](@id karoubian-categories)

Let $\mathcal C$ be an additive category. An endomorphism
$e\colon X\to X$ is idempotent if $e^2=e$. It **splits** if there are an
object $Y$ and morphisms

```math
Y\xrightarrow{i}X\xrightarrow{p}Y
```

such that

```math
\label{eq:idempotent-splitting}
p\circ i=\operatorname{id}_Y,
\qquad
i\circ p=e.
```

An additive category is **Karoubian**, or **idempotent complete**, if every
idempotent splits. In an abelian category this is automatic: one may take
$Y=\operatorname{im}(e)$. For an arbitrary additive category, its Karoubi
envelope formally adjoins objects $(X,e)$ for idempotents
$e\in\operatorname{End}(X)$; a morphism $(X,e)\to(Y,d)$ is a morphism
$f\colon X\to Y$ satisfying $f=d\circ f\circ e$.

TensorCategories.jl does not currently provide a general Karoubi-envelope type
for every category. Abelian implementations split idempotents through their
image algorithms. The public function `karoubian_envelope` is currently
specialized to center and centralizer categories, where this completion is
needed by the corresponding construction.

For example, the projection onto the first coordinate of
$\mathbb Q^2$ splits through its one-dimensional image:

```@example split_idempotent
using TensorCategories, Oscar
C = vector_spaces(QQ)
X = VectorSpaceObject(C, 2)
e = morphism(X, X, QQ[1 0; 0 0])
Y, i = image(e)
p = left_inverse(i) ∘ e
@assert int_dim(Y) == 1
@assert p ∘ i == id(Y)
@assert i ∘ p == e
Y
```

Here `image(e)` returns the image object $Y$ and its inclusion
$i\colon Y\to X$. Restricting $e$ to this image gives the projection
$p\colon X\to Y$, and the assertions verify both equations in
\eqref{eq:idempotent-splitting}. The image alone is enough for the generic
decomposition algorithms described below.

## Krull–Schmidt decomposition

An additive category is a **Krull–Schmidt category** if every object is a finite
direct sum of objects having local endomorphism rings. These summands are
indecomposable, and the resulting decomposition is unique up to permutation and
isomorphism. Equivalently, an additive category is Krull–Schmidt when it has
split idempotents and every endomorphism ring is semiperfect
[krause2015krull; Corollary 4.4](@cite). In particular, a Hom-finite
$k$-linear additive category is Krull–Schmidt precisely when it is
idempotent complete.

The distinction between *indecomposable* and *simple* matters outside a
semisimple category. An indecomposable object has a local endomorphism ring,
while Schur's lemma only says that a simple object has a division endomorphism
ring. A representation in modular characteristic can be indecomposable without
being simple.

The package uses:

| Operation | Meaning |
|:---|:---|
| `is_indecomposable(X)` | test whether $X$ is indecomposable |
| `decompose(X)` | return pairs `(Y,m)` of indecomposable summands and multiplicities |
| `TensorCategories.is_krull_schmidt(C)` | report that the implementation treats $\mathcal C$ as Krull–Schmidt |

The generic decomposition algorithm is not uniform over all fields. Over a
finite field it decomposes the regular right module of
$\operatorname{End}(X)$, obtains primitive idempotents, and forms their
images. In a semisimple category over supported exact fields it instead uses
the semisimple endomorphism algebra. Other nonsemisimple coefficient fields
require a category-specific backend. Central primitive idempotents only split
blocks; they need not give the individual indecomposable summands.

Continue with [finite and semisimple categories](@ref semisimple-categories).
