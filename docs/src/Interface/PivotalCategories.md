# [Pivotal and spherical structures](@id pivotal-braided)

Fix the chosen left duality from the preceding pages. A **pivotal structure** is
a monoidal natural isomorphism

```math
\label{eq:pivotal-component}
j_X\colon X\longrightarrow X^{**}.
```

The method `pivotal(X)` returns the component in equation
\eqref{eq:pivotal-component}. The generic right-duality methods transport the
chosen left duality through this isomorphism. Thus changing the pivotal
structure can change right traces and dimensions without changing the
underlying left dual objects.

The predicate `is_pivotal(C)` reports a stored or category-specific pivotal
structure. In supported finite semisimple models, `is_pivotal(C; check=true)`
checks invertibility and monoidality on the chosen simple representatives. It
does not infer naturality merely from a list of components; a backend with
non-scalar simple endomorphisms must account for that naturality itself.

## Traces and dimensions

For an endomorphism $f\colon X\to X$, the chosen duality and pivotal
structure define left and right pivotal traces. The package calls them
`left_trace(f)` and `right_trace(f)`; `tr(f)` abbreviates the left trace. These
are endomorphisms of the tensor unit. Converting them to scalars therefore
requires the relevant endomorphism to be a scalar multiple of
$\operatorname{id}_{\mathbb 1}$.

The package convention is

```math
\label{eq:package-dimension-and-norm}
\dim(X)=\dim_L(X),
\qquad
|X|^2=\dim(X)\dim(X^*).
```

The corresponding calls are `dim(X)` and
`TensorCategories.squared_norm(X)`. For a multifusion category with simple
representatives $S_i$, the generic category dimension is

```math
\label{eq:package-category-dimension}
\dim(\mathcal C)=\sum_i |S_i|^2.
```

This is distinct from the Frobenius–Perron dimension, which depends only on the
Grothendieck ring.

A pivotal structure is **spherical** when its left and right traces agree. In
the split semisimple setting it is enough to compare the left and right
dimensions on simple objects [EGNO; Definition 4.7.14 and Theorem 4.7.15](@cite).
The generic `is_spherical(C; check=true)` first checks the pivotal structure and
then performs this comparison. For a non-split model, a category-specific
method is required; the generic checked predicate does not claim sphericality.

```@example pivotal_vector_spaces
using TensorCategories, Oscar
C = vector_spaces(QQ)
X = VectorSpaceObject(C, 3)
@assert is_pivotal(C; check=true)
@assert is_spherical(C; check=true)
@assert dim(X) == 3
(dim(X), dim(dual(X)))
```

Continue with [braided and symmetric categories](@ref braided-categories).
