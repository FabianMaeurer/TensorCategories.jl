# [Pivotal and spherical structures](@id pivotal-braided)

Fix the chosen left duality from the preceding pages. A **pivotal structure** is
a monoidal natural isomorphism

```math
\label{eq:pivotal-component}
j_X\colon X\longrightarrow X^{**}.
```

Changing the pivotal structure can change right traces and dimensions without
changing the underlying left dual objects.

## Traces and dimensions

For an endomorphism $f\colon X\to X$, the chosen duality and pivotal
structure define left and right pivotal traces. These are endomorphisms of the
tensor unit. Converting them to scalars therefore requires the relevant
endomorphism to be a scalar multiple of $\operatorname{id}_{\mathbb 1}$.

The package convention is

```math
\label{eq:package-dimension-and-norm}
\dim(X)=\dim_L(X),
\qquad
|X|^2=\dim(X)\dim(X^*).
```

For a multifusion category with simple representatives $S_i$, the category
dimension is

```math
\label{eq:package-category-dimension}
\dim(\mathcal C)=\sum_i |S_i|^2.
```

This is distinct from the Frobenius–Perron dimension, which depends only on the
Grothendieck ring.

A pivotal structure is **spherical** when its left and right traces agree. In
the split semisimple setting it is enough to compare the left and right
dimensions on simple objects [EGNO; Definition 4.7.14 and Theorem 4.7.15](@cite).

## The interface

| Operation | Meaning |
|:---|:---|
| `pivotal(X)` | the component $j_X\colon X\to X^{**}$ |
| `left_trace(f)`, `right_trace(f)` | the two pivotal traces of an endomorphism |
| `tr(f)` | the left pivotal trace |
| `dim(X)` | the left pivotal dimension of $X$ |
| `TensorCategories.squared_norm(X)` | the squared norm $|X|^2$ |
| `dim(C)` | the category dimension |
| `is_pivotal(C)` | report that $\mathcal C$ has a pivotal structure |
| `is_spherical(C)` | report that the pivotal structure is spherical |

The generic right-duality methods transport the chosen left duality through
`pivotal(X)`. In supported finite semisimple models,
`is_pivotal(C; check=true)` checks invertibility and monoidality on the chosen
simple representatives. It does not infer naturality merely from a list of
components; a backend with non-scalar simple endomorphisms must account for
that naturality itself.

The generic `is_spherical(C; check=true)` first checks the pivotal structure and
then performs this comparison. For a non-split model, a category-specific
method is required; the generic checked predicate does not claim sphericality.

## Example: Vector spaces

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
