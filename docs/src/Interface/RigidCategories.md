# [Rigidity](@id rigid-categories)

A chosen left dual of an object $X$ consists of an object $X^*$ and
morphisms

```math
\label{eq:left-duality-morphisms}
\operatorname{ev}_X\colon X^*\otimes X\longrightarrow\mathbb 1,
\qquad
\operatorname{coev}_X\colon\mathbb 1\longrightarrow X\otimes X^*.
```

This is the convention of [EGNO; Definition 2.10.1](@citet).

With the associator direction in equation \eqref{eq:monoidal-associator}, the
two triangle identities are

```math
\label{eq:left-duality-triangle-x}
(\operatorname{id}_X\otimes\operatorname{ev}_X)
\circ a_{X,X^*,X}
\circ(\operatorname{coev}_X\otimes\operatorname{id}_X)
=\operatorname{id}_X.
```

and

```math
\label{eq:left-duality-triangle-dual}
(\operatorname{ev}_X\otimes\operatorname{id}_{X^*})
\circ a^{-1}_{X^*,X,X^*}
\circ(\operatorname{id}_{X^*}\otimes\operatorname{coev}_X)
=\operatorname{id}_{X^*}.
```

A chosen right dual ${}^*X$ has morphisms

```math
\label{eq:right-duality-morphisms}
\widetilde{\operatorname{ev}}_X\colon
X\otimes{}^*X\longrightarrow\mathbb 1,
\qquad
\widetilde{\operatorname{coev}}_X\colon
\mathbb 1\longrightarrow{}^*X\otimes X.
```

This agrees with [EGNO; Definition 2.10.2](@citet).

A monoidal category is **rigid** if every object has left and right duals.

## The interface

| Operation | Meaning |
|:---|:---|
| `dual(X)`, `left_dual(X)` | the chosen left dual $X^*$ |
| `ev(X)` | the left evaluation $X^*\otimes X\to\mathbb 1$ |
| `coev(X)` | the left coevaluation $\mathbb 1\to X\otimes X^*$ |
| `right_dual(X)` | the chosen right dual ${}^*X$ |
| `right_ev(X)` | the right evaluation $X\otimes{}^*X\to\mathbb 1$ |
| `right_coev(X)` | the right coevaluation $\mathbb 1\to{}^*X\otimes X$ |
| `TensorCategories.is_rigid(C)` | report that every object has left and right duals |

The dual objects alone are insufficient: evaluation, coevaluation, and their
normalizations are part of the implemented data. A chosen left dual is not
automatically a chosen right dual at the level of the interface. A model must
implement the right-duality data or provide additional structure from which
the generic methods can construct it. The predicate `is_rigid(C)` is a
structural declaration; it does not check the triangle identities for
arbitrary objects.

## Example: Vector spaces and representations

For vector spaces and group representations, `dual(X)` uses the dual vector
space. In row coordinates, the action on the dual representation is
$\rho_{X^*}(g)=\rho_X(g^{-1})^{\mathsf T}$. Evaluation and coevaluation are
the usual contraction and coevaluation tensors in the chosen bases.

```@example concrete_rigidity
using TensorCategories, Oscar
C = vector_spaces(QQ)
X = VectorSpaceObject(C, 2)
triangle = (id(X) ⊗ ev(X)) ∘ associator(X,dual(X),X) ∘
    (coev(X) ⊗ id(X))
@assert triangle == id(X)
triangle == id(X)
```

Continue with [ring and tensor categories](@ref ring-and-tensor-categories).
