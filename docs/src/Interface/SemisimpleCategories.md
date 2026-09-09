# [Finite and semisimple categories](@id semisimple-categories)

An abelian category $\mathcal C$ is **semisimple** if every object is a direct
sum of simple objects [EGNO; Definition 1.5.1](@cite). If $\mathcal C$ is
locally finite, these direct sums are finite. It is then a Krull--Schmidt
category, and its simple objects are precisely its indecomposable objects.
Thus every object $X$ admits a decomposition

```math
\label{eq:semisimple-decomposition}
X\cong\bigoplus_{S\in\operatorname{Irr}(\mathcal C)}S^{\oplus m_S}.
```

If $D_S=\operatorname{End}_{\mathcal C}(S)$, then
$\operatorname{Hom}_{\mathcal C}(S,X)$ is a right $D_S$-module by
precomposition and

```math
\label{eq:nonsplit-hom-multiplicity}
m_S=\dim_{D_S}\operatorname{Hom}_{\mathcal C}(S,X).
```

Schur's lemma says that $D_S$ is a division algebra
[EGNO; Lemma 1.5.2](@cite). Outside a semisimple category, the converse in
Schur's lemma fails: an object with division endomorphism algebra need not be
simple, and an indecomposable object need not be simple.

## The interface

The predicate `is_semisimple(C)` records that this property is known for the
implemented category. When the backend can enumerate the simple objects,
`simples(C)` returns chosen representatives of their isomorphism classes. The
order of this list belongs to the implementation and must not be treated as a
mathematical invariant.

For an object $X$, `decompose(X)` returns pairs `(S,m)` representing the
decomposition in equation \eqref{eq:semisimple-decomposition}. The generic
implementation of `decompose(X,S)` computes the multiplicity using equation
\eqref{eq:nonsplit-hom-multiplicity}; it does not replace the $D_S$-dimension
by the dimension over $k$.

The generic `is_simple(X)` uses indecomposability only when the parent category
is semisimple. Otherwise a category-specific test is required.

## Example: Finite-group representations

For a finite group $G$, Maschke's theorem gives

```math
\label{eq:maschke-condition}
\operatorname{Rep}_k(G)\text{ is semisimple}
\quad\Longleftrightarrow\quad
\operatorname{char}(k)\nmid |G|.
```

```@example representation_semisimplicity
using TensorCategories, Oscar
G = cyclic_group(5)
C0 = representation_category(QQ, G)
C5 = representation_category(GF(5), G)
@assert is_semisimple(C0)
@assert !is_semisimple(C5)
(is_semisimple(C0), is_semisimple(C5))
```

Equation \eqref{eq:maschke-condition} concerns semisimplicity only. Even when it
holds, the simple representations need not remain simple after extending the
coefficient field.

When a semisimple category has only finitely many simple isomorphism classes,
it is finite: every object is projective, so the direct sum of representatives
of the simples is a projective generator.

Continue with [scalar extension and splitting](@ref splitting-and-scalars).
