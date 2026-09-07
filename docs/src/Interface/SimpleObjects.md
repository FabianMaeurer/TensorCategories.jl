# [Simple objects and finite-length categories](@id simple-objects)

Let $\mathcal C$ be an abelian category. A nonzero object $S$ is *simple* if
its only subobjects are $0$ and $S$. An object $X$ has *finite length* if it
admits a composition series

```math
\label{eq:composition-series}
0=X_0\subset X_1\subset\cdots\subset X_n=X
```

whose successive quotients $X_i/X_{i-1}$ are simple. The Jordan--Hölder
theorem states that their isomorphism classes and multiplicities do not depend
on the chosen composition series; see [EGNO; §1.5](@cite).

Simple composition factors are subquotients. They need not be subobjects or
direct summands. An indecomposable object is likewise not necessarily simple:
a nonsplit extension of two simple objects is indecomposable but has length
two. These distinctions matter especially in positive characteristic.

## The interface

The relevant public functions have different computational requirements:

| Function | Meaning |
|:---|:---|
| `is_simple(X)` | decide whether the nonzero object $X$ has no proper nonzero subobject |
| `composition_factors(X)` | return pairs `(S,m)` recording the simple factors and their Jordan--Hölder multiplicities |
| `simple_subobjects(X)` | return the simple isomorphism types occurring in the socle, when supported |
| `simples(C)` | enumerate representatives of all simple isomorphism classes, when this is finite and computable |

There is no general algorithm for `composition_factors(X)` from the bare
abelian interface. A category of finite length must supply a
category-specific method if it supports this computation. In particular, the
endomorphism algebra $\operatorname{End}(X)$ alone does not determine a
composition series. The generic `is_simple` method applies to supported
semisimple categories; a nonsemisimple model must provide its own valid test.

Enumeration by `simples(C)` is a further problem. Knowing how to test one
given object for simplicity does not provide an algorithm that finds every
simple object of a category.

## Example: Modular group representations

For a finite group $G$ and a field $k$, finite-dimensional
$k$-representations form a finite-length abelian category. When the
characteristic of $k$ divides $|G|$, this category need not be semisimple.
TensorCategories.jl converts representations over finite fields to GAP modules.
It uses GAP's MeatAxe routines `MTX.IsIrreducible` and
`MTX.CollectedFactors` to test simplicity and compute composition factors; see
[gapmanual2026; §§69.5 and 69.7](@cite). These algorithms use the action
matrices of the module, not merely its endomorphism algebra.

Let $G=C_5$ and $k=\mathbb F_5$. The matrix

```math
\label{eq:modular-jordan-block}
J=\begin{pmatrix}1&1\\0&1\end{pmatrix}
```

satisfies $J^5=1$ in characteristic $5$ and therefore defines a
two-dimensional representation of $C_5$. It is a nonsplit self-extension of
the trivial representation.

```@example modular_composition_factors
using TensorCategories, Oscar
F = GF(5)
G = cyclic_group(5)
C = representation_category(F, G)
J = matrix(F, [1 1; 0 1])
X = Representation(C, gens(G), [J]; check=true)

@assert J^5 == identity_matrix(F, 2)
@assert !is_simple(X)
factors = composition_factors(X)
@assert length(factors) == 1
S, multiplicity = only(factors)
@assert is_simple(S) && multiplicity == 2
(int_dim(X), int_dim(S), multiplicity)
```

Although both composition factors are trivial, $X$ is not their direct sum.
Its simple subobject belongs to the socle, and its simple quotient belongs to
the top. Neither splits off as a direct summand:

```@example modular_composition_factors
socle_types = simple_subobjects(X)
@assert length(socle_types) == 1
S = only(socle_types)
inclusions = basis(Hom(S, X))
projections = basis(Hom(X, S))
@assert length(inclusions) == 1 && length(projections) == 1
@assert is_zero(only(projections) ∘ only(inclusions))
nothing # hide
```

Direct-sum decompositions into indecomposable objects are governed by the
Krull--Schmidt property rather than by composition series. We turn to that
next.

Continue with [idempotents and Krull--Schmidt categories](@ref
karoubian-categories).
