# [Ring and tensor categories](@id ring-and-tensor-categories)

We now combine the linear, abelian, monoidal, and rigid structures. Following
[EGNO; Definition 4.1.1](@citet), a **multitensor category** over an
algebraically closed field $k$ is a locally finite $k$-linear abelian rigid
monoidal category whose tensor product is bilinear on morphisms. It is a
**tensor category** when

```math
\label{eq:tensor-category-unit}
\operatorname{End}_{\mathcal C}(\mathbb 1)\cong k.
```

Because tensoring with an object in a rigid category has both adjoints, it is
exact. Dropping rigidity gives the more general terminology of
[EGNO; Definition 4.2.3](@citet): a **multiring category** is a locally finite
$k$-linear abelian monoidal category whose tensor product is bilinear and
exact in each variable, and it is a **ring category** when equation
\eqref{eq:tensor-category-unit} holds. Thus every multitensor category is a
multiring category, and every tensor category is a ring category.

A ring category is a category with specified linear and monoidal properties; it
should not be confused with the Grothendieck ring constructed from such a
category later in this chapter.

Over a general coefficient field, TensorCategories.jl also supports the
non-split convention of [maurer2024computing; §2.1](@citet), in which a simple
tensor unit need not have endomorphism algebra $k$. Algorithms which require
the split condition must impose it separately.

## The interface

The corresponding package predicates are `is_multiring`, `is_ring`,
`is_multitensor`, and `is_tensor`. Stronger structural declarations imply the
weaker ones through generic fallbacks. In particular, every multiring category
reports `is_locally_finite(C) == true`. The predicates report what a backend
has declared or established; they do not reconstruct the axioms from the
available methods.

In the non-split convention, a category with a simple but non-scalar tensor
unit can be reported as
`is_tensor(C) == true` and `is_ring(C) == true`, even though equation
\eqref{eq:tensor-category-unit} fails. Algorithms that require the split EGNO
condition must additionally test `is_split_semisimple(C)`. This broader
behavior is part of the current package interface.

## Example: Group representations

Both $\operatorname{Vec}_k$ and $\operatorname{Rep}_k(G)$ are tensor
categories in the sense above, whether or not $\operatorname{Rep}_k(G)$ is
semisimple. The representation category therefore illustrates why *tensor
category* and *fusion category* are not synonyms.

```@example tensor_predicates
using TensorCategories, Oscar
G = cyclic_group(3)
C = representation_category(GF(3), G)
@assert is_tensor(C) && is_ring(C)
@assert !is_semisimple(C)
(is_tensor(C), is_semisimple(C))
```

Continue with [fusion and multifusion categories](@ref tensor-conventions).
