# [Fusion and multifusion categories](@id tensor-conventions)

We now combine finiteness and semisimplicity with the monoidal and rigid
structures. Over an algebraically closed field, a **multifusion category** is a
finite semisimple multitensor category, and a **fusion category** is a finite
semisimple tensor category [EGNO; Definition 4.1.1](@cite). Here finite has the
meaning fixed for [finite abelian categories](@ref semisimple-categories).
Equivalently, a
fusion category has a simple tensor unit; over an algebraically closed field
Schur's lemma then identifies its endomorphism algebra with the coefficient
field.

No choice of simple representatives, bases, or matrix coordinates is part of
this definition. In particular, a fusion category need not be presented by
fusion rules or by structural matrices.

## Arbitrary coefficient fields

TensorCategories.jl also works over fields that are not algebraically closed.
It uses the terminology of [maurer2024computing; §2.1](@citet):

- a **weak multifusion category** is a finite semisimple $k$-linear rigid
  monoidal category, without requiring the tensor unit to be simple;
- a **weak fusion category** is a weak multifusion category whose tensor unit
  is simple;
- a **multifusion category** is a split weak multifusion category; and
- a **fusion category** is a split weak fusion category.

Here split means that every simple object $S$ is scalar:

```math
\label{eq:fusion-split-simple}
\operatorname{End}_{\mathcal C}(S)\cong k.
```

The adjective *weak* refers only to splitness over the chosen coefficient
field. It does not weaken rigidity, semisimplicity, or the monoidal coherence
axioms. Over an algebraically closed field the weak and split notions agree.

The package predicates are:

| Predicate | Implemented meaning |
|:---|:---|
| `is_weak_multifusion(C)` | finite, semisimple, rigid, with possibly decomposable unit and non-scalar simples |
| `is_weak_fusion(C)` | weak multifusion with simple unit |
| `is_multifusion(C)` | split weak multifusion |
| `is_fusion(C)` | split weak fusion |

Generic implications between these predicates are implemented in
`FrameworkChecks.jl`; individual category types can provide more direct
methods. The predicates report properties of the category over its current
field. Extending the field can change both the simple objects and which of the
four predicates applies.

## Representation categories

Finite-group representations show separately the effects of semisimplicity and
splitness. Let $G=C_3$. Over $\mathbb F_2$, Maschke's theorem gives a weak
fusion category, but the irreducible polynomial $x^2+x+1$ produces a
two-dimensional simple object with endomorphism field $\mathbb F_4$, so the
category is not split. Over $\mathbb F_3$, Maschke's condition fails and the
category is not weak fusion at all.

```@example characteristic_fusion
using TensorCategories, Oscar
G = cyclic_group(3)
C2 = representation_category(GF(2), G)
C3 = representation_category(GF(3), G)
@assert is_weak_fusion(C2) && !is_fusion(C2)
@assert !is_weak_fusion(C3)
(is_weak_fusion(C2), is_fusion(C2), is_weak_fusion(C3))
```

This example uses concrete representations and intertwiners throughout.

## Separability in positive characteristic

Over an imperfect field, some authors include separability in the definition
of a multifusion category. For a finite semisimple category this condition is
automatic over a perfect field, including number fields and finite fields, but
can be stronger over an imperfect field
[sanford2025fusion; Definition 2.9 and pp. 3--4](@cite). The package predicates
use the semisimple convention above and do not independently test this
separability condition.

Continue with [Grothendieck rings](@ref grothendieck-rings).
