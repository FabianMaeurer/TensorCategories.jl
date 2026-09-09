# [Monoidal functors](@id monoidal-functors)

Let $F\colon\mathcal C\to\mathcal D$ be a functor between monoidal
categories. A strong monoidal structure on $F$ includes natural
isomorphisms

```math
\label{eq:monoidal-functor-tensorator}
J_{X,Y}\colon F(X)\otimes F(Y)
\longrightarrow F(X\otimes Y)
```

and a compatible unit isomorphism, satisfying the associator and unit
coherence equations [EGNO; §2.4](@cite). The direction in equation
\eqref{eq:monoidal-functor-tensorator} is the direction used by
TensorCategories.jl.

A **tensor functor** in the terminology of [EGNO; Definition 4.2.5](@citet) is
an exact faithful $k$-linear monoidal functor between multiring categories.
Some literature uses this term for a strong monoidal functor without the
exactness and faithfulness hypotheses; the manual uses the EGNO meaning.

Monoidal structure is additional data. A functor that preserves tensor-product
classes in a Grothendieck ring does not thereby acquire tensorators, and an
objectwise identification $F(X\otimes Y)\cong F(X)\otimes F(Y)$ does not
supply a coherent choice of them.

## The interface

The implemented solvers assume additive $k$-linear behavior, strict
preservation of the represented unit, and normalized unit tensorators. The
constructor checks the source and target and the image of the unit, but it does
not prove that the supplied functor is additive or linear.

`monoidal_structure_candidates(F; check=false)` searches for tensorators in
supported split fusion categories. The solution scheme can have positive
dimension, in which case sampling candidates is not a classification and an
empty sample is not a proof of nonexistence. The option `check=true`
re-evaluates the monoidal-functor equation on each returned candidate.
`monoidal_structures(F)` currently restricts its generic complete solver to the
normalized case in which the source has one simple object.

For monoidal functors $F,G$, a monoidal natural transformation is an ordinary
natural transformation compatible with the two tensorators. The method
`monoidal_natural_transformations(F,G)` solves these additional equations for
the supported finite additive models.

Continue with [fiber functors](@ref fiber-functors).
