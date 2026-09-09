# [Functors and natural transformations](@id linear-functors)

## Functors

All functor types in TensorCategories.jl are subtypes of
`TensorCategories.AbstractFunctor`. A functor must provide its domain and
codomain categories and its action on objects and morphisms. The general
constructor

```julia
F = functor(C, D, object_map, morphism_map)
```

returns a `Functor` from `C` to `D`. It can then be applied directly as `F(X)`
or `F(f)`. The functions `domain(F)` and `codomain(F)` return its source and
target categories.

The constructor stores the two supplied maps; it does not verify the functor
laws. In particular, the user must ensure that

```math
\label{eq:implemented-functor-laws}
F(\operatorname{id}_X)=\operatorname{id}_{F(X)},
\qquad
F(g\circ f)=F(g)\circ F(f).
```

Likewise, linearity, additivity, and exactness are properties of the supplied
maps rather than consequences of the constructor. A reusable category-specific
implementation can instead define a new subtype of
`TensorCategories.AbstractFunctor`, store the data it needs, and implement its
action on objects and morphisms by Julia call methods.

Functor composition follows the same order as morphism composition:
`compose(F,G)` and `G ∘ F` represent the functor $G\circ F$. The identity
functor of a category is `id(C)`.

### Example: Defining a functor

Here is the identity functor on finite-dimensional rational vector spaces,
entered through the general constructor.

```@example linear_functors
using TensorCategories, Oscar
C = vector_spaces(QQ)
F = functor(C, C, X -> X, f -> f)

X = VectorSpaceObject(C, 2)
f = morphism(X, X, matrix(QQ, [1 1; 0 1]))
@assert domain(F) == C && codomain(F) == C
@assert F(id(X)) == id(F(X))
@assert F(f ∘ f) == F(f) ∘ F(f)
F(f)
```

## Natural transformations

The abstract type `NaturalTransformation` is a subtype of `Morphism`. For
additive functors between Krull–Schmidt categories, the concrete type
`TensorCategories.AdditiveNaturalTransformation` represents a natural
transformation by its components on a chosen family of indecomposable
objects.

Suppose $F,G\colon\mathcal C\to\mathcal D$ are additive functors and
$S_1,\ldots,S_r$ are indecomposable representatives from which the objects
under consideration can be decomposed. A transformation is entered as

```julia
eta = TensorCategories.AdditiveNaturalTransformation(
    F, G, indecomposables, components
)
```

where `indecomposables` is the list $[S_1,\ldots,S_r]$ and `components`
specifies morphisms $\eta_{S_i}\colon F(S_i)\to G(S_i)$. Components may be
given in the same order as the objects, or as pairs `S => eta_S`; omitted
pairs are interpreted as zero components.

For a decomposition $X\cong\bigoplus_j S_j$ with inclusions
$i_j\colon S_j\to X$ and projections $p_j\colon X\to S_j$, the implementation
reconstructs the component on $X$ as

```math
\label{eq:additive-natural-transformation-extension}
\eta_X
=
\sum_j G(i_j)\circ\eta_{S_j}\circ F(p_j).
```

Thus `eta(X)` is available for decomposable objects as soon as the category can
compute `direct_sum_decomposition(X, indecomposables)`.

The generic `direct_sum_decomposition` method uses decomposition into simple
objects and therefore applies to semisimple categories. A nonsemisimple
Krull–Schmidt model can use the same representation of natural transformations
after providing its own `direct_sum_decomposition` method.

The constructor does not check naturality. The supplied components define a
natural transformation only if

```math
\label{eq:natural-transformation-naturality-linear}
G(f)\circ\eta_S=\eta_T\circ F(f)
\qquad
\text{for every }f\colon S\to T
```

between the chosen indecomposable objects. This includes all endomorphisms of
each indecomposable, not only morphisms between distinct objects.

The function `Nat(F,G; indecomposables=objects)` computes the vector space of
additive natural transformations by solving these linear equations. If the
keyword is omitted, it calls `indecomposables(domain(F))`; the generic fallback
can enumerate indecomposables only in the semisimple case. A nonsemisimple
Krull–Schmidt model must provide its own method or pass a suitable list
explicitly. The solver requires additive functors with the same domain and
codomain, finite bases for the relevant Hom spaces, and a common coefficient
field. Evaluating the resulting transformations away from the supplied
representatives additionally requires effective direct-sum decompositions. The
solver returns a `NaturalTransformations` Hom space; its basis elements can be
obtained with `basis`.

### Example: Components on indecomposable objects

Over $\mathbb F_3$, the category $\operatorname{Rep}_{\mathbb F_3}(C_2)$ has two
simple, hence indecomposable, objects. Since there are no morphisms between
the two simple objects, a natural endomorphism of the identity functor may act
on them by independent scalars.

```@example additive_natural_transformations
using TensorCategories, Oscar
C = representation_category(GF(3), cyclic_group(2))
S = simples(C)
F = id(C)

eta = TensorCategories.AdditiveNaturalTransformation(
    F, F, S, [S[1] => id(S[1]), S[2] => -id(S[2])]
)

X, inclusions, projections = direct_sum(S[1], S[2])
@assert eta(X) ∘ inclusions[1] == inclusions[1] ∘ eta(S[1])
@assert eta(X) ∘ inclusions[2] == inclusions[2] ∘ eta(S[2])

N = Nat(F, F; indecomposables=S)
@assert int_dim(N) == 2
int_dim(N)
```

The two-dimensional result records the independent scalar components on the
two simple objects. In a nonsemisimple Krull–Schmidt category, morphisms
between distinct indecomposables and non-scalar endomorphisms impose further
equations through \eqref{eq:natural-transformation-naturality-linear}.

Continue with [tensor products and associators](@ref monoidal-categories).
