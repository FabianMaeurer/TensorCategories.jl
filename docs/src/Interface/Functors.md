# [Functors and natural transformations](@id linear-functors)

A functor must specify its action on objects and morphisms. In a linear or
abelian setting, additional adjectives impose familiar compatibility:

| Functor | Additional requirement |
|:---|:---|
| additive | preserves finite direct sums |
| $k$-linear | acts $k$-linearly on every Hom space |
| exact | preserves short exact sequences |

The generic constructor `functor(C,D,obj_map,mor_map)` stores an object map and
a morphism map. It does not infer or check the functor laws or any of the
additional properties in the table.

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

A custom subtype of `TensorCategories.AbstractFunctor` can instead implement
call methods on objects and morphisms. Functor composition has the same order
as morphism composition: `compose(F,G)` represents $G\circ F$.

## Natural transformations

For functors $F,G\colon\mathcal C\to\mathcal D$, a natural transformation
$\eta\colon F\Rightarrow G$ has components
$\eta_X\colon F(X)\to G(X)$ satisfying

```math
\label{eq:natural-transformation-naturality-linear}
G(f)\circ\eta_X=\eta_Y\circ F(f)
\qquad
\text{for every }f\colon X\to Y.
```

An `AdditiveNaturalTransformation` stores components on specified
indecomposable objects and extends them using direct-sum decompositions.
`Nat(F,G)` solves equation
\eqref{eq:natural-transformation-naturality-linear} in supported finite
additive models. The solver requires finite Hom bases, effective
Krull–Schmidt decompositions and coordinates, additive functors, and a common
base field. Semisimplicity is not a formal requirement.

If $S$ is an indecomposable with a nontrivial endomorphism algebra, its
component must satisfy naturality with every $d\in\operatorname{End}(S)$:

```math
\label{eq:natural-transformation-endomorphisms}
G(d)\circ\eta_S=\eta_S\circ F(d).
```

Choosing an arbitrary matrix on each indecomposable therefore does not in
general define a natural transformation.

Continue with [tensor products and associators](@ref monoidal-categories).
