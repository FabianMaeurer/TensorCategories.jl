# [Braided and symmetric categories](@id braided-categories)

A braiding in TensorCategories.jl has direction

```math
\label{eq:package-braiding}
c_{X,Y}\colon X\otimes Y\longrightarrow Y\otimes X.
```

This agrees with the convention of [EGNO; Definition 8.1.1](@cite).
Naturality and the two hexagon equations use the associator direction fixed in
equation \eqref{eq:monoidal-associator}.

Explicitly, the two hexagon equations are

```math
\label{eq:package-positive-hexagon}
a_{Y,Z,X}\circ c_{X,Y\otimes Z}\circ a_{X,Y,Z}
=
(\operatorname{id}_Y\otimes c_{X,Z})
\circ a_{Y,X,Z}
\circ(c_{X,Y}\otimes\operatorname{id}_Z)
```

and

```math
\label{eq:package-negative-hexagon}
a^{-1}_{Z,X,Y}\circ c_{X\otimes Y,Z}\circ a^{-1}_{X,Y,Z}
=
(c_{X,Z}\otimes\operatorname{id}_Y)
\circ a^{-1}_{X,Z,Y}
\circ(\operatorname{id}_X\otimes c_{Y,Z}).
```

The first has source $(X\otimes Y)\otimes Z$, while the second has source
$X\otimes(Y\otimes Z)$. These source bracketings distinguish the two
equations even in a concrete model whose represented associators are identity
matrices.

A braided category is **symmetric** when

```math
\label{eq:symmetric-braiding}
c_{Y,X}\circ c_{X,Y}=\operatorname{id}_{X\otimes Y}
```

for all objects $X,Y$. The standard braiding on vector spaces is the flip,
and the same flip is equivariant for the diagonal action on
$\operatorname{Rep}_k(G)$ in every characteristic.

Given a braiding and pivotal structure, the package uses the twist convention

```math
\label{eq:package-twist}
\theta_X=u_X^{-1}\circ j_X,
```

where $u_X\colon X\to X^{**}$ is the Drinfeld isomorphism and
$j_X\colon X\to X^{**}$ is the pivotal component. This is the convention of
[EGNO; §8.10](@cite).

## The interface

| Operation | Meaning |
|:---|:---|
| `braiding(X,Y)` | the braiding $c_{X,Y}\colon X\otimes Y\to Y\otimes X$ |
| `is_braided(C)` | report that $\mathcal C$ supplies a braiding |
| `hexagon_axiom(C)` | check both hexagon equations on all listed simple triples |
| `twist(X)` | the twist $\theta_X$ |
| `twist_scalar(X)` | return the scalar of $\theta_X$ when it is a scalar endomorphism |

Providing `braiding(X,Y)` does not establish naturality or the hexagon axioms.
In supported finite semisimple models, `hexagon_axiom(C)` performs the
exhaustive check on listed simple objects. A successful ball-valued check means
that the equations hold at the chosen working precision. The function
`twist_scalar(X)` requires the twist to be represented by a unique scalar
multiple of the identity.

## Example: Group representations

The implemented categories $\operatorname{Rep}_k(G)$ are symmetric.

```@example representation_braiding
using TensorCategories, Oscar
G = cyclic_group(3)
C = representation_category(QQ, G)
X = one(C)
c = braiding(X,X)
@assert c ∘ c == id(X ⊗ X)
matrix(c)
show(stdout, MIME"text/plain"(), matrix(c)); println() # hide
```

Continue with [unitarity and modularity](@ref unitary-modular-categories).
