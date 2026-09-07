# [Braided and symmetric categories](@id braided-categories)

A braiding in TensorCategories.jl has direction

```math
\label{eq:package-braiding}
c_{X,Y}\colon X\otimes Y\longrightarrow Y\otimes X.
```

The call `braiding(X,Y)` returns the morphism in equation
\eqref{eq:package-braiding}. This agrees with the convention of
[EGNO; Definition 8.1.1](@cite). Naturality and the two hexagon equations use
the associator direction fixed in equation \eqref{eq:monoidal-associator}.

Explicitly, the two equations checked by `hexagon_axiom` are

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
$\operatorname{Rep}_k(G)$ in every characteristic. These implemented
representation categories are therefore symmetric.

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

The predicate `is_braided(C)` records that the category supplies a braiding.
It does not establish naturality or the hexagon axioms from the existence of a
method. In supported finite semisimple models, `hexagon_axiom(C)` checks all
triples of listed simple objects. A successful ball-valued check means that the
two sides agree at the chosen working precision.

Given a braiding and pivotal structure, the package uses the twist convention

```math
\label{eq:package-twist}
\theta_X=u_X^{-1}\circ j_X,
```

where $u_X\colon X\to X^{**}$ is the Drinfeld isomorphism and
$j_X\colon X\to X^{**}$ is the pivotal component. This is the convention of
[EGNO; §8.10](@cite). The call `twist(X)` returns the resulting endomorphism;
`twist_scalar(X)` additionally requires it to be represented by a unique scalar
multiple of the identity.

Continue with [unitarity and modularity](@ref unitary-modular-categories).
