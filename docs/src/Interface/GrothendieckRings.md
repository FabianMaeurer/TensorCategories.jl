# [Grothendieck rings](@id grothendieck-rings)

The Grothendieck ring is a decategorification of a multiring category: objects are
replaced by their classes, short exact sequences give additive relations, and
tensor products give multiplication. In a semisimple category the additive
relations are precisely the direct-sum relations. For a fusion category, the
multiplication table is exactly the table of fusion rules.

[EGNO; Chapter 3](@citet) develops the resulting combinatorics as the theory of
$\mathbb Z_+$-rings; §§4.5 and 4.9 explain how these rings arise from tensor
and fusion categories. TensorCategories.jl implements both passage from a
category to its ring and direct computation with the resulting based ring.

## Classes of objects

For an abelian category of finite-length objects, the Grothendieck group
$\operatorname{Gr}(\mathcal C)$ has a generator $[X]$ for each isomorphism
class, with relations
```math
\label{eq:grothendieck-relation}
[Y]=[X]+[Z]
\quad\text{for every exact sequence }0\longrightarrow X\longrightarrow Y
\longrightarrow Z\longrightarrow0.
```
Its basis consists of the classes of simple objects, and
$[X]=\sum_i[X:S_i][S_i]$, where $[X:S_i]$ is a composition multiplicity
[EGNO; Definition 1.5.8, p. 5](@cite).

When tensor product is exact in both variables, it induces multiplication
$[X][Y]=[X\otimes Y]$, with identity $[\mathbb 1]$. Associativity follows from the
associator, since isomorphic objects have equal classes
[EGNO; §4.5, pp. 71--72](@cite).

The **split Grothendieck ring** instead imposes only the direct-sum relations
$[X\oplus Y]=[X]+[Y]$. In a Krull–Schmidt category its basis consists of classes
of indecomposable objects. For semisimple categories the two constructions
agree, and the coefficients of $[X]$ are its simple-summand multiplicities.

The function `split_grothendieck_ring(C)` implements the direct-sum
construction in models with finite indecomposable enumeration and tensor-product
decomposition. Here “split” refers to the relations, not to a splitting field:
the function also applies to supported non-split semisimple categories.
For nonsemisimple input it does not impose relations from nonsplit exact
sequences.

## Positive bases and fusion rings

Write $b_i=[S_i]$. In the semisimple case the fusion rules are
```math
\label{eq:fusion-ring-product}
b_i b_j=\sum_l N_{ij}^{\,l}b_l,\qquad
N_{ij}^{\,l}=[S_i\otimes S_j:S_l]\in\mathbb Z_{\geq0}.
```
Thus $\operatorname{Gr}(\mathcal C)$ is a ring over $\mathbb Z$ with a
distinguished basis whose structure constants are nonnegative integers. A
**$\mathbb Z_+$-ring** has such a basis and an identity whose coordinates in
that basis are nonnegative; it is called *unital* when the identity itself is
a basis element [EGNO; Definition 3.1.1, p. 49](@cite).

The positive cone $\sum_i\mathbb Z_{\geq0}b_i$ records actual objects of a
semisimple category.
The ring also contains virtual classes with negative coefficients. Its
coefficient ring is $\mathbb Z$ even when the category is defined over a field of
positive characteristic.

The package type is `ZPlusRing`; `ℤ₊Ring` and `ℕRing` are aliases for the same
type. The last alias reflects nonnegative structure constants, but the
mathematical terminology used in this manual is EGNO's
$\mathbb Z_+$-ring.

For a split semisimple rigid category, duality induces a basis permutation
$b_i^*=[S_i^*]$ and an anti-involution $(xy)^*=y^*x^*$. Let $I_0$ be the
indices occurring in the unit and set
$\tau(\sum_i a_ib_i)=\sum_{i\in I_0}a_i$. Then
```math
\label{eq:based-ring-trace}
\tau(b_i b_j)=\delta_{i,j^*}.
```
These are the **based-ring** conditions. A based ring of finite rank is a
**multifusion ring**; it is a **fusion ring** when the unit is a basis element
[EGNO; Definitions 3.1.3 and 3.1.7, pp. 49--50](@cite). In particular, the
Grothendieck ring of a fusion category is a fusion ring
[EGNO; Proposition 4.9.1, pp. 76--77](@cite).

With simple unit $b_u=1$, the last equation reads
$N_{ij}^{\,u}=\delta_{i,j^*}$.
A braiding makes the ring commutative, but commutativity alone does not specify
a braiding.

## Fusion matrices and Frobenius–Perron dimensions

For fixed $i$, the slice `N[i,:,:]` has entries
```math
\label{eq:fusion-matrix}
(M_i)_{j,l}=N_{ij}^{\,l}.
```
It is the matrix of left multiplication by $b_i$ in row coordinates:
if $x$ has coefficient row $v$, then $b_ix$ has coefficient row $vM_i$.

For a fusion ring, the **Frobenius–Perron dimension** of $b_i$ is the
spectral radius of $M_i$. Extending additively gives the unique ring
homomorphism to $\mathbb R$ taking strictly positive values on the distinguished
basis [EGNO; Definition 3.3.3 and Proposition 3.3.6, pp. 53--54](@cite).
Consequently,
```math
\label{eq:fpdim-multiplicativity}
\operatorname{FPdim}(X\otimes Y)
=\operatorname{FPdim}(X)\operatorname{FPdim}(Y).
```
These dimensions depend only on the fusion rules. They need not lie in the
category's base field and do not require a pivotal structure. In contrast,
`dim(X)` uses the chosen duality and pivotal data.

The package computes `fpdim(r)` from the integer multiplication matrices
and returns an exact algebraic number. For a fusion ring $R$,
```math
\label{eq:based-ring-fpdim}
\operatorname{FPdim}(R)=\sum_i\operatorname{FPdim}(b_i)^2,
```
which is also `fpdim(C)` for a split fusion category with Grothendieck ring $R$.

## Non-split categories

If $D_l=\operatorname{End}(S_l)$ is larger than the base field $k$, the
multiplicity is
```math
\label{eq:nonsplit-fusion-multiplicity}
N_{ij}^{\,l}
=\frac{\dim_k\operatorname{Hom}(S_l,S_i\otimes S_j)}
       {\dim_k D_l}.
```
This is the division used by `multiplication_table(C)` and
`coefficients(X,simples(C))` in the generic semisimple implementation.

Let $u$ denote the simple unit and put
$D_u=\operatorname{End}(\mathbb 1)$. Rigidity gives
```math
\label{eq:nonsplit-unit-coefficient}
N_{ij}^{\,u}=\frac{\dim_kD_i}{\dim_kD_u}\,\delta_{i,j^*}.
```
Thus the coefficient of the unit in $[S_i][S_i^*]$ can exceed one.
The resulting ring is a **weak fusion ring**
[EGNO; §3.8, p. 63](@cite); see also
[maurer2024computing; §2.1](@citet).

For such a category the package uses
```math
\label{eq:nonsplit-category-fpdim}
\operatorname{FPdim}(\mathcal C)
=\sum_i
 \frac{\dim_kD_u}{\dim_kD_i}\operatorname{FPdim}(S_i)^2.
```
Under the separability convention for categories over arbitrary fields, this
is the categorical Frobenius--Perron dimension of
[sanford2025fusion; Theorem 3.19 and Definition 3.21](@citet). When the unit is
scalar, it is equation (5.3) of
[maurer2024computing; §5.3](@citet). The category method `fpdim(C)` includes
these endomorphism-algebra factors for every supported semisimple input. The
ring method `fpdim(R)` always sums the squares without these factors, so it does
not give `fpdim(C)` in the non-split case.

## What decategorification forgets

Decategorification retains tensor-product multiplicities but forgets the actual
objects, morphisms, associator, and any additional monoidal structure. A
**categorification** of a based ring is a category whose Grothendieck ring is
the given based ring. A based ring can have inequivalent categorifications or
none at all [EGNO; §§4.9--4.10](@cite).

For example, $\operatorname{Gr}(\operatorname{Vec}_k(G))=\mathbb Z[G]$, with
one basis element for each group element.
A 3-cocycle twist changes the associator but leaves this ring unchanged.
For $\operatorname{Rep}_k(G)$, decategorification gives the representation
ring, computed from tensor products of representations. These computations use
the concrete categorical models introduced above.

An exact tensor functor induces a ring homomorphism
$[X]\mapsto[F(X)]$. The converse is false: a homomorphism between
Grothendieck rings does not determine a functor between the categories.

Continue with [Computing with fusion rings](@ref computing-fusion-rings) for
the package interface and examples, including a non-split representation ring.
