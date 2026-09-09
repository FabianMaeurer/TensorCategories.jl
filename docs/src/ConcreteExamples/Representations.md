# [Representations of finite groups](@id representations)

Finite-group representations are a fundamental source of symmetric tensor
categories and the model behind ordinary finite-group symmetry. The
constructor realizes $\operatorname{Rep}_K(G)$ as in
[EGNO; Examples 2.3.4 and 2.10.13, pp. 26 and 43](@citet): the tensor product
uses the diagonal $G$-action, the unit is the trivial representation, and the
symmetry is the ordinary flip. The package uses row coordinates, so action
matrices and intertwiners multiply on the right. Explicitly, a stored row
vector transforms by $v\mathbin{\cdot}g=v\rho(g)$, and an intertwiner
$M:X\to Y$ satisfies
$\rho_X(g)M=M\rho_Y(g)$. This right-action convention is equivalent to the
usual left-action convention after replacing $g$ by $g^{-1}$.

`representation_category(K,G)` models finite-dimensional representations of a
finite group over the specified field. Objects store a group homomorphism into
a matrix group. Morphisms are intertwiners in the row-vector convention.
Use either `Representation(C,generators,matrices; check=false)` for an existing
category $C$, or `Representation(G,generators,matrices; check=false)` to infer
the field and parent category from the matrices. The corresponding overloads
with a Julia function evaluate that function on the package's generators.

## Supplying action matrices

```@example representations
using TensorCategories, Oscar
G = cyclic_group(3)
C = representation_category(QQ,G)
A = matrix(QQ,[0 1; -1 -1])
X = Representation(C,gens(G),[A]; check=true)
@assert A^3 == identity_matrix(QQ,2)
@assert int_dim(End(X)) == 2
@assert is_simple(X)
f = morphism(X,X,A; check=true)
@assert matrix(f ∘ f) == A^2
int_dim(X ⊗ X)
```

The keyword defaults to `check=false` for both constructors.
`check=true` on `Representation` checks the group relations, while on
`morphism` it checks equivariance. Parent categories, endpoint dimensions, and
coefficient fields are checked regardless. Supply a generating set of $G$ and
one square action matrix for each generator, with a common size and coefficient
field. A generator image may have smaller order than the generator.

Tensor products use diagonal group actions; duals use contragredient actions;
the usual flip is symmetric in every characteristic. The canonical spherical
structure is implemented. Tensor-product matrices use Kronecker products, with
the coordinate from the right tensor factor varying fastest. No $F$-symbols are
required.

Equality of representations compares their parent categories and their action
matrices on the package's generators; representations related by a nontrivial
change of basis are generally only isomorphic. Use `is_isomorphic(X,Y)` for
that comparison. Two `GroupRepresentationCategory` values compare equal when
their stored groups and coefficient fields compare equal.

## Fields and enumeration

Maschke's theorem gives semisimplicity when the characteristic does not divide
$|G|$. Splitting is a further condition. The rational example above is simple
with a two-dimensional endomorphism field, so it is not absolutely simple.
Finiteness of $G$ is a hypothesis of this model and is not checked by the
constructor.

The no-field constructor `representation_category(G)` uses OSCAR's abelian
closure of `QQ`, which is a splitting field for every finite group. Ordinary
irreducible representations can therefore be enumerated:

```@example representations
D = representation_category(symmetric_group(3))
int_dim.(simples(D))
```

Over a nonsplitting characteristic-zero field, the absolutely irreducible
matrices returned by GAP need not be defined over the requested field. The
implementation checks their entries before constructing any objects and fails
rather than silently changing the coefficient field. Explicit representations,
Hom spaces, and exact tests such as `is_simple(X)` remain available in the
cases described below.

## Representation algorithms and backends

The functions `simples`, `is_simple`, `composition_factors`,
`simple_subobjects`, and `decompose` accept a `backend` keyword. The available
values are:

| Value | Behaviour |
|:---|:---|
| `:auto` | use the established default for the coefficient field and operation |
| `:gap` | use GAP's ordinary representation routines or its finite-field MeatAxe |
| `:hecke` | use Hecke's matrix-module MeatAxe routines |

For finite fields, `:auto` uses GAP. It supports irreducibility tests, simple
enumeration, composition factors, and indecomposable decomposition. GAP's
ordinary routine supplies absolutely irreducible representations in
characteristic zero. This works directly over the default abelian closure and
over another field when GAP's chosen matrices can be converted to that field.
It does not by itself construct the simple objects over an arbitrary
nonsplitting field; Galois orbits and Schur indices enter that problem.

The Hecke backend works directly with the stored generator matrices. It is
available over finite fields and for one-generator modules over supported
infinite exact fields. For several generators over an infinite field, the
current Hecke algorithm does not provide a dependable backend and the package
reports this limitation. In a semisimple category, the default implementation
can compute composition multiplicities from Hom spaces once `simples(C)` is
available. In modular characteristic, composition factors and indecomposable
direct summands are distinct; `decompose(X; backend=:hecke)` is therefore not
used there. GAP's relevant matrix-module operations are documented in
[gapmanual2026; Chapters 69.5 and 69.7](@cite), and the Hecke matrix-module
implementation is part of the Nemo/Hecke system described by
[fieker2017nemo](@cite).

The following rational example has one generator, so the Hecke backend can
recover the two rational simple modules of $C_3$ from the regular
representation:

```@example representations
Cq = representation_category(QQ, cyclic_group(3))
int_dim.(simples(Cq; backend=:hecke))
```

## Splitting fields

The predicate `is_split_semisimple(C)` tests both semisimplicity and whether
the endomorphism algebra of every simple is the coefficient field. Over
$\mathbb Q$, the implementation instead uses the equivalent ordinary-character
criterion: all irreducible characters must be rational-valued and have Schur
index one.

If $m$ is the exponent of $G$, Brauer's splitting theorem states that a field
containing a primitive $m$-th root of unity is a splitting field for $G$. In
characteristic $p>0$, only the prime-to-$p$ part of $m$ contributes a
nontrivial root of unity. The function `splitting_field(C)` constructs the
splitting field of the corresponding polynomial $x^{m'}-1$ over the current
exact coefficient field, where $m'=m$ in characteristic zero and $m'$ is the
prime-to-$p$ part of $m$ in characteristic $p$; compare
[webb2016representations; Theorem 9.2.7](@cite).

!!! note
    `splitting_field(C)` returns **a** splitting field for $G$ over the current
    coefficient field. It does not claim to return a minimal splitting field.

```@example representations
C2 = representation_category(GF(2), cyclic_group(3))
L = splitting_field(C2)
@assert degree(L) == 2
@assert !is_split_semisimple(C2)
@assert is_split_semisimple(representation_category(L, cyclic_group(3)))
L
```

## Indecomposable representations

The predicate `is_finite_representation_type(C)` records whether there are
only finitely many indecomposable isomorphism classes. In characteristic zero
and in nonmodular characteristic this follows from Maschke's theorem. If the
coefficient field has characteristic $p$ dividing $|G|$, Higman's theorem says
that the representation type is finite exactly when a Sylow $p$-subgroup of
$G$ is cyclic [higman1954indecomposable](@cite).

When the category is semisimple, `indecomposables(C)` is the same list as
`simples(C)`. The installed GAP and Hecke backends do not provide a general
enumeration of every indecomposable module in the modular finite-type case, so
the predicate is available there but enumeration is not.

```@example representations
@assert is_finite_representation_type(
    representation_category(GF(5), cyclic_group(5)))
@assert !is_finite_representation_type(
    representation_category(GF(2), symmetric_group(4)))
nothing # hide
```

Semisimplicity, splitting, and the distinction between simple and absolutely
simple objects follow the conventions of [EGNO; §§4.2 and 4.16](@citet).
