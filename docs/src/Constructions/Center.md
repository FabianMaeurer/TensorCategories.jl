# [Drinfeld centers and half-braidings](@id center)

The center is implemented in terms of objects and morphisms of the input
category. It does not require an $F$-symbol presentation. For the construction
see [EGNO; §7.13](@citet), with the convention translation stated below. The
structural theorem, algorithm, and splitting step over arbitrary fields are
developed in
[maurer2024computing; §§2.5 and 4--5](@citet) and
[maeurer2026thesis; Chapter 2](@citet).

## The direction of a half-braiding

A central object in the package is a pair $(Z,\gamma)$, with natural
isomorphisms
```math
\label{eq:center-half-braiding-component}
\gamma_X:Z\otimes X\longrightarrow X\otimes Z.
```
With $a_{X,Y,Z}:(X\otimes Y)\otimes Z\to X\otimes(Y\otimes Z)$, the
compatibility is
```math
\label{eq:center-half-braiding-hexagon}
\gamma_{X\otimes Y}
=a^{-1}_{X,Y,Z}\circ(\mathrm{id}_X\otimes\gamma_Y)
 \circ a_{X,Z,Y}\circ(\gamma_X\otimes\mathrm{id}_Y)
 \circ a^{-1}_{Z,X,Y}.
```
The unit condition is $\gamma_{\mathbb 1}=\operatorname{id}_Z$ under the
package's normalized unit identifications. This is the convention of
[muger2003subfactorsII; Definition 3.1](@citet),
[maurer2024computing; Equation (1.2)](@citet) and
[maeurer2026thesis; Definition 2.1.1](@citet): the central object $Z$ occurs
first in the source. By contrast, [EGNO; Definition 7.13.1](@citet) put the
ambient object first and write
```math
\label{eq:egno-half-braiding-component}
\widetilde\gamma_X:X\otimes Z\longrightarrow Z\otimes X.
```
The conventions are related by
$\widetilde\gamma_X=\gamma_X^{-1}$. This componentwise inversion also reverses
the displayed braiding. The package uses
$c_{(Z,\gamma),(W,\delta)}=\gamma_W$, whereas
[EGNO; Proposition 8.5.1, Equation (8.15)](@citet) use
$c_{(Z,\widetilde\gamma),(W,\widetilde\delta)}=\widetilde\delta_Z=(\delta_Z)^{-1}$.
Thus, under the direct identification by
inverse half-braidings, the two conventions give reverse braided structures.

A morphism $f:(Z,\gamma)\to(Z',\gamma')$ is a morphism $f:Z\to Z'$ satisfying
```math
\label{eq:center-morphism-condition}
(\mathrm{id}_X\otimes f)\circ\gamma_X
 =\gamma'_X\circ(f\otimes\mathrm{id}_X).
```

`center(C)` constructs the center parent without enumerating its simple
objects. Explicit `CenterObject`s can therefore be formed whenever the ambient
category supplies the operations needed for their half-braidings.
When the input is already modular, `simples(center(C))` uses its braiding to
construct the center simples directly. Otherwise it invokes the substantially
stronger induction algorithm and decomposes endomorphism algebras. For input
$X$, central induction has underlying object

```math
\label{eq:center-induction-object}
\bigoplus_{S\in\operatorname{Irr}(\mathcal C)}
(S\otimes X)\otimes S^*.
```

Here $\operatorname{Irr}(\mathcal C)$ denotes the chosen representatives of
the simple objects.

The enumeration algorithm applies this construction to induction generators
and obtains simple central summands from their endomorphism algebras. The
mathematical algorithm is proved for a pivotal fusion category in the split
sense used in this manual, with $\dim(\mathcal C)\ne0$
[maurer2024computing; Assumption 4.1](@cite). The present implementation's
general induction branch uses chosen spherical morphisms in its projection
formulas. Its supported input is therefore the narrower class of
[split spherical fusion categories](@ref tensor-conventions), again with
$\dim(\mathcal C)\ne0$ in the coefficient field. It is not an enumeration
algorithm for arbitrary monoidal categories.

In positive characteristic these hypotheses require separate verification.
For a pivotal fusion category over any field, the center is weak fusion exactly
when $\dim(\mathcal C)\ne0$; if the center is split, it is modular
[maurer2024computing; Theorem 2.1](@cite). Thus semisimplicity of the input alone
does not imply semisimplicity of its center, and the center need not be split
even when the input is split.

Once the simple objects of $\mathcal Z(\mathcal C)$ over $k$ have been
computed, one chooses a common splitting field for their endomorphism algebras.
After scalar extension, primitive idempotents in these algebras give the split
simple central summands. This is an application of the general
[splitting procedure](@ref splitting-and-scalars), and is
Algorithm 4 of [maurer2024computing; §5.2](@citet). It is a separate stage from
computing the center over $k$: a center can be fully computed while some of
its simple objects remain non-split over that field.

For supported `CenterCategory` models, `split(Z)` searches for a common
splitting field and returns the extended center together with the field
embedding; the field found need not be minimal. If a target field $K$ and an
embedding have already been chosen, `extension_of_scalars(Z,K;
embedding=...)` extends the center and decomposes the known non-split simple
central objects.

[Numerical center computations](@ref numerical-centers) use the same
public interface, with structural equations interpreted in the working ball
field.

## Constructing half-braidings by hand

For $\operatorname{Vec}_{k}(C_2)$ in characteristic zero, the four simple
central objects are indexed by a degree $g\in C_2$ and a character
$\chi:C_2\to k^\times$. On an object homogeneous of degree $h$, the
half-braiding is $\chi(h)$ times the usual interchange map.

Here is this construction over $\mathbb Q$. We do not enumerate the center to
construct these objects.

The call to `braiding` below supplies the ordinary interchange map of this
particular symmetric ambient category. The center construction itself does not
require the ambient category to be braided.

```@example halfbraiding
using TensorCategories, Oscar
G = cyclic_group(2)
C = graded_vector_spaces(QQ,G)
S = simples(C)
Z = center(C)
central = [
    CenterObject(Z, X,
        [braiding(X,Y) * (Y == one(C) ? QQ(1) : QQ(epsilon)) for Y in S])
    for X in S for epsilon in (1,-1)
]
@assert length(central) == 4
for X in central
    @assert is_central(X)
    @assert all(is_invertible, half_braiding(X))
    @assert half_braiding(X, one(C)) == id(object(X))
end
[int_dim(Hom(X,Y)) for X in central, Y in central]
```

The displayed matrix is the identity: these objects are simple and pairwise
nonisomorphic. Their completeness follows independently from the classification
by the two degrees and two characters. We can also test the equation on a
non-simple object:

```@example halfbraiding
X = central[end]
Y = S[1] ⊕ S[2]
W = S[2]
rhs = inv_associator(Y,W,object(X)) ∘
      (id(Y) ⊗ half_braiding(X,W)) ∘ associator(Y,object(X),W) ∘
      (half_braiding(X,Y) ⊗ id(W)) ∘ inv_associator(object(X),Y,W)
@assert half_braiding(X,Y⊗W) == rhs
matrix(half_braiding(X,Y))
show(stdout, MIME"text/plain"(), matrix(half_braiding(X,Y))); println() # hide
```

`CenterObject(Z,X,gamma)` stores components in the order `simples(C)`;
`object` forgets the half-braiding, `half_braiding(X)` returns the stored list,
and `half_braiding(X,Y)` extends it to a general object.
The constructor stores the components without validation. `is_central` checks
the coherence equations; invertibility and the normalized unit condition
$\gamma_{\mathbb 1}=\operatorname{id}_X$ are additional requirements on the supplied components.
For split semisimple input, naturality on simple objects is automatic from
scalar simple endomorphisms. With non-split simples, the current check does not
separately test naturality against their non-scalar endomorphisms, so manually
entered data require that verification as well.

The braiding in the center is
$c_{(X,\gamma),(Y,\delta)}=\gamma_Y$. Thus two central objects with the same
underlying object can be different. The forgetful functor is faithful and
tensor, but its target is $\mathcal C$, not necessarily vector spaces.

## Extracting structural data

Once a center is split, use `six_j_category(Z)` to construct a skeleton,
or `six_j_symbols(Z)` to compute its associators. Use the same binary Hom
bases when computing associators and braiding; independently chosen bases need
not give a compatible pair of $F$- and $R$-symbols.

## [Numerical centers and skeletonization](@id numerical-centers)

A numerical center computation follows the same categorical construction as an
exact one. The scalar comparisons and linear algebra are instead performed in
the category's arbitrary-precision ball field, with the interpretation given
under [Numerical fusion categories](@ref numerical-fusion-categories). The
orthonormal-basis construction described below and its numerical application
to the Haagerup center are developed by
[maeurer2026thesis; §§4.2.4--4.3 and 5.2.3](@citet).

Under the pivotal-fusion hypotheses above, semisimplicity of the center is
equivalent to $\dim(\mathcal C)\ne0$. The center-specific predicate
`is_semisimple(Z)` currently tests this condition by structural comparison when
$Z=\mathcal Z(\mathcal C)$.
If this dimension is a ball containing zero, test
`!contains_zero(dim(category(Z)))` before relying on that predicate. The
half-braiding and skeletal coherence checks use ball-aware comparisons; the
dimension comparison in `is_fusion(Z)` does so for `ArbField` and `AcbField`.

```julia
C = numeric(anyonwiki(3,1,0,1,1,1,1),256)
Z = center(C)
simples(Z)
S = skeletonize(Z)
```

For unitary input, numerical skeletonization chooses an orthonormal basis in
each binary fusion space with respect to the dagger inner product. It uses the
same chosen bases to compute both the associator and the braiding, so the
resulting $F$- and $R$-symbols belong to one common gauge. The skeletal category
can then be checked with the ordinary structural predicates:

```julia
@assert is_unitary(S)
@assert is_modular(S)
@assert pentagon_axiom(S)
@assert hexagon_axiom(S)
```

Its numerical symbol dictionaries can be written with the
[numerical data-exchange format](@ref symbol-data).

The complete pentagon test examines every quadruple of simple objects and can
be expensive. `randomized_pentagon_axiom(S,n)` checks $n$ randomly selected
quadruples and is useful during a long computation, but it is not a replacement
for the complete test when final validation is feasible.

For a pivotal fusion category with $\dim(\mathcal C)\ne0$, the center satisfies

```math
\label{eq:center-global-dimension}
\dim\mathcal Z(\mathcal C)=(\dim\mathcal C)^2.
```

The method `dim(Z)` for a `CenterCategory` returns the right-hand side by
construction. Once the center is split, check an enumeration by comparing this
value with
$\sum_{T\in\operatorname{Irr}(Z)}\dim(T)\dim(T^*)$; the center-specific
`is_fusion(Z)` method makes this comparison after checking that the simples are
split. Further checks are the fusion rules, unitarity when applicable,
nondegeneracy of the $S$-matrix, and the pentagon and hexagon equations. These
checks probe different parts of the construction. The displayed identity uses
the pivotal categorical dimension; sphericality is not required
[maurer2024computing; Theorem 2.1](@cite). Independently, in the usual
characteristic-zero fusion setting,
$\operatorname{FPdim}(\mathcal Z(\mathcal C))=\operatorname{FPdim}(\mathcal C)^2$
[EGNO; Theorem 7.16.6](@cite).

Continue with [The Ising center over two fields](@ref ising-center), which
follows a non-split center through scalar extension and extraction of a single
compatible skeleton. Precomputed database centers are described under
[AnyonWiki](../F-symbols/AnyonWiki.md); a saved skeletal center does not retain
the original explicit half-braidings.
