# [Anyons, CFT, and tensor-category language](@id physics-bridge)

!!! note "Optional terminology bridge"
    This page is for readers coming from anyons or rational conformal field
    theory (RCFT). It translates the categorical structures and skeletal
    coordinates introduced above into common physics terminology. Other readers
    may proceed directly to numerical fusion categories.

The categorical language follows [EGNO](@citet). The anyon conventions and
fusion-tree language follow [bonderson2007thesis; Chapter 2](@citet).

## Translation of terminology

In a unitary anyon model one usually works with a unitary braided fusion
category equipped with compatible duality and spherical or ribbon data. The
general package interface does not assume that every category has all these
structures.

| Anyon or CFT terminology | Tensor-category and package terminology |
|:---|:---|
| topological charge, particle type, primary label | simple object $a$ |
| vacuum charge | tensor unit $\mathbb 1$, returned by `one(C)` |
| antiparticle or conjugate charge | dual object $a^*$, returned by `dual(a)` |
| fusion | tensor product $a\otimes b$ |
| fusion multiplicity | $N_{ab}^{c}=\dim_k\operatorname{Hom}(a\otimes b,c)$ in the split case |
| fusion channel $c$ | simple summand $c$ of $a\otimes b$ |
| fusion space | $V_{ab}^{c}=\operatorname{Hom}(a\otimes b,c)$ in the package's projection convention |
| splitting space | $\operatorname{Hom}(c,a\otimes b)$, composition-dual to a chosen projection basis |
| fusion tree | basis obtained by iterating binary fusion spaces |
| $F$-move | associator written in two fusion-tree bases |
| exchange or elementary braid | braiding; its binary fusion-space matrix gives the $R$-symbols |
| topological spin | twist eigenvalue $\theta_a$ |
| quantum dimension | categorical dimension, returned by `dim(a)` for the chosen pivotal structure |
| fusion-rule dimension | Frobenius--Perron dimension, returned by `fpdim(a)` |
| total quantum dimension $\mathcal D$ | $\sqrt{\dim(\mathcal C)}$ for a unitary spherical fusion category |
| modular $S$ and $T$ data | `normalized_smatrix(C)` and `tmatrix(C)`, with the normalization below |

For a split fusion category, the integers $N_{ab}^{c}$ are the structure
constants of its [Grothendieck ring](@ref grothendieck-rings). Over a
non-splitting field, simple endomorphism division algebras enter the
multiplicity formula. The distinction between categorical and
Frobenius--Perron dimensions also matters when the chosen spherical dimensions
are not the positive Frobenius--Perron dimensions.

## First anyon computation: Ising

The default constructor uses $\mathbb Q(\sqrt2)$, which contains the
associator and pivotal coefficients needed below but not the phases of a
braided realization. The [Ising catalogue entry](../F-symbols/TambaraYamagami.md)
describes the coefficient-field choices. The package labels its three simple
objects $(\mathbb 1,\chi,X)$; in the example we rename these as the conventional
charges $(\mathbb 1,\psi,\sigma)$:

```@example physics
using TensorCategories, Oscar
C = ising_category()
vacuum, psi, sigma = simples(C)
@assert vacuum == one(C)
@assert psi ⊗ psi == vacuum
@assert psi ⊗ sigma == sigma ⊗ psi == sigma
@assert sigma ⊗ sigma == vacuum ⊕ psi
@assert dim(sigma)^2 == 2
@assert dim(C) == 4
dim.(simples(C))
```

Thus the package's global dimension is

```math
\label{eq:physics-ising-global-dimension}
\dim(\mathcal C)=\sum_a d_a^2=4.
```

Physicists usually call $\mathcal D=2$ the total quantum dimension. The method
`dim(C)` returns $\mathcal D^2$, whereas `dim(sigma)` returns the quantum
dimension $d_\sigma$.

## Fusion spaces are not object dimensions

For three Ising anyons with total charge $\sigma$,

```math
\label{eq:physics-ising-four-anyon-space}
\mathcal H^{\sigma}_{\sigma\sigma\sigma}
=\operatorname{Hom}((\sigma\otimes\sigma)\otimes\sigma,\sigma)
```

has dimension two. Its usual fusion-tree basis is indexed by the intermediate
charges $\mathbb 1$ and $\psi$:

```@example physics
H = Hom((sigma ⊗ sigma) ⊗ sigma, sigma)
@assert int_dim(H) == 2
@assert size(matrix(id(sigma))) == (1,1)
int_dim(H)
```

The $1\times1$ matrix representing `id(sigma)` does not say that
$d_\sigma=1$. This matrix records the multiplicity coordinates of the identity
morphism in the package's split semisimple realization; quantum dimension is
categorical trace data. The distinction is developed under
[Fiber functors and semisimple coordinates](@ref fiber-functors).

## $F$- and $R$-symbols

An $F$-matrix represents the associator between two fusion-tree bases. In the
package's projection convention, an $R$-matrix represents the pullback by the
braiding

```math
\label{eq:physics-braiding-pullback}
c_{a,b}^{*}:V_{ba}^{c}\longrightarrow V_{ab}^{c},
\qquad p\longmapsto p\circ c_{a,b}.
```

Equivalently, the braiding itself maps the composition-dual splitting space
$\operatorname{Hom}(c,a\otimes b)$ to
$\operatorname{Hom}(c,b\otimes a)$. Consequently, the entries of $F$- and
$R$-matrices depend on the chosen projection and splitting bases, on the
direction of the structural maps, and on the ordering of multiplicity indices.
They are not determined by the fusion rules alone.

The [skeletal fusion model](@ref skeletal-fusion) defines these symbols from
the structural morphisms and chosen fusion bases. Those definitions, rather
than the names of array entries, must be used when translating formulas or data
from the physics literature.

## $S$, $T$, topological spins, and CFT conventions

The formal hypotheses and definitions appear under
[Premodular and modular categories](@ref premodular-categories). In those
hypotheses, `smatrix(C)` returns the unnormalized pivotal trace of double
braiding, while the package uses

```math
\label{eq:physics-modular-data-normalization}
S^{\mathrm{norm}}=\frac{1}{\sqrt{\dim(\mathcal C)}}S,
\qquad
T^{\mathrm{cat}}=\operatorname{diag}(\theta_a).
```

These are the categorical twist and $S$-matrix conventions of
[EGNO; §§8.10 and 8.13, especially Eq. (8.46)](@citet).

The square root in `normalized_smatrix(C)` is a choice over a general
coefficient field. In a unitary realization one chooses the positive total
quantum dimension.

For RCFT characters, a common convention is

```math
\label{eq:physics-rcft-t-matrix}
T^{\mathrm{RCFT}}_{aa}
=\exp\!\left(2\pi i\left(h_a-\frac{c}{24}\right)\right)
=e^{-2\pi i c/24}\theta_a,
\qquad \theta_a=e^{2\pi i h_a}.
```

Thus `tmatrix(C)` is the categorical twist matrix, without the overall
central-charge phase. A unitary modular category constrains the topological
central charge through its modular data, but it does not determine a full RCFT
central charge or conformal weights beyond their categorical phase data. The
relation between RCFT fusing and braiding data and categorical modular data is
described by [moore1989classical](@citet).

## Exact and numerical use

Fusion data may be stored exactly over a number field or evaluated over complex
balls at a chosen working precision. Exact coefficients can have several
complex embeddings, and different embeddings can change unitarity or phases.
The [splitting and scalar-extension chapter](@ref splitting-and-scalars)
explains these choices. The
[numerical-computation chapter](@ref numerical-computations) explains ball
arithmetic, and the later
[numerical fusion-category chapter](@ref numerical-fusion-categories) applies
it to coherence, unitarity, and modularity. Numerical center computations are
introduced in the [Drinfeld-center chapter](@ref numerical-centers), after
half-braidings have been defined.

## Scope and current limitations

A modular tensor category records topological and chiral categorical data; it
is not by itself a complete two-dimensional CFT. TensorCategories.jl does not
construct characters, position-dependent conformal blocks,
operator-product-expansion (OPE) coefficients,
or modular-invariant partition functions. Full RCFT constructions require
additional data; for example, [fuchs2002tft](@citet) use a symmetric special
Frobenius algebra in a modular tensor category.

The package supplies associators, braidings, and local $F$- and $R$-symbols, from
which braid actions can be assembled. It currently has no public high-level API
that accepts a braid word and returns its matrix in a chosen multi-anyon
fusion-tree basis. The ordinary `braiding(X,Y)` method is a categorical
structural map, not such a braid-word interface.

Continue with [numerical fusion categories](@ref numerical-fusion-categories).
