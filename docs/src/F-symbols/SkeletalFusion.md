# [Skeletal fusion categories and symbol conventions](@id skeletal-fusion)

Let $\mathcal C$ be a split fusion category over a field $k$, and choose
representatives $S_1,\ldots,S_r$ of its simple objects. Decompositions into
these simples give a skeletal linear model: objects become multiplicity
vectors, morphisms become blocks of matrices over $k$, and the tensor structure
is described by fusion rules and associator matrices in chosen bases.
TensorCategories.jl implements this model as `SixJCategory`. We use the
definitions of semisimple and fusion
categories in [EGNO; Chapters 1, 2, and 4](@citet). The implementation of this
model and its use in the package are described by
[maurer2024computing; §2](@citet).

The name `SixJCategory` refers to this general skeletal data model. It allows
arbitrary fusion multiplicities, and its associator entries need not be literal
Wigner $6j$-symbols.

## Objects and morphisms

An object

```math
\label{eq:skeletal-object-decomposition}
X=\bigoplus_{i=1}^r S_i^{\oplus m_i}
```

is represented by the vector $(m_1,\ldots,m_r)$ of nonnegative integers. In
the implementation this vector is `X.components`, and `C[i]` denotes the
chosen representative $S_i$. The zero object has all multiplicities zero, and
direct sums add multiplicity vectors.

To keep formulas and code legible, the manual also uses a simple label such as
$a$ for its position in `simples(C)` when it occurs inside an array access.
Thus `C.ass[a,b,c,d]` means the block indexed by the positions of the four
simples $a,b,c,d$; Julia code must of course supply the corresponding integers.

A `SixJCategory` is mutable, and its objects retain that particular category as
their parent. Two independently constructed categories can carry identical
arrays without being the same parent; transport objects explicitly between
them.

If

```math
\label{eq:skeletal-target-decomposition}
Y=\bigoplus_{i=1}^r S_i^{\oplus n_i},
```

then

```math
\label{eq:skeletal-homspace}
\operatorname{Hom}_{\mathcal C}(X,Y)
\cong\bigoplus_{i=1}^r\operatorname{Mat}_{m_i\times n_i}(k).
```

A morphism $f:X\to Y$ is therefore stored as one
$m_i\times n_i$ matrix for each simple $S_i$. These matrices act on row
coordinates, so the block representing $g\circ f$ is the block for $f$
multiplied by the block for $g$.

This description uses
$\operatorname{End}_{\mathcal C}(S_i)=k$. In a non-split semisimple category,
the blocks have coefficients in the division algebras
$\operatorname{End}_{\mathcal C}(S_i)$ instead. `SixJCategory` implements the
split case; see [scalar extension and splitting](@ref splitting-and-scalars)
for the non-split setting.

The matrix blocks describe maps between multiplicity spaces. Their existence
does not require a monoidal fiber functor
$\mathcal C\to\operatorname{Vec}_k$; see
[Fiber functors and semisimple coordinates](@ref fiber-functors).

## Fusion rules

The fusion multiplicities are the structure constants of the
[Grothendieck ring](@ref grothendieck-rings) in its simple-object basis:

```math
\label{eq:skeletal-fusion-rule}
S_i\otimes S_j\cong
\bigoplus_l S_l^{\oplus N_{ij}^{\,l}},
\qquad
N_{ij}^{\,l}
=\dim_k\operatorname{Hom}_{\mathcal C}(S_i\otimes S_j,S_l).
```

The integer array `C.tensor_product[i,j,l]` stores $N_{ij}^{\,l}$. By
bilinearity,

```math
\label{eq:skeletal-tensor-multiplicity}
[X\otimes Y:S_l]=\sum_{i,j}m_i n_jN_{ij}^{\,l}.
```

The simple unit has a distinguished index $u$ and satisfies

```math
\label{eq:skeletal-unit-fusion}
N_{ui}^{\,j}=N_{iu}^{\,j}=\delta_{ij}.
```

The program stores this index separately through `set_one!`.

For morphisms, tensor products use Kronecker products of matrix blocks, with
one copy for each binary fusion channel. Coordinates in these copies require a
choice of basis in every space
$\operatorname{Hom}_{\mathcal C}(S_i\otimes S_j,S_l)$.

Associativity of the fusion rules requires

```math
\label{eq:skeletal-associativity-dimensions}
\sum_eN_{ab}^{\,e}N_{ec}^{\,d}
=
\sum_fN_{bc}^{\,f}N_{af}^{\,d}.
```

The equality of these dimensions is not yet an associator. To write its
matrix, one must first choose and order bases in the binary fusion spaces.

## [Fusion bases and structural matrices](@id f-conventions)

Scalar structural data arise only after bases have been chosen in the binary
fusion spaces. We now fix those bases, use them to define the matrices of the
associator, braiding, and pivotal structure, and only then name their scalar
entries.

### Binary fusion bases

The implementation uses the fixed order $(S_1,\ldots,S_r)$ returned by
`simples(C)`. For simple objects $a,b,e$, put

```math
\label{eq:binary-fusion-space}
V_{ab}^{e}=\operatorname{Hom}_{\mathcal C}(a\otimes b,e),
\qquad
N_{ab}^{e}=\dim_k V_{ab}^{e},
```

and choose an ordered projection basis

```math
\label{eq:binary-projection-basis}
p^{ab}_{e,\mu}:a\otimes b\longrightarrow e,
\qquad 1\leq \mu\leq N_{ab}^{e}.
```

Let

```math
\label{eq:binary-splitting-basis}
s_{ab}^{e,\mu}:e\longrightarrow a\otimes b
```

be the composition-dual splitting basis. Thus

```math
\label{eq:composition-dual-bases}
p^{ab}_{e,\mu}\circ s_{ab}^{e,\nu}
=\delta_{\mu,\nu}\operatorname{id}_e,
\qquad
\sum_{e,\mu}s_{ab}^{e,\mu}\circ p^{ab}_{e,\mu}
=\operatorname{id}_{a\otimes b}.
```

These binary bases determine bases for the two projection spaces of a triple
tensor product. For a fixed output simple $d$, define

```math
\label{eq:left-triple-projection}
L_{e,\mu,\nu}
=p^{ec}_{d,\nu}\circ
  \bigl(p^{ab}_{e,\mu}\otimes\operatorname{id}_c\bigr)
:(a\otimes b)\otimes c\longrightarrow d
```

and

```math
\label{eq:right-triple-projection}
R_{f,\rho,\sigma}
=p^{af}_{d,\sigma}\circ
  \bigl(\operatorname{id}_a\otimes p^{bc}_{f,\rho}\bigr)
:a\otimes(b\otimes c)\longrightarrow d.
```

The left paths are ordered lexicographically as $(e,\mu,\nu)$, with $\nu$
varying fastest. The right paths are ordered as $(f,\rho,\sigma)$, with
$\sigma$ varying fastest. Intermediate simples follow the order of
`simples(C)`.

### The associator matrix

Let

```math
\label{eq:associator-map}
\alpha_{a,b,c}:(a\otimes b)\otimes c
   \longrightarrow a\otimes(b\otimes c)
```

be the associator. Its block with output $d$ is the matrix $A=A^{abc}_d$
defined by

```math
\label{eq:associator-projection-bases}
R_{f,\rho,\sigma}\circ\alpha_{a,b,c}
=\sum_{e,\mu,\nu}
 A^{abc}_d[(e,\mu,\nu),(f,\rho,\sigma)]
 L_{e,\mu,\nu}.
```

Thus rows are left paths and columns are right paths. In the documentation's
index notation, the stored block is `C.ass[a,b,c,d]`. Its size is

```math
\label{eq:associator-block-size}
\left(\sum_e N_{ab}^{e}N_{ec}^{d}\right)
\times
\left(\sum_f N_{bc}^{f}N_{af}^{d}\right).
```

Associativity of the fusion rules makes the two numbers equal.

Let $L^{e,\mu,\nu}$ and $R^{f,\rho,\sigma}$ be the triple-product splittings
induced by the composition-dual binary bases. Equation
$\eqref{eq:associator-projection-bases}$ is equivalent to

```math
\label{eq:associator-splitting-bases}
\alpha_{a,b,c}\circ L^{e,\mu,\nu}
=\sum_{f,\rho,\sigma}
 A^{abc}_d[(e,\mu,\nu),(f,\rho,\sigma)]
 R^{f,\rho,\sigma}.
```

The entries of $A^{abc}_d$ are the **$F$-symbols** in this manual. Explicitly,

```math
\label{eq:F-symbol-definition}
\left[F^{abc}_d\right]_{(e,\mu,\nu),(f,\rho,\sigma)}
:=A^{abc}_d[(e,\mu,\nu),(f,\rho,\sigma)].
```

Equation $\eqref{eq:associator-splitting-bases}$ replaces a left-associated
splitting tree by a linear
combination of right-associated splitting trees; this change of basis is the
**$F$-move**. The first multi-index in $\eqref{eq:F-symbol-definition}$ labels
the input tree and the second labels the output tree. At the coordinate-free
level, this is the associativity isomorphism on the splitting spaces
$H^e_{ab}=\operatorname{Hom}(e,a\otimes b)$ in
[EGNO; §4.9, Eqs. (4.12)--(4.13)](@citet). The same direction and coefficient
convention are used in
[bonderson2008interferometry; Eq. (2.14)](@citet) and
[barkeshli2019symmetry; Eq. (10)](@citet). Both references work with unitary
anyon models and orthonormal splitting bases. The same coefficient convention
makes sense for the composition-dual bases used here over an arbitrary
splitting field.

The function `F_symbols(C; convention=...)` returns all admissible
coefficients, including zeros, as a Julia dictionary. The two accepted
dictionary conventions are:

| `convention` | Keys |
|:---|:---|
| `:column_major_packing` | `[a,b,c,d,f,e]` without multiplicities and `[a,b,c,d,f,sigma,rho,e,mu,nu]` in general |
| `:bonderson` | `[a,b,c,d,e,f]` without multiplicities and `[a,b,c,d,e,mu,nu,f,rho,sigma]` in general |

With `convention=:bonderson`, the value is precisely the entry in
$\eqref{eq:F-symbol-definition}$ named by the two fusion paths in the key. The
default is `convention=:column_major_packing`, the layout used by
TensorCategories data files. In that layout the suffix records the order used
to pack the entries of $A$ in Julia column-major order: for fixed $a,b,c,d$,
the keys are traversed in the nested order
$(e,f,\nu,\mu,\rho,\sigma)$ and matched with successive entries of `vec(A)`.
Programmatically, this pairs the column-major list of matrix entries with a
nested enumeration of the admissible basis morphisms. The keys record that
enumeration rather than a pair of direct fusion-path indices. Since the
admissible path lists depend on $(a,b,c,d)$ and need not coincide, conversion
to direct path indices is block-dependent and cannot be given by a fixed
permutation of the key positions. The [data exchange page](@ref symbol-data)
gives the equivalent entry-by-entry reconstruction rule. Changing
`convention` changes the dictionary keys and their interpretation, not the
matrix $A$, the chosen bases, or the gauge.
`numeric_F_symbols` uses the same keyword.

The package stores matrices for row coordinates: whenever the represented
composition is defined,

```julia
matrix(g ∘ f) == matrix(f) * matrix(g)
```

A reader who instead represents splitting-space vectors by columns therefore
uses $A^{\mathsf T}$ as the matrix of the associator. This transpose is only a
coordinate convention. It is neither an inverse nor a complex conjugate, and
composition-dual bases need not be Hermitian-adjoint bases.

### The braiding matrix

For simple objects $a,b,d$, let $B=B^{ab}_d$ be defined by

```math
\label{eq:braiding-projection-bases}
p^{ba}_{d,\nu}\circ c_{a,b}
=\sum_{\mu}
 B^{ab}_d[\mu,\nu]p^{ab}_{d,\mu}.
```

In composition-dual splitting bases this is

```math
\label{eq:braiding-splitting-bases}
c_{a,b}\circ s_{ab}^{d,\mu}
=\sum_{\nu}B^{ab}_d[\mu,\nu]s_{ba}^{d,\nu}.
```

The entries of $B^{ab}_d$ are the **$R$-symbols**:

```math
\label{eq:R-symbol-definition}
\left[R^{ab}_d\right]_{\mu,\nu}:=B^{ab}_d[\mu,\nu].
```

Equation $\eqref{eq:braiding-splitting-bases}$ is the **$R$-move** induced by
the braiding. The row index $\mu$ labels the input splitting vector
$s_{ab}^{d,\mu}$ and its composition-dual projection $p^{ab}_{d,\mu}$; the
column index $\nu$ similarly labels the output splitting vector and dual
projection for $b\otimes a$. This is precisely the convention of
[bonderson2008interferometry; Eqs. (2.30)--(2.32)](@citet), abbreviated BSS
below: their operator $R_{ab}$ sends the splitting basis for $a\otimes b$ to
that for $b\otimes a$, and

```math
\label{eq:bonderson-R-identification}
B^{ab}_d[\mu,\nu]
=\left[R^{ab}_d\right]^{\mathrm{BSS}}_{\mu,\nu}
```

after identifying the chosen splitting bases. With the same index notation,
the structural matrix is `C.braiding[a,b,d]`.

The order in the superscript is not uniform in the physics literature.
[barkeshli2019symmetry; Eqs. (27)--(28)](@citet), abbreviated BBCW below,
define $R^{ab}_d:V^{ba}_d\to V^{ab}_d$. With matching bases, their notation is
therefore related to the package's by

```math
\label{eq:barkeshli-R-identification}
B^{ab}_d[\mu,\nu]
=\left[R^{ba}_d\right]^{\mathrm{BBCW}}_{\mu,\nu}.
```

The function `R_symbols(C; convention=...)` returns all admissible
coefficients, including zeros. The two accepted dictionary conventions are:

| `convention` | Key with multiplicities | Value stored at that key |
|:---|:---|:---|
| `:column_major_packing` | `[a,b,d,mu,nu]` | $B^{ab}_d[\nu,\mu]$ |
| `:bonderson` | `[a,b,d,mu,nu]` | $B^{ab}_d[\mu,\nu]$ |

Without multiplicities, both conventions use the key `[a,b,d]` and return the
single entry of $B^{ab}_d$. Thus the two multiplicity indices are transposed
in the default dictionary layout. This is a packing rule, not an inverse
braiding or a change of gauge. `numeric_R_symbols` uses the same keyword.

Both functions return an ordinary `Dict`, which does not retain the selected
convention. Code that stores or passes the dictionary separately from the
category must therefore retain the convention as accompanying metadata.

The order of $a$ and $b$ matters. The inverse of $B^{ab}_d$ represents the
inverse map from $b\otimes a$ to $a\otimes b$; it is not generally
$B^{ba}_d$.

### [Pivotal coefficients](@id pivotal-coefficients)

A pivotal structure is a monoidal natural isomorphism

```math
\label{eq:pivotal-structure-map}
j:\operatorname{id}_{\mathcal C}\Longrightarrow(-)^{**}
```

[EGNO; Definition 4.7.7](@cite). In the skeletal model the chosen duality has
$S_i^{**}=S_i$. Since $S_i$ is split simple, the component of $j$ is a scalar:

```math
\label{eq:P-symbol-definition}
j_{S_i}=P_i\operatorname{id}_{S_i},\qquad P_i\in k^\times.
```

TensorCategories.jl calls the scalars $P_i$ the **$P$-symbols**. They are stored as
`C.pivotal[i]`, and `P_symbols(C)` returns the dictionary
`Dict([i] => P_i)`. This function name refers specifically to the components
in $\eqref{eq:P-symbol-definition}$; the term “pivotal symbols” is also used elsewhere for different data
attached to trivalent fusion spaces.

The coefficients $P_i$ are additional structure: they are not determined by
the $F$- or $R$-symbols. They must make $j$ monoidal. If

```math
\label{eq:double-dual-tensorator}
\phi_{X,Y}:(X\otimes Y)^{**}\longrightarrow X^{**}\otimes Y^{**}
```

is the package's monoidal structure of the double-dual functor, the required
identity is

```math
\label{eq:pivotal-monoidality}
j_X\otimes j_Y=\phi_{X,Y}\circ j_{X\otimes Y}.
```

The initializer `six_j_category` sets $P_i=1$ for every simple object, and
`set_pivotal!` replaces these components. Validation is explicit: use
`is_pivotal(C; check=true)` to check
$\eqref{eq:pivotal-monoidality}$. Sphericality is the
additional equality of the left and right pivotal traces and can be checked
with `is_spherical(C; check=true)`. Since a $P$-symbol dictionary has one
scalar per simple object rather than matrix entries indexed by fusion paths,
`P_symbols` has no `convention` keyword.

When `skeletonize(C)` constructs a skeletal model from a concrete category,
the ordered simple representatives and the bases of
$\operatorname{Hom}(S_i\otimes S_j,S_k)$ determine the skeletal associator,
braiding, and duality. The same choices must therefore be used when the
pivotal structure is transported. Let $S_i$ be the source representative,
let $\widehat S_i$ be the corresponding skeletal simple, and let
$j^{(0)}_{\widehat S_i}=\operatorname{id}_{\widehat S_i}$ denote the reference
double-dual identification before its pivotal coefficient is changed. The
transported coefficient is

```math
\label{eq:skeletal-pivotal-transport}
P_i=
\frac{\operatorname{Tr}_{L}(j_{S_i})}
     {\operatorname{Tr}_{L}(j^{(0)}_{\widehat S_i})}.
```

This trace ratio is the coordinate of the transported map in the
one-dimensional split-simple double-dual Hom space. Its denominator is
nonzero: in a semisimple tensor category the trace of an isomorphism from a
simple object to its double dual is nonzero
[EGNO; Proposition 4.8.4](@cite). This remains true in positive
characteristic under the split semisimple hypotheses required by
`skeletonize`. With ball-valued coefficients, an enclosure containing zero
means that the chosen working precision does not resolve this nonzero trace;
skeletonization then asks for higher precision instead of retaining default
coefficients.

The values $P_i$ are coordinates in the resulting skeletal gauge. Changing
the multiplicity-space bases can change them, even though the transported
pivotal category remains equivalent. In particular, an all-one pivotal vector
is a statement about a chosen presentation, not an invariant of the category.

### [Unit normalization](@id unit-normalization)

`SixJCategory` uses strict unit constraints in its skeletal coordinates. The
normalized convention requires every associator block with a unit input to be
an identity matrix, and the public `associator` function treats such inputs as
strict. The low-level setters accept `check=false` for prevalidated input. With
`check=true`, `set_one!` and `set_associator!` check unit normalization. The
full pentagon is checked separately by `pentagon_axiom(C)`.

### Pentagon and hexagon equations

In a multiplicity-free category, write
$\mathcal F^{abc}_d[e,f]$ for the coefficient in
$\eqref{eq:associator-splitting-bases}$, and write $\mathcal R^{ab}_d$ for the
coefficient in $\eqref{eq:braiding-splitting-bases}$. With inadmissible fusion paths
interpreted as zero, the pentagon equation is

```math
\label{eq:pentagon-multiplicity-free}
\mathcal F^{fcd}_e[g,l]\,\mathcal F^{abl}_e[f,k]
=\sum_h
  \mathcal F^{abc}_g[f,h]\,
  \mathcal F^{ahd}_e[g,k]\,
  \mathcal F^{bcd}_k[h,l].
```

This is the multiplicity-free specialization of the standard indexed pentagon
equation in [barkeshli2019symmetry; Eq. (12)](@citet). The two hexagon
equations, in the same convention, are

```math
\label{eq:hexagon-positive-multiplicity-free}
\mathcal R^{ca}_e\,
\mathcal F^{acb}_d[e,g]\,
\mathcal R^{cb}_g
=\sum_f
\mathcal F^{cab}_d[e,f]\,
\mathcal R^{cf}_d\,
\mathcal F^{abc}_d[f,g]
```

and

```math
\label{eq:hexagon-negative-multiplicity-free}
(\mathcal R^{ac}_e)^{-1}\,
\mathcal F^{acb}_d[e,g]\,
(\mathcal R^{bc}_g)^{-1}
=\sum_f
\mathcal F^{cab}_d[e,f]\,
(\mathcal R^{fc}_d)^{-1}\,
\mathcal F^{abc}_d[f,g].
```

These agree with [bonderson2007thesis; Eqs. (2.57)--(2.58)](@citet), with the label
order and $R$-symbol superscripts used here.

Equations $\eqref{eq:hexagon-positive-multiplicity-free}$ and
$\eqref{eq:hexagon-negative-multiplicity-free}$ result by suppressing the
multiplicity indices in those full equations. The compact formulas printed as
[bonderson2007thesis; Eqs. (2.78)--(2.79)](@citet) reverse the ordered
superscripts of the $R$-symbols relative to the defining $R$-move in Eq. (2.54) and
the full hexagon equations. The convention here follows Eq. (2.54),
Eqs. (2.57)--(2.58), and the categorical hexagon.

Equations $\eqref{eq:pentagon-multiplicity-free}$--$\eqref{eq:hexagon-negative-multiplicity-free}$
display the multiplicity-free scalar form. With fusion
multiplicities, the coherence conditions are the same pentagon and hexagon
identities between structural morphisms, using the full matrices and all four
binary-basis indices. The methods `pentagon_axiom(C)` and `hexagon_axiom(C)`
evaluate those morphism equations without a multiplicity-free assumption.

### Changes of fusion bases

A gauge transformation changes the basis of every binary fusion space
$V_{ab}^{e}$. It consequently changes both triple-product bases. If

```math
\label{eq:gauge-basis-change}
L'_u=\sum_s L_sU_{s,u},
\qquad
R'_v=\sum_t R_tV_{t,v},
```

then the associator block in the new bases is

```math
\label{eq:associator-gauge-change}
A'=U^{-1}AV.
```

The matrices $U$ and $V$ are assembled from the binary basis changes, block by
block over the intermediate channels. The same binary basis choices determine
the transformation of $B$. Consequently, raw $F$- and $R$-symbol entries are not
gauge invariants. The corresponding entrywise formulas for changes of binary
splitting bases are
[bonderson2007thesis; Eqs. (2.75)--(2.76)](@citet).

### Relation to other published conventions

The definitions in $\eqref{eq:F-symbol-definition}$ and
$\eqref{eq:R-symbol-definition}$ are the published anyon conventions cited above.
Other sources may instead use projection trees, reverse the associator,
or place the output coordinate first. Those choices change the displayed
matrix without changing the underlying structural morphism.

The fusion spaces in [osborne2019h3; §2, p. 2, and §4, pp. 3--4](@citet) are
$\operatorname{Hom}(d,a\otimes b)$. Its associator map is written in column
coordinates with the output channel as the first matrix index. With
composition-dual bases and matching labels, that matrix is $A^{\mathsf T}$.
This explains the transpose between the package's structural matrices and the
convention used for the $H_3$ formulas in that reference.

A different published convention expresses the $F$-move on projection rather
than splitting trees. The diagrammatic equations in
[ardonne2010clebsch; Eq. (2)](@citet) and
[barter2022associators; Eq. (3)](@citet) suppress associators. Restoring the
parenthesized sources with
$\alpha:(a\otimes b)\otimes c\to a\otimes(b\otimes c)$ gives

```math
\label{eq:inverse-associator-projection-matrix}
L_u\circ\alpha^{-1}=\sum_vM^{F}_{u,v}R_v.
```

Here $u=(e,\mu,\nu)$ and $v=(f,\rho,\sigma)$ have the path order fixed above.
Solving $\eqref{eq:associator-projection-bases}$ for
$L_u\circ\alpha^{-1}$ shows that this matrix is related to the package's
structural matrix by

```math
\label{eq:projection-inverse-conversion}
M^{F}=(A^{-1})^{\mathsf T}.
```

Thus, in the package's row-coordinate realization, the associator block
corresponding to the projection-inverse matrix $M^F$ is
$A=(M^F)^{-\mathsf T}$. The dictionary conventions described above change
the association of keys with entries of $A$; they do not apply the
inverse-transpose operation in
$\eqref{eq:projection-inverse-conversion}$.

There is an analogous distinction for the braiding. The package follows the
direct splitting-space $R$-move in
$\eqref{eq:braiding-splitting-bases}$. If instead one represents the map on
projection spaces induced by the inverse braiding and defines $M^R$ by

```math
\label{eq:inverse-braiding-projection-matrix}
p^{ab}_{d,\mu}\circ c_{a,b}^{-1}
=\sum_\nu M^R_{\mu,\nu}p^{ba}_{d,\nu},
```

then solving $\eqref{eq:braiding-projection-bases}$ gives

```math
\label{eq:projection-inverse-braiding-conversion}
M^R=(B^{-1})^{\mathsf T}.
```

The package stores $B$, not $M^R$. This inverse transpose is independent of
the transpose used by `convention=:column_major_packing`, which only assigns
the entries of $B$ to dictionary keys. The same published source explicitly
distinguishes the braiding from its inverse and defines its $R$-symbols using
the direct move [bonderson2008interferometry; Eqs. (2.30)--(2.32)](@cite).

These comparisons concern mathematical matrix conventions. The dictionary
packing used by the database is a separate software format,
described on the [data exchange page](@ref symbol-data).

## Constructing and checking skeletal data

Conversely, arrays of plausible dimensions do not yet define a fusion
category. The fusion multiplicities must give an associative unital fusion
ring with duals. Every associator block must be an invertible matrix of the
prescribed size, the unit blocks must have the [normalization fixed
above](@ref unit-normalization), and all pentagon equations must hold. Braiding requires invertible blocks
of the prescribed sizes satisfying both hexagon equations. Pivotal and
spherical structures require their own coherence conditions.

The call `six_j_category(K,N,names)` creates a mutable container with fusion
array $N$, identity associator blocks of the required sizes, and all-one
pivotal components. The names argument may be omitted, in which case the
labels are `X1`, `X2`, and so on. The shorter call
`six_j_category(K,names)` sets only the coefficient ring, rank, labels, and
all-one pivotal components; a subsequent
`set_tensor_product!` call installs the fusion array and initializes the
identity associator blocks. Neither form sets the tensor unit or certifies any
coherence axiom. After entering data, use `pentagon_axiom(C)`,
`hexagon_axiom(C)`, and the relevant checked structural predicates.

These initializer methods accept a Julia `Ring`, but the split fusion-category
interpretation on this page and algorithms that use dimensions of Hom spaces
require $K$ to be a field. The initializer does not check that $N$ is a
nonnegative, associative, unital fusion table or that the number of supplied
names matches its rank.

The [worked examples](@ref working-with-fusion-data) now apply these
definitions to explicit categories. The [data exchange](@ref symbol-data)
section specifies the exact key packing and serialization metadata.
