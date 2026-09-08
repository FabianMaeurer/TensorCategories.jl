```@meta
DocTestSetup = :(using TensorCategories, Oscar)
```

# [Coefficient fields and numeric computations](@id base-fields)

For a linear category, the field $k$ is part of the mathematical input, and
`base_ring(C)` returns the corresponding coefficient parent of a category
`C`. For exact computations, this parent may in principle be any field
provided by OSCAR. For numerical computations over $\mathbb R$ or $\mathbb C$,
it may instead be a parent whose elements are rigorous ball enclosures. In
either case, this does not imply that every theorem or algorithm applies: its
hypotheses may restrict the characteristic, require an algebraically closed or
splitting field, or depend on operations that have only been implemented for
certain coefficient types.

Most fields used in TensorCategories.jl are represented exactly. Their
elements are symbolic algebraic objects, and arithmetic and equality are exact.
Numerical computations are also possible using arbitrary-precision real or
complex ball arithmetic, as described below.

## Exact symbolic computations

OSCAR provides several exact fields that commonly occur in tensor-category
computations. For further coefficient fields, constructors, and conversion
functions, see the [Fields chapter of the OSCAR
manual](https://docs.oscar-system.org/stable/Fields/intro/).

| Field | OSCAR constructor | Description |
|:---|:---|:---|
| $\mathbb Q$ | `QQ` | the rational numbers |
| a number field $K/\mathbb Q$ | `number_field(f, "a")`, `quadratic_field(d)` | a finite extension given by algebraic generators and relations |
| an embedded real number field $K\subset\mathbb R$ | `embedded_number_field(f,r)` | a number field together with a chosen real embedding |
| $\mathbb Q^{\mathrm{ab}}$ | `abelian_closure(QQ)` | the maximal abelian extension of $\mathbb Q$, containing all roots of unity |
| $\overline{\mathbb Q}$ | `algebraic_closure(QQ)` | the field of algebraic numbers |
| $\mathbb F_p$ | `GF(p)` | the prime field of characteristic $p$ |
| $\mathbb F_{p^n}$ | `GF(p,n)` | the finite field with $p^n$ elements |
| $\overline{\mathbb F}_p$ | `algebraic_closure(GF(p))` | the algebraic closure of the prime field $\mathbb F_p$ |

Ordinary division of Julia integers produces a floating-point number. Use the
OSCAR field `QQ` when an exact rational number and its mathematical parent
field are required:

```jldoctest
julia> a = QQ(1)/3;

julia> 3*a == 1
true

julia> parent(a) == QQ
true
```

A number field is a finite extension of $\mathbb Q$. OSCAR presents it by
generators and polynomial relations. The following constructs
$K=\mathbb Q(s)$ with $s^2=2$:

```@example coefficient_fields
using TensorCategories, Oscar
K, s = quadratic_field(2)
@assert s^2 == 2
(K, minpoly(s))
```

The element $s$ is an exact algebraic element. The abstract field does not
declare that it is the positive real square root of $2$. An embedding
$\sigma\colon K\hookrightarrow\mathbb C$ specifies which complex root the
generator represents. In this example there are two real embeddings,
$\sigma_+(s)=\sqrt{2}$ and $\sigma_-(s)=-\sqrt{2}$, exchanged by the
nontrivial element of $\operatorname{Gal}(K/\mathbb Q)$. Applying different
embeddings to algebraic coefficients produces their Galois-conjugate
realizations. When an ordered realization inside $\mathbb R$ is required,
OSCAR records the chosen real embedding with `embedded_number_field`.

Computations are often more efficient over a small number field containing the
required coefficients than over a large algebraic closure. Working over the
smaller field also retains arithmetic information that disappears after
choosing one complex realization. It can, however, prevent objects from
decomposing into absolutely simple summands. The resulting questions of
splitness and scalar extension, together with the categorical effect of
embeddings and Galois conjugacy, are treated in
[Scalar extension and splitting](@ref splitting-and-scalars).

The algebraic and abelian closures are exact fields rather than numerical
approximations to $\mathbb C$:

```jldoctest
julia> Qab, z = abelian_closure(QQ);

julia> Qbar = algebraic_closure(QQ);

julia> F25 = GF(5, 2);

julia> order(F25)
25

julia> Fbar = algebraic_closure(GF(5))
Algebraic closure of prime field of characteristic 5
```

Here $\mathbb Q^{\mathrm{ab}}$ is the union of the cyclotomic fields, whereas
$\overline{\mathbb Q}$ contains every algebraic number. Algebraically closed
coefficient fields remove division-algebra phenomena for finite-dimensional
endomorphism algebras of simple objects, but they are not automatically the
best computational choice. The positive-characteristic constructor currently
takes a prime field `GF(p)` and represents
$\overline{\mathbb F}_p$ as the union of its finite extensions.

TensorCategories.jl supports positive-characteristic coefficient fields as a
first-class use case. Field extension within one characteristic can split
objects, but it cannot repair a failure of semisimplicity caused by modular
representation theory.

OSCAR can work exactly with a finite real number field $K\subset\mathbb R$ by
choosing a real embedding with `embedded_number_field`, and every finite
collection of real algebraic numbers is contained in such a field. OSCAR can
also represent real algebraic numbers as elements of $\overline{\mathbb Q}$ and
test whether an element is real. It does not, however, provide the real closed
field $\mathbb R_{\mathrm{alg}}=\overline{\mathbb Q}\cap\mathbb R$ as a
separate exact coefficient field. The parent of a real element represented in
`algebraic_closure(QQ)` is still all of $\overline{\mathbb Q}$, so this does not
model a category whose coefficient field is $\mathbb R_{\mathrm{alg}}$. Exact
computations over a chosen embedded real number field are possible in
principle, but a universal exact real closed coefficient field is not currently
available in OSCAR. Investigating exact and numerical computations with fusion
categories over real fields is ongoing work in TensorCategories.jl.

## [Numerical computations](@id numerical-computations)

A numerical computation begins with an intended coefficient field
$\mathbb R$ or $\mathbb C$ and a computational parent representing its scalars,
just as an exact computation begins with an exact coefficient field.
TensorCategories.jl uses arbitrary-precision ball arithmetic: the user chooses
a working precision, and arithmetic propagates rigorous enclosures for the
represented quantities.

OSCAR provides two interfaces to Arb's real and complex ball arithmetic.
The constructors `real_field()` and `complex_field()` return parents of types
`RealField` and `ComplexField`; their precision is controlled through the
global `Balls` precision. In OSCAR examples, `RR` and `CC` are conventional
variable names assigned with `RR = real_field()` and `CC = complex_field()`.
They are not exact symbolic models of $\mathbb R$ and $\mathbb C$.

These parents belong to AbstractAlgebra's `Field` type hierarchy because they
implement its computational field interface. This is not an assertion that
the ball objects themselves form a field in the algebraic sense: a ball is an
enclosure, and a ball containing zero cannot be inverted. Exactness is a
separate part of the interface. For an exact field such as `QQ`,
`is_exact_type(elem_type(QQ))` is `true`; for `RealField`, `ComplexField`,
`ArbField`, and `AcbField` elements it is `false`. Thus a category with
`base_ring(C) == real_field()` is a numerical model of a category over
$\mathbb R$, computed using rigorous enclosures at the selected precision.

The constructors `ArbField(p)` and `AcbField(p)` instead return real and
complex ball fields whose working precision of $p$ bits is stored in the
parent. TensorCategories.jl currently uses these explicit-precision parents
for its principal numerical workflows, and `numeric(E,p)` produces a category
over an `AcbField`. A ball records a midpoint and an error radius. The
precision may be chosen as large as needed, but a particular field, category,
and computation use one chosen precision. The implementation does not
automatically increase it unless an algorithm explicitly says so. See
[johansson2017arb](@citet) for the arithmetic model.

Unlike `Float64`, ball fields are not restricted to 53 binary digits and retain
uncertainty as part of every scalar. They provide controlled numerical
computation rather than point-valued floating-point arithmetic with untracked
rounding error. In `ArbField(p)` and `AcbField(p)`, the precision belongs to the
field and hence to a category over that field. For reproducible category
computations, prefer these explicit-precision parents.

Here is a scalar computation over a real ball field:

```@example numerical_scalars
using Oscar
R = ArbField(128)
x = sqrt(R(2))
@assert contains_zero(x^2 - R(2))
precision(R)
```

These numerical coefficient parents can also be used directly by a supported
category. For example, the usual coordinate model of vector spaces and its
matrix operations work over a real ball field:

```@example numerical_vector_spaces
using TensorCategories, Oscar
R = ArbField(128)
C = vector_spaces(R)
X = VectorSpaceObject(C, 2)
A = matrix(R, [sqrt(R(2)) 0; 0 R(1)/3])
f = morphism(X, X, A)
g = f ∘ f
@assert base_ring(C) == R
(int_dim(X), contains_zero(matrix(g)[1,1] - R(2)))
```

This returns `(2,true)`: the first diagonal entry of $g$ is a ball containing
the exact value $2$.

!!! warning "Current numerical limitations"
    TensorCategories.jl and its OSCAR/Nemo backends do not yet provide every
    algorithm required by the numerical side of the package. Basic category
    models can often be constructed over a real or complex ball field, and
    specified objects and morphisms can be manipulated using the available
    matrix operations. This does not imply that every higher algorithm is
    implemented for that field type.

    In particular, simple-object enumeration for group representations over
    ball fields is not currently implemented. Direct computation and
    skeletonization of Drinfeld centers over real ball fields also require
    further support for scalar extraction, numerical decomposition and linear
    dependence, and the recovery of fusion multiplicities. Extending these
    interfaces and algorithms is ongoing work.

Exact field types support algebraic equality and symbolic field operations.
Inexact ball types support high-precision analytic and linear-algebra
computations and retain an enclosure for every result. When exact source data
are available, keep them and construct a numerical model for the required
computation. The exact source records the algebraic object; the numerical
model records one realization at one chosen working precision.

For a supported object `E`, `numeric(E,p)` requests a numerical realization at
approximately $p$ bits. A conversion may use guard bits, so the precision of
the resulting base field is authoritative. For coefficients in a number field,
the conversion also requires a complex embedding. Increasing the precision
cannot recover information already lost through decimal or low-precision
input.

Structural equality `==` remains an equivalence relation. It is not replaced by
ball overlap, because overlap is not transitive. Numerical algorithms instead
use explicit tests such as `overlaps(x,y)` and `contains_zero(x)` where the
mathematics asks whether two enclosures are compatible or whether a quantity
may vanish.

Suppose that a ball $B$ contains an exact quantity $b$.

- If $0\notin B$, then $b\ne0$.
- Disjoint balls containing $b$ and $c$ prove $b\ne c$.
- If $0\in B$, or if two balls overlap, this does not prove equality.

Thus a numerical calculation can give an exact mathematical conclusion. For
example, a determinant ball excluding zero certifies nonvanishing. This does
not turn every ball-valued calculation into symbolic algebra: its scalars are
still enclosures rather than exact algebraic expressions.

Categorical predicates over a ball field use the corresponding numerical
criterion at the field's working precision. They do not reject an input merely
because it lacks an exact symbolic certificate. Unless a predicate explicitly
documents a rigorous certificate, a successful result means that the defining
equations hold to the chosen working precision. Ball arithmetic can
nevertheless certify particular conclusions, such as nonvanishing or strict
separation, when the computed enclosures imply them. The pages introducing
each categorical structure state the equations tested by its predicate.

Continue with [simple objects and finite-length categories](@ref simple-objects).
