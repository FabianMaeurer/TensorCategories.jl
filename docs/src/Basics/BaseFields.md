```@meta
DocTestSetup = :(using TensorCategories, Oscar)
```

# [Coefficient fields and numeric computations](@id base-fields)

For a linear category, the coefficient field is part of the mathematical input.
It determines the available scalars, the meaning of equality, and the algebra
algorithms that can be used. TensorCategories.jl calls it `base_ring(C)`, even
when a particular construction requires a field. A Julia method accepting a
value of type `Ring` does not by itself assert that its algorithm is valid over
every ring.

For example, the category of finite-dimensional vector spaces remembers its
coefficient field:

```@example coefficient_fields
using TensorCategories, Oscar
C = vector_spaces(QQ)
X = VectorSpaceObject(C, 2)
@assert base_ring(C) == QQ
@assert base_ring(X) == QQ
base_ring(C)
```

## Exact scalars

Ordinary division of Julia integers produces a floating-point number. Use
OSCAR's rational field `QQ` when exact rational arithmetic and a parent field
are required:

```jldoctest
julia> a = QQ(1)/3;

julia> 3*a == 1
true

julia> parent(a) == QQ
true
```

Other useful exact coefficient fields include prime fields `GF(p)`, number
fields, and the algebraic closure `algebraic_closure(QQ)`. Characteristic is
part of the input, rather than an implementation detail. For example,
representation categories can be semisimple in characteristic zero and
nonsemisimple when the characteristic divides the group order.

## Number fields

A number field is a finite extension of $\mathbb Q$. OSCAR presents it by
generators and polynomial relations. The following constructs
$K=\mathbb Q(s)$ with $s^2=2$:

```@example coefficient_fields
K, s = quadratic_field(2)
@assert s^2 == 2
(K, minpoly(s))
```

The element $s$ is an exact algebraic element. The abstract field does not
declare that it is the positive real square root of $2$; that interpretation
requires a chosen embedding of $K$ into the complex numbers.

Computations are often more efficient over a small number field containing the
required coefficients than over a large algebraic closure. Working over the
smaller field also retains arithmetic information that disappears after
choosing one complex realization. It can, however, prevent objects from
decomposing into absolutely simple summands. The resulting questions of
splitness, scalar extension, embeddings, and Galois conjugacy are treated after
semisimple categories have been introduced.

## Algebraic closures

OSCAR can work exactly over the field of algebraic numbers:

```jldoctest
julia> Qbar = algebraic_closure(QQ)
Algebraic closure of rational field
```

This is an exact algebraically closed field, not a numerical approximation to
$\mathbb C$. Algebraically closed coefficient fields remove division-algebra
phenomena for finite-dimensional endomorphism algebras of simple objects, but
they are not automatically the best computational choice. Support for a
particular package algorithm can also be narrower than the collection of
coefficient fields available in OSCAR.

TensorCategories.jl supports positive-characteristic coefficient fields as a
first-class use case. Field extension within one characteristic can split
objects, but it cannot repair a failure of semisimplicity caused by modular
representation theory.

## [Numerical computations](@id numerical-computations)

A numerical computation begins with a choice of coefficient field, just as an
exact computation does. TensorCategories.jl uses arbitrary-precision ball
arithmetic: the user chooses a working precision, and arithmetic propagates
rigorous enclosures for the represented quantities.

### Arbitrary-precision ball arithmetic

`ArbField(p)` and `AcbField(p)` provide real and complex ball arithmetic at a
user-selected working precision of $p$ bits. A ball records a midpoint and an
error radius. The precision may be chosen as large as needed, but a particular
field, category, and computation use one chosen precision. The implementation
does not automatically increase it unless an algorithm explicitly says so.
See [johansson2017arb](@citet) for the arithmetic model and
[maeurer2026thesis; §5.2.3](@citet) for its use in numerical center and symbol
computations in TensorCategories.jl.

```@example numerical_scalars
using Oscar
R = ArbField(128)
x = sqrt(R(2))
@assert contains_zero(x^2 - R(2))
precision(R)
```

Unlike `Float64`, these fields are not restricted to 53 binary digits and
retain uncertainty as part of every scalar. In `ArbField(p)` and `AcbField(p)`,
the precision belongs to the field and hence to a category over that field.
`ComplexField()` instead uses Nemo's mutable global ball precision. For
reproducible category computations, prefer `AcbField(p)`; conversions produced
by `numeric` use it.

### Exact and numerical models

Exact fields support algebraic equality and symbolic field operations.
Numerical fields support high-precision analytic and linear-algebra
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

### Equality and overlap

Structural equality `==` remains an equivalence relation. It is not replaced by
ball overlap, because overlap is not transitive. Numerical algorithms instead
use explicit tests such as `overlaps(x,y)` and `contains_zero(x)` where the
mathematics asks whether two enclosures are compatible or whether a quantity
may vanish.

### Rigorous enclosures and mathematical conclusions

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
because it lacks an exact symbolic certificate. The pages introducing each
categorical structure state the equations tested by its predicate.

Continue with [simple objects and finite-length categories](@ref simple-objects).
