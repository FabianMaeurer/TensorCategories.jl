```@meta
DocTestSetup = :(using TensorCategories, Oscar)
```

# Julia and OSCAR

The **Basics** chapter explains how mathematical categories are represented and
implemented in TensorCategories.jl. We assume familiarity with elementary
category theory and generally follow the conventions of [EGNO](@citet). No
previous experience with Julia is required; the language features needed in
the manual are introduced as they occur.

!!! note "For readers from mathematical physics"
    Readers primarily interested in anyons or conformal field theory may begin
    with [skeletal fusion categories](@ref skeletal-fusion). They should
    nevertheless read [coefficient fields and numeric computations](@ref
    base-fields), because the choice between exact and numerical scalars is
    part of the category being computed.

[Julia](https://julialang.org/) is a general-purpose programming language
designed especially for numerical, scientific, and technical computing. Its
dynamic, parametric type system organizes all values in one explicit type
hierarchy, and its central programming paradigm is multiple dispatch: a single
function may have different methods selected from the types of all its
arguments. These features make Julia particularly suitable for mathematical
software. Categories, objects, and morphisms can be represented by mathematical
data types, while operations such as composition can retain one name across
many different concrete models. See the Julia manual on
[types](https://docs.julialang.org/en/v1/manual/types/) and
[methods](https://docs.julialang.org/en/v1/manual/methods/) for the full language
description.

## Starting a session

After following the [installation instructions](../index.md#Installation),
start Julia. At the interactive Julia prompt (the REPL), load
TensorCategories.jl and OSCAR with

```julia
using TensorCategories, Oscar
```

Installing a package and loading it are different operations; you need not
install it again each session.

!!! note "First computations"
    Julia compiles code when it is first needed. Loading packages and running a
    computation for the first time can therefore take longer than repeating it.
    Keep the session open while working through the manual.

## Values and types

Julia can be used as a calculator:

```jldoctest
julia> 1 + 1
2

julia> 2^64
0

julia> typeof(2)
Int64

julia> Int
Int64
```

Every Julia value has a type, returned by `typeof`. On a 64-bit system, the
literal `2` has the concrete type `Int64`, and `Int` is an alias for `Int64`.
On a 32-bit system, `Int` instead means `Int32`. Arithmetic on these machine
integer types wraps on overflow, so the second computation evaluates $2^{64}$
in `Int64` arithmetic and returns `0`. Use `BigInt` for integers of unbounded
size:

```jldoctest
julia> BigInt(2)^64
18446744073709551616
```

The types themselves form a hierarchy. `Int64` is a subtype of the abstract
type `Signed`, which is a subtype of `Integer`, then `Real`, `Number`, and
finally `Any`. The operator `<:` tests this relation:

```jldoctest
julia> (Int64 <: Signed, Signed <: Integer, Integer <: Real, Real <: Number, Number <: Any)
(true, true, true, true, true)

julia> Integer <: AbstractFloat
false
```

Concrete types such as `Int64` describe values that can be created and are
final: concrete types cannot have subtypes. Abstract types such as `Integer`
cannot themselves be instantiated; they collect related types under a common
interface. Variables are merely names bound to values and ordinarily need no
type declaration. TensorCategories.jl similarly uses the abstract types
`Category`, `Object`, and `Morphism`, with concrete subtypes for each
implemented model.

The Julia expression `1//2` constructs the exact rational number one half. Its
type is `Rational{Int64}` on a 64-bit system; arbitrary-size numerators and
denominators give values of type `Rational{BigInt}`.

## Multiple dispatch

A function in Julia is a collection of methods. A method signature may restrict
the types of its arguments with `::`, and a call uses the most specific method
applicable to the complete tuple of argument types:

```jldoctest
julia> combine(x::Integer, y::Integer) = "two integers";

julia> combine(x::Integer, y::AbstractString) = "an integer and text";

julia> (combine(2, 3), combine(2, "three"))
("two integers", "an integer and text")
```

Thus `compose(f,g)`, `Hom(X,Y)`, and later `tensor_product(X,Y)` can have
category-specific methods while generic algorithms use the same mathematical
names. A new category is implemented by defining data types for its categories,
objects, and morphisms and adding methods for the operations it supports. In
particular, a binary operation such as composition need not be assigned to one
of its arguments: dispatch can inspect both morphisms.

## Computer algebra

[OSCAR](https://docs.oscar-system.org/stable/) is an open-source computer
algebra system for research in algebra, number theory, geometry, and related
areas. It connects several established mathematical systems through Julia and
supplies the rings, fields, matrices, groups, and algebra algorithms used by
TensorCategories.jl.

For computer algebra, use OSCAR's integers `ZZ` and rational field `QQ` rather
than Julia's `BigInt` and `Rational{BigInt}` values. The OSCAR elements belong
to the actual mathematical parents $\mathbb Z$ and $\mathbb Q$ and participate
in the common ring and field interfaces:

```@example julia
using TensorCategories, Oscar
R, x = polynomial_ring(ZZ, "x")
f = x^2 + 2*x + 1
f^2
show(stdout, MIME"text/plain"(), f^2); println() # hide
```

Here `R, x = ...` assigns two returned values to two variables. The name `x`
denotes an element of `R`, not an unspecified complex number. `ZZ` denotes the
integers and `QQ` the rational field.

The type and the mathematical *parent* answer different questions. The type of
`x` determines which Julia methods can act on its representation, while
`parent(x) == R` records the particular polynomial ring containing it.
Elements of different polynomial rings may have the same Julia type but
different parents. Likewise, an object in TensorCategories.jl has a parent
category, and a morphism has a domain and codomain.

## Reading Julia examples

Indices start at **1**. A vector is written `[a, b, c]`; a matrix is written
`[a b; c d]`. A semicolon separates rows. Use `matrix(K, ...)` to construct
an OSCAR matrix over a specified field:

```@example julia
A = matrix(QQ, [1 2; 3 4])
@assert A[1, 2] == 2
size(A)
```

`@assert condition` checks that `condition` is true and raises an error if it
is not. A successful assertion produces no output. The manual uses assertions
to record the mathematical result expected from an example.

An exclamation mark, as in `sort!`, conventionally indicates mutation. A dot
applies an operation elementwise: `sqrt.([1,4,9])` applies `sqrt` to every
entry. A trailing semicolon suppresses the display of a result. In a function
call, arguments following a semicolon are keyword arguments; for example,
`sort([3,1,2]; rev=true)` requests descending order.

Enter `?` at the REPL to switch to help mode, then type a name to see its
documentation. Type `\circ` followed by Tab to enter `∘`.

Continue with [implementing categories](@ref category-interface).
