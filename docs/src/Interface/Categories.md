# [Implementing categories](@id category-interface)

TensorCategories.jl represents categories, objects, and morphisms by Julia
values. The abstract types `Category`, `Object`, and `Morphism` provide the
common language used by generic algorithms, while each concrete category model
chooses its own stored data.

## [The interface](@id basic-category-interface)

To implement a category, one must first represent its three kinds of data. A
typical model defines

1. a concrete subtype of `Category` for the category itself;
2. a concrete subtype of `Object` for its objects; and
3. a concrete subtype of `Morphism` for its morphisms.

There may be several object or morphism types in one category. The category
value usually stores data shared by all its objects, such as parameters or a
coefficient field. An object must determine its parent category, and a morphism
must determine its domain and codomain.

The represented types must provide the following behavior:

| Function | Required meaning |
|:---|:---|
| `parent(X)` | the category containing the object $X$ |
| `domain(f)` | the domain of the morphism $f$ |
| `codomain(f)` | the codomain of the morphism $f$ |
| `id(X)` | the identity morphism of $X$ |
| `compose(f,g)` | the composite $g\circ f$, defined when `codomain(f) == domain(g)` |
| `==` | equality of represented categories, objects, and morphisms |

These are functions for which the new model must supply methods. There is one
convenience: the generic implementations of `parent(X)`, `domain(f)`, and
`codomain(f)` read fields named `parent`, `domain`, and `codomain`,
respectively. If the new types use those field names, no additional accessor
methods are needed. A different internal representation must implement the
accessors explicitly.

Constructors should enforce the invariants of the representation. In
particular, a morphism constructor should verify that its stored data define a
map between the stated endpoints and that both endpoints lie in compatible
parent categories.

For composable morphisms

```math
X\xrightarrow{f}Y\xrightarrow{g}Z,
```

TensorCategories.jl uses `compose(f,g)` for $g\circ f$. The infix expression
`g ∘ f` has the usual mathematical order. Defining the methods above does not
by itself prove that the represented data form a category. The implementer is
responsible for ensuring

```math
\label{eq:category-interface-axioms}
\operatorname{id}_Y\circ f=f=f\circ\operatorname{id}_X,
\qquad
h\circ(g\circ f)=(h\circ g)\circ f,
```

with the correct domains and codomains. These identities should be tested on
examples that exercise the representation rather than inferred from the mere
existence of Julia methods.

The operation `==` means equality in the chosen representation. An
implementation must decide when two category values are equal, then compare
objects and morphisms together with their parents and endpoints. Mathematically,
equality is already part of the ambient set- or class-theoretic language in
which a category is defined. Computationally, it must be made effective for
the chosen representation: generic operations use it, for example, to compare
the middle endpoints of two morphisms. It is not an additional categorical
structure, and equality of represented objects is different from isomorphism.

### Extending the interface

The methods above provide only the category structure. Everything else is
added by defining further Julia methods for the represented types. A
categorical construction can be exposed by `product(X,Y)` or `kernel(f)`. A
chosen monoidal structure is exposed by methods such as `tensor_product(X,Y)`,
`one(C)`, and `associator(X,Y,Z)`. Properties and relations are reported by
predicates such as `is_invertible(f)`.

Thus the Julia mechanism is uniform even though the mathematics is not: some
functions return chosen structure, some compute universal constructions, and
some test properties. Each function has a mathematical contract that its
methods must satisfy. The package does not require all these operations to be
stored in the category type or encoded in a single inheritance hierarchy.
The implementer must supply every construction that the model claims to
support. Computing a kernel, for example, may require a substantial
category-specific algorithm; TensorCategories.jl cannot derive one merely
from the basic category interface.

Category-level predicates such as `is_abelian(C)` and `is_finite(C)` are
conservative capability queries. A `true` result means that the implementation
declares or establishes the property, and generic fallbacks include its logical
consequences. A `false` result can mean either that the property fails or that
the implementation has not established it. Thus an algorithm may use a true
result to select a method whose hypotheses apply, but it should not interpret
every false result as a mathematical counterexample.

## [Example: Finite sets](@id finite-sets)

We first implement the category of finite sets. Its objects are finite sets,
its morphisms are functions between them, identities are identity functions,
and composition is ordinary composition of functions. This is the category
described in [EGNO; Example 2.3.1, p. 26](@citet).

TensorCategories.jl already contains a finite-set model constructed by
`Sets()`. Here we independently reimplement the category as a tutorial, using
a Julia `Set` for each object and a complete table of values for each
morphism. The following blocks form one continuous Julia session.

```@example finite_set_tutorial
using TensorCategories

struct FinSetCategory <: Category end

struct FinSetObject <: Object
    parent::FinSetCategory
    elements::Set{Any}
end

FinSetObject(C::FinSetCategory, xs) = FinSetObject(C, Set{Any}(xs))

struct FinSetMorphism <: Morphism
    domain::FinSetObject
    codomain::FinSetObject
    values::Dict{Any,Any}
end
```

The field names already provide the three accessor functions. The morphism
constructor evaluates a Julia function on every domain element and checks that
every result belongs to the codomain:

```@example finite_set_tutorial
function FinSetMorphism(X::FinSetObject, Y::FinSetObject, f::Function)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    values = Dict{Any,Any}()
    for x in X.elements
        y = f(x)
        y in Y.elements || throw(ArgumentError("value outside the codomain"))
        values[x] = y
    end
    FinSetMorphism(X, Y, values)
end
```

Next we define equality and the elementary Julia operations needed to inspect
objects and evaluate morphisms. The names `==`, `length`, `iterate`, and `in`
belong to Julia's module `Base`, whereas `id` and `compose` are exported by
`TensorCategories`. This is why the method definitions below begin with
`Base.:(==)` but `TensorCategories.id` and `TensorCategories.compose`.
Qualifying a definition tells Julia which existing generic function is being
extended. It does not affect how an exported function is called after the
module has been loaded.

```@example finite_set_tutorial
Base.:(==)(::FinSetCategory, ::FinSetCategory) = true
Base.:(==)(X::FinSetObject, Y::FinSetObject) =
    parent(X) == parent(Y) && X.elements == Y.elements
Base.:(==)(f::FinSetMorphism, g::FinSetMorphism) =
    domain(f) == domain(g) && codomain(f) == codomain(g) && f.values == g.values

Base.length(X::FinSetObject) = length(X.elements)
Base.iterate(X::FinSetObject) = iterate(X.elements)
Base.iterate(X::FinSetObject, state) = iterate(X.elements, state)
Base.in(x, X::FinSetObject) = x in X.elements

function (f::FinSetMorphism)(x)
    x in domain(f) || throw(ArgumentError("argument outside the domain"))
    f.values[x]
end
```

Only identity and composition remain for the basic category interface:

```@example finite_set_tutorial
TensorCategories.id(X::FinSetObject) = FinSetMorphism(X, X, x -> x)

function TensorCategories.compose(f::FinSetMorphism, g::FinSetMorphism)
    codomain(f) == domain(g) || throw(ArgumentError("incompatible endpoints"))
    FinSetMorphism(domain(f), codomain(g), x -> g(f(x)))
end
```

We can now construct and compose morphisms through the common interface:

```@example finite_set_tutorial
C = FinSetCategory()
X = FinSetObject(C, [1, 2, 3])
Y = FinSetObject(C, [0, 1])
f = FinSetMorphism(X, Y, x -> x % 2)
g = FinSetMorphism(Y, Y, x -> 1 - x)
h = g ∘ f
@assert domain(h) == X && codomain(h) == Y
@assert h(1) == 0 && h(2) == 1 && h(3) == 0
@assert id(Y) ∘ f == f && f ∘ id(X) == f
nothing # hide
```

### Extending the tutorial category

The code above is a complete implementation of the bare category. We can now
add optional algorithms and constructions. As an example, we implement binary
products. The method returns the Cartesian product together with its two
projections:

```@example finite_set_tutorial
function TensorCategories.product(X::FinSetObject, Y::FinSetObject)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    P = FinSetObject(parent(X), ((x, y) for x in X for y in Y))
    p = [FinSetMorphism(P, X, z -> z[1]), FinSetMorphism(P, Y, z -> z[2])]
    P, p
end

P, p = product(X, Y)
@assert all(p[1](z) == z[1] && p[2](z) == z[2] for z in P)
(length(X), length(Y), length(P))
```

The built-in types `Sets`, `SetObject`, and `SetMorphism` use the same basic
representation and additionally provide products and coproducts:

```@example finite_sets
using TensorCategories
C = Sets()
X = SetObject([1, 2, 3])
Y = SetObject([0, 1])
f = SetMorphism(X, Y, x -> x % 2)
Z, projections = product(X, Y, true)
@assert f(3) == 1
@assert length(Z) == 6
nothing # hide
```

## Example: Positive integers ordered by divisibility

For a second example, regard the positive integers as a category: there is a
unique morphism $m\to n$ precisely when $m$ divides $n$. This is a thin
category, so a morphism needs to store only its endpoints. Given the unique
morphisms $m\to n$ and $n\to r$, their composite is the unique morphism
$m\to r$; it exists because divisibility is transitive. In the implementation,
`DivMorphism(domain(f),codomain(g))` constructs precisely this composite.
Associativity follows because there is at most one morphism between any two
objects.

```@example divisibility_category
using TensorCategories

struct DivisibilityCategory <: Category end

struct DivObject <: Object
    parent::DivisibilityCategory
    n::Int
    function DivObject(C::DivisibilityCategory, n::Int)
        n > 0 || throw(ArgumentError("the integer must be positive"))
        new(C, n)
    end
end

struct DivMorphism <: Morphism
    domain::DivObject
    codomain::DivObject
    function DivMorphism(X::DivObject, Y::DivObject)
        parent(X) == parent(Y) || throw(ArgumentError("different categories"))
        Y.n % X.n == 0 || throw(ArgumentError("the domain does not divide the codomain"))
        new(X, Y)
    end
end

Base.:(==)(::DivisibilityCategory, ::DivisibilityCategory) = true
Base.:(==)(X::DivObject, Y::DivObject) =
    parent(X) == parent(Y) && X.n == Y.n
Base.:(==)(f::DivMorphism, g::DivMorphism) =
    domain(f) == domain(g) && codomain(f) == codomain(g)

TensorCategories.id(X::DivObject) = DivMorphism(X, X)

function TensorCategories.compose(f::DivMorphism, g::DivMorphism)
    codomain(f) == domain(g) || throw(ArgumentError("incompatible endpoints"))
    DivMorphism(domain(f), codomain(g))
end
```

The common interface now composes divisibility relations as morphisms:

```@example divisibility_category
C = DivisibilityCategory()
X, Y, Z, W = (DivObject(C, n) for n in (2, 6, 30, 210))
f = DivMorphism(X, Y)
g = DivMorphism(Y, Z)
h = DivMorphism(Z, W)
@assert domain(g ∘ f) == X && codomain(g ∘ f) == Z
@assert h ∘ (g ∘ f) == (h ∘ g) ∘ f
@assert f ∘ id(X) == f && id(Y) ∘ f == f
(domain(g ∘ f).n, codomain(g ∘ f).n)
```

In this category the categorical product of two objects is their greatest
common divisor, while their coproduct is their least common multiple. Those
constructions could be added with methods for `product` and `coproduct`, just as
the Cartesian product was added to the finite-set model.

Continue with [linear and abelian categories](@ref linear-categories).
