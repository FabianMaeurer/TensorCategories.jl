# Multitensor Categories

The main idea of TensorCategories is to provide abstract methods for categorical computations in (multi)tensor, respectively (multi)fusion categories. To use this framework with custom structs implementing such a category one has to follow some rules.

# Necessities

For basic functionality of `YourCategory <: Category` together with Objects
`YourObject <: Object` and morphisms `YourMorphism <: Morphism` you have to minimally build

```julia
struct YourCategory <: Category
  base_ring::Field
end

struct YourObject <: Object
  parent::YourCategory
end

struct YourMorphism <: Morphism
  domain::YourObject
  codomain::YourObject
end
```

For objects you need to provide the following methods:

```
dsum(X::YourObject, Y::YourObject) ::YourObject
tensor_product(X::YourObject, Y::YourObject) ::YourObject
id(X::YourObject) ::YourMorphism
```

For morphisms you need

```
dsum(f::YourMorphism, g::YourMorphism) ::YourMorphism
tensor_product(f::YourMorphism, g::YourMorphism) ::YourMorphism
compose(f::YourMorphism, g::YourMorphism) ::YourMorphism
```

Last you should give constructors for the one and zero

```
one(C::YourCategory) ::YourObject
zero(C::YourCategory) ::YourObject
```

# Groethendieck Ring

TensorCategories can compute the Groethendieck ring of a semisimple multitensor category,
if provided with methods

```
# Return the simple objects
simples(C::YourCategory) ::Vector{YourObject}

# Decompose into simple summands with multiplicity
decompose(X::YourObject) ::Vector{Tuple{YourObject, Int}}

# Check for isomorphy
isisomorphic(X::YourObject, Y::YourObject) ::Tuple{Bool, Union{YourMorphism, Nothing}}
```

Provided these methods exist you can call the function `grothendieck_ring` on your
category.

```@docs
grothendieck_ring
```

Since this is a generic methods it is by nature not the fastest. Thus if you know more
about your category you should specify `grothendieck_ring` on your type manually.

## Example

Take a look at the convolution category:

```@setup Ex
using TensorCategories, Oscar
```


```@example Ex
G = symmetric_group(2)
X = gset(G, [1,2,3])
F,_ = FiniteField(5)
Conv = ConvolutionCategory(X,F)

R,f = grothendieck_ring(Conv)
```

# The S-Matrix

For a braided fusion category with spherical structure and simple objects ``X_i``
the ``S``-matrix is defined as

```math
\begin{aligned}
(S)_{ij} = Tr(c_{X_iX_j} \circ c_{X_jX_i})
\end{aligned}
```

Hence you see you need to provide a trace function and a braiding. The braiding

```
braiding(X::YourObject, Y::YourObject) ::YourMorphism
```

There are two ways to archieve the trace: either you provide the trace by yourself

```
tr(f::YourMorphism) ::YourMorphism
```

or you provide a spherical structure, duals, evaluation and coevaluation

```
spherical(X::YourObject) ::YourMorphism
dual(X::YourObject) ::YourObject
ev(X::YourObject) ::YourMorphism
coev(X::YourObject) ::YourMorphism
```

with which the generic left categorical trace will be computed.

```@docs
tr
```

## Example

We can compute now the ``S``-matrix of for the representations.

```@example Ex
G = symmetric_group(3)
F,_ = FiniteField(5)
Rep = RepresentationCategory(G,F)

smatrix(Rep)
```
