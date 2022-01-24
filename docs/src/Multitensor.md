# Multitensor Categories

The main idea of JuCat is to provide abstract methods for categorical computations.
The abstract types `MultiTensorCategory` and `TensorCategory` are intended for
this purpose. To use these abstract methods you need to provide some functionality
for your categories.

# Necessities

For minimal function of `YourCategory <: MultiTensorCategory` together with Objects
`YourObject <: Object` and morphisms `YourMorphism <: Morphism` you have to minimally build

```julia
struct YourCategory <: MultiTensorCategory
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

For basic functionality on objects you need to provide some methods:

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

JuCat can compute the Groethendieck ring of a semisimple `MultiTensorCategory`,
if provided with methods

```
# Return the simple objects
simples(C::YourCategory) ::Vector{YourObject}

# Decompose into simple summands with multiplicity
decompose(X::YourObject) ::Vector{Tuple{YourObject, Int}}

# Check for isomorphy
isisomorphic(X::YourObject, Y::YourObject) ::Tuple{Bool, Union{YourMorphism, Nothing}}
```

Provided these methods exist you can call the function `groethendieck_ring` on your
category.

```@docs
groethendieck_ring
```

Since this is a generic methods it is by nature not the fastest. Thus if you know more
about your category you should specify `groethendieck_ring` on your type manually.

## Example

Take a look at the convolution category:

```@setup Ex
using JuCat, Oscar
```


```@example Ex
G = symmetric_group(2)
X = gset(G, [1,2,3])
F = FiniteField(5)
Conv = ConvolutionCategory(X,F)

R,f = groethendieck_ring(Conv)
```
