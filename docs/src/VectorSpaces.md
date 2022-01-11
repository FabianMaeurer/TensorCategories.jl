# Vector Space Categories

Vector spaces in JuCat are of the abstract type

```
VectorSpaceObject{T} <: Object
```

All objects with vector space structure like hom-spaces are and should be implemented as a
subtype of this type. They always need the following fields:

```
basis::Vector{Any}
parent::Category
```

## The Category of Finite Dimensional VectorSpaces

The simplest example to provide are the finite dimensional vector spaces over a field.
This category has type

```@doc
VectorSpaces{T} <: TensorCategory{T}
```

and can be constructed like so:

```@example
using Oscar,JuCat # hide
F = FiniteField(5,2)
Vec = VectorSpaces(F)
```

Objects of this category are of the type

```
VSObject{T} <: VectorSpaceObject{T}
```

Every vector space object is defined by a basis and a base field provided by the
parent category.

```@docs
VectorSpaceObject
VectorSpaceObject(::VectorSpaces,::Int)
```

## The Category of Graded Vector Spaces

Very similar we have the category of (twisted) ``G``-graded vector spaces for a finite group ``G``.
We have the type

```@docs
GradedVectorSpaces{T,G} <: VectorSpaces{T,G}
```
and they are constructed in straightforward manner

```@repl
using JuCat, Oscar # hide
G = symmetric_group(6)
F,a = FiniteField(2,3)
VecG = GradedVectorSpaces(G,F)
```

To add a nontrivial associator (twist) construct a Cocycle{3} object coding a 3-cocycle
of the group ``G``. By now no checking of this condition happens.

```@repl
using JuCat, Oscar #hide
C = #TODO CoCycle(G, )
VecG = #VectorSpaces(G,QQ,C)
```

## Functionality

(Graded) vector spaces form a semisimple tensor category. Thus the methods for
direct sums, standard tensor products, one and zero object are all implemented.

```@autodocs
Modules = [JuCat]
pages = ["VectorSpaces.jl"]
Order = [:function]
```
