# Vector Space Categories

Vector spaces in JuCat are of the abstract type

```
VectorSpaceObject{T} <: Object
```

They always need the following fields:

```
basis::Vector{Any}
parent::Category
```

## The Category of Finite Dimensional VectorSpaces

The simplest example to provide are the finite dimensional vector spaces over a field.
This category has type

```
VectorSpaces{T} <: TensorCategory{T}
```

and can be constructed like so:

```@example
using Oscar #hide
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
