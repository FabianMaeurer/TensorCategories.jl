```@meta
CurrentModule = JuCat
```

# Representations

We provide a simple abstract type hirarchy for representation categories:

```
abstract type RepresentationCategory{T} <: TensorCategory{T}
```

A representation category always requires the field

```
base_ring::Field
```

## Representations of Finite groups

Let ``G`` be a finite group. We consider the category of finite dimensinal
``k``-representations of ``G``.

```
GroupRepresentationCategory{T,G} <: RepresentationCategory{T}
```

Build it with the constructor


```@docs
RepresentationCategory(::GAPGroup, ::Field)
```

A group representation is defined by a group homomorphism from ``G`` into a
finite dimensional vector space ``k^n``. These objects are of type

```
GroupRepresentation{T,G} <: Representation
```

They are constructed in one of two ways, either by images of generators or by a function

```@docs
Representation(::GAPGroup,::Vector,::Vector)
Representation(::GAPGroup,::Function)
```

where in both cases the images are required to be fitting MatrixElem objects.
