```@meta
CurrentModule = TensorCategories
```

# Representations

We provide a simple abstract type hierarchy for representation categories:

```
abstract type RepresentationCategory <:Category
```

A representation category always requires the field

```
base_ring::Field
```

## Representations of Finite groups

Let ``G`` be a finite group. We consider the category of finite dimensional
``k``-representations of ``G``.

```
GroupRepresentationCategory <: RepresentationCategory
```

Build it with the constructor


```@docs
RepresentationCategory(::GAPGroup, ::Field)
```

A group representation is defined by a group homomorphism from ``G`` into a
finite dimensional vector space ``k^n``. These objects are of type

```
GroupRepresentationCategoryObject <: RepresentationCategoryObject
```

They are constructed in one of two ways, either by imag"es of generators or by a function

```@docs
Representation(::GAPGroup,::Vector,::Vector)
Representation(::GAPGroup,::Function)
```

where in both cases the images are required to be fitting MatrixElem objects.

Since group representation categories are tensor categories, we again have methods
for the important operations

```@autodocs
Modules = [TensorCategories]
Pages = ["GroupRepresentations.jl"]
```
