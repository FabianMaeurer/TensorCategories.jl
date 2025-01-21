```@meta
CurrentModule = TensorCategories
```

# Representations

We provide a simple abstract type hierarchy for representation categories:

```julia
abstract type RepresentationCategory <:Category
```


## Representations of Finite groups

Let ``G`` be a finite group. We consider the category of finite dimensional
``k``-representations of ``G``.

```
GroupRepresentationCategory <: RepresentationCategory
```

Build it with the constructor

```@docs; canonical = false
representation_category(::Field, ::Group)
```

A group representation is defined by a group homomorphism from ``G`` into a
finite dimensional vector space ``k^n``. These objects are of type

```
GroupRepresentationObject <: RepresentationObject
```

They are constructed in one of two ways, either by images of generators or by a function

```@docs; canonical = false
Representation(::Group,::Vector,::Vector)
Representation(::Group,::Function)
```

where in both cases the images are required to be fitting MatrixElem objects.

Since group representation categories are tensor categories, we again have methods
for the important operations

```@autodocs
Modules = [TensorCategories]
Pages = ["GroupRepresentations.jl"]
Private = false
```
