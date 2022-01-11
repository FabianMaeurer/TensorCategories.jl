# Representations

We provide a simple abstract type hirarchy for representation categories:

```
abstract type RepresentationCategory{T} <: TensorCategory{T}
abstract type GroupRepresentationCategory{T,G} <: RepresentationCategory{T}
```
