```@setup VS
using JuCat, Oscar
```

```@meta
CurrentModule = JuCat
```

# Vector Space Categories

Vector spaces in JuCat are of the abstract type

```
abstract type VectorSpaceObject{T} <: Object end
```

All objects with vector space structure like hom-spaces are and should be implemented as a
subtype of this type. They always need the following fields:

```
basis::Vector{Any}
parent::Category
```

## Finite Dimensional VectorSpaces

The simplest example to provide are the finite dimensional vector spaces over a field.
This category has type

```
VectorSpaces{T} <: TensorCategory{T}
```

and can be constructed like so:

```@example VS
F,a = FiniteField(5,2)
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

Morphisms in this Category are defined only by matrices of matching dimensions.
They are typed as

```julia
VSMorphism{T} <: Morphism
```

and constructed giving a domain, codomain and matrix element.

```@docs
Morphism(::VectorSpaceObject{T}, ::VectorSpaceObject{T}, ::MatElem{T}) where T
```

## Graded Vector Spaces

Very similar we have the category of finite dimensional (twisted) ``G``-graded vector spaces for a finite group ``G``.
We have the type

```
GradedVectorSpaces{T,G} <: VectorSpaces{T,G}
```
and they are constructed in straightforward manner

```@example VS
G = symmetric_group(6)
F,a = FiniteField(2,3)
VecG = GradedVectorSpaces(F,G)
```

To add a nontrivial associator (twist) construct a Cocycle{3} object coding a 3-cocycle
of the group ``G``. By now no checking of this condition happens.

```@example
C = #TODO CoCycle(G, )
VecG = #VectorSpaces(G,QQ,C)
```

Graded vector spaces decompose into direct sums of vector spaces for each element in
``G``.

```
GVSObject{T,G} <: VectorSpaceObject{T}
```

```@example VS
G = symmetric_group(5)
g,s = gens(G)
V1 = VectorSpaceObject(QQ,5)
V2 = VectorSpaceObject(QQ, [:v, :w])
W = VectorSpaceObject(g => V1, s => V2, g*s => V1⊗V2)
```

Morphisms are implemented analogously by pairs of group elements and vectorspace objects.

```
GVSMorphism{T,G} <: Morphism
```

The constructors are

```@docs
Morphism(::GVSObject{T,G},::GVSObject{T,G},::Dict{G,S}) where {T,S<:VectorSpaceMorphism,G}
Morphism(::GVSObject{T,G}, ::GVSObject{T,G},::Pair{G,S}...) where {T,G, S <: VectorSpaceMorphism{T}}
```


## Functionality

(Graded) vector spaces form a semisimple tensor category. Thus the methods for
direct sums, standard tensor products, one and zero object are all implemented.

```@autodocs
Modules = [JuCat]
Pages = ["VectorSpaces.jl"]
Order = [:function]
```
