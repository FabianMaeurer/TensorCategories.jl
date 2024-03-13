```@setup VS
using TensorCategories, Oscar
```


# Vector Space Categories

Vector spaces in TensorCategories are of the abstract type

```
abstract type VectorSpaceObject <: Object end
```

All objects with vector space structure like hom-spaces are and should be implemented as a
subtype of this type. They always need the following fields:

```
basis::Vector
parent::Category
```

## Finite Dimensional VectorSpaces

The simplest example to provide are the finite dimensional vector spaces over a field.
This category has type

```
VectorSpaces <: TensorCategory
```

and can be constructed like so:

```@example VS
F = GF(5,2)
Vec = VectorSpaces()
```

Objects of this category are of the type

```
VSObject <: VectorSpaceObject
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
VSMorphism <: Morphism
```

and constructed giving a domain, codomain and matrix element.

```@docs VS
Morphism(::VectorSpaceObject, ::VectorSpaceObject, ::MatElem)
```

## Graded Vector Spaces

Very similar we have the category of finite dimensional (twisted) ``G``-graded vector spaces for a finite group ``G``.
We have the type

```
GradedVectorSpaces <: VectorSpaces
```
and they are constructed in straightforward manner

```@example VS
G = symmetric_group(6)
VecG = GradedVectorSpaces(G)
```

To add a non-trivial associator (twist) there is another constructor. 

```@docs
TwistedGradedVectorSpaces
```

Graded vector spaces decompose into direct sums of vector spaces for each element in
``G``.

```
GVSObject <: VectorSpaceObject
```

```@example VS
G = symmetric_group(5)
g,s = gens(G)
V1 = VectorSpaceObject(QQ,5)
V2 = VectorSpaceObject(QQ, [:v, :w])
W = VectorSpaceObject(g => V1, s => V2, g*s => V1⊗V2)
```

Morphisms are implemented analogously by pairs of group elements and vector space objects.

```
GVSMorphism <: Morphism
```

The constructor is given by 

```@docs
Morphism(::GVSObject, ::GVSObject,::MatElem) where {G, S <: VectorSpaceMorphism}
```


## Functionality

(Graded) vector spaces form a fusion category. Thus the methods for
direct sums, tensor products, dual, one and zero object are all implemented.

```@autodocs
Modules = [TensorCategories]
Pages = ["VectorSpaces.jl"]
Order = [:function]
```
