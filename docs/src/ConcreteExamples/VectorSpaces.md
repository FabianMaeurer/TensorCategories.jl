```@meta 
DocTestSetup = quote 
    using TensorCategories, Oscar
end
```


# Vector Space Categories

Vector spaces in TensorCategories are of the abstract type

```julia
abstract type VectorSpaceObject <: Object end
```

All objects with vector space structure like hom-spaces are and should be implemented as a
subtype of this type. They always need the following fields:

```julia
basis::Vector
parent::Category
```

## Finite Dimensional VectorSpaces

The simplest example to provide are the finite dimensional vector spaces over a field.
This category has type

```julia
VectorSpaces <: TensorCategory
```

and can be constructed like so:

```jldoctest
F = GF(5,2)
Vec = VectorSpaces(F)

# output
Category of finite dimensional VectorSpaces over Finite field of degree 2 and characteristic 5
```

Objects of this category are of the type

```julia
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

```@docs; canonical = false
morphism(::VectorSpaceObject, ::VectorSpaceObject, ::MatElem)
```

## Graded Vector Spaces

Very similar we have the category of finite dimensional (twisted) ``G``-graded vector spaces for a finite group ``G``.
We have the type

```
GradedVectorSpaces <: VectorSpaces
```
and they are constructed in straightforward manner

```jldoctest; output = false
G = symmetric_group(6)
VecG = graded_vector_spaces(G)

# output
Category of G-graded vector spaces over Rational field where G is Sym(6)
```

To add a non-trivial associator (twist) there is another constructor. 

```@docs
twisted_graded_vector_spaces
```

Graded vector spaces decompose into direct sums of vector spaces for each element in
``G``.

```julia
GVSObject <: VectorSpaceObject
```

```jldoctest; output = false
G = symmetric_group(5)
g,s = gens(G)
V1 = VectorSpaceObject(QQ,5)
V2 = VectorSpaceObject(QQ, [:v, :w])
W = VectorSpaceObject(g => V1, s => V2, g*s => V1⊗V2)

# output
Graded vector space of dimension 17 with grading
PermGroupElem[(1,2,3,4,5), (1,2,3,4,5), (1,2,3,4,5), (1,2,3,4,5), (1,2,3,4,5), (1,2), (1,2), (2,3,4,5), (2,3,4,5), (2,3,4,5), (2,3,4,5), (2,3,4,5), (2,3,4,5), (2,3,4,5), (2,3,4,5), (2,3,4,5), (2,3,4,5)]
```

Morphisms are implemented analogously by pairs of group elements and vector space objects.

```julia
GVSMorphism <: Morphism
```

The constructor is given by 

```@docs
morphism(::GVSObject, ::GVSObject,::MatElem)
```


## Functionality

(Graded) vector spaces form a fusion category. Thus the methods for
direct sums, tensor products, dual, one and zero object are all implemented.

```@autodocs
Modules = [TensorCategories]
Pages = ["VectorSpaces.jl"]
Order = [:function]
```
