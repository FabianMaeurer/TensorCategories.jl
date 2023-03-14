```@setup VS
using TensorCategories, Oscar
```


# Vector Space Categories

Vector spaces in TensorCategories are of the abstract type

```
abstract type VectorSpaceCategoryObject{T} <: CategoryObject end
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
VectorSpaces <: TensorCategory
```

and can be constructed like so:

```@example VS
F,a = FiniteField(5,2)
Vec = VectorSpaces(F)
```

CategoryObjects of this category are of the type

```
VSCategoryObject <: VectorSpaceCategoryObject
```

Every vector space object is defined by a basis and a base field provided by the
parent category.

```@docs
VectorSpaceCategoryObject
VectorSpaceCategoryObject(::VectorSpaces,::Int)
```

Morphisms in this Category are defined only by matrices of matching dimensions.
They are typed as

```julia
VSCategoryMorphism <: CategoryMorphism
```

and constructed giving a domain, codomain and matrix element.

```@docs VS
CategoryMorphism(::VectorSpaceCategoryObject, ::VectorSpaceCategoryObject, ::MatElem)
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
F,a = FiniteField(2,3)
VecG = GradedVectorSpaces(F,G)
```

To add a nontrivial associator (twist) construct a Cocycle{3} object coding a 3-cocycle
of the group ``G``. By now no checking of this condition happens.

```@example VS
function cyclic_group_3cocycle(G, F, ξ)
	g = G[1]
	n = order(G)
	D = Dict((g^i,g^j,g^k) => ξ^(div(i*(j+k - rem(j+k,n)),n)) for i ∈ 0:n-1, j ∈ 0:n-1, k ∈ 0:n-1)
	return Cocycle(G,D)
end

F,ξ = CyclotomicField(5, "ξ")
c = cyclic_group_3cocycle(G,F,ξ)
VecG = GradedVectorSpaces(F,G,c)
```

Graded vector spaces decompose into direct sums of vector spaces for each element in
``G``.

```
GVSCategoryObject <: VectorSpaceCategoryObject
```

```@example VS
G = symmetric_group(5)
g,s = gens(G)
V1 = VectorSpaceCategoryObject(QQ,5)
V2 = VectorSpaceCategoryObject(QQ, [:v, :w])
W = VectorSpaceCategoryObject(g => V1, s => V2, g*s => V1⊗V2)
```

Morphisms are implemented analogously by pairs of group elements and vector space objects.

```
GVSCategoryMorphism <: CategoryMorphism
```

The constructor is given by 

```@docs
CategoryMorphism(::GVSCategoryObject, ::GVSCategoryObject,::MatElem) where {G, S <: VectorSpaceCategoryMorphism}
```


## Functionality

(Graded) vector spaces form a fusion category. Thus the methods for
direct sums, tensor products, dual, one and zero object are all implemented.

```@autodocs
Modules = [TensorCategories]
Pages = ["VectorSpaces.jl"]
Order = [:function]
```
