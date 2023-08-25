```@meta
CurrentModule = TensorCategories
```

# Coherent sheaves on Finite Sets

Coherent sheaves on a finite set ``X`` are characterized by the stalks of the
the elements ``x \in X``. I.e. a coherent sheaf on ``X`` is a collection of
vector spaces ``V_x`` for each ``x\in X``. More generally for a group ``G`` and
a G-set ``X`` we can define consider ``G``-equivariant sheaves on ``X``. In this
case a sheaf is specified by representations of the stabilizers of representatives
for the orbits.

## Equivariant Coherent Sheaves

We provide the datatype

```
CoherentSheafObject <: Object
```

The category of equivariant coherent sheaves has type

```
CohSheaves <: MultiTensorCategory
```

and can be constructed via

```@docs
CohSheaves
```

Morphisms are given by morphisms of representations of the stalks and are of type

```
CohSheafMorphism{T,G} <: Morphism
```

```@autodocs
Modules = [TensorCategories]
Pages = ["CoherentSheaves.jl"]
```

## Convolution Category

The objects of this category are again ``G``-equivariant coherent sheaves on a
finite ``G``-set ``X\times X``. But we endow them with a different monoidal product.

Let ``p_{ij}: X\times X\times X \to X \times X`` be the canonical projections.
Then we define the monoidal product of two coherent sheaves ``x`` and ``y``

```math
\begin{aligned}
x\otimes y = p_{13}_\ast(p_{12}^\ast(x)\otimes' p_{23}^\ast(y))
\end{aligned}
```

where ``\otimes'``is the monoidal product of ``Coh(X\times X\times X)``. Similar for
morphisms.

Objects in this category are of type

```
ConvolutionObject <: Object
```

while the convolution category is of type

```
ConvolutionCategory <: MultiTensorCategory
```

and can be constructed by

```@docs
ConvolutionCategory
```

Morphisms are just morphisms of coherent sheaves with the new tensor product.

```
ConvolutionMorphism <: Morphism
```

```@autodocs
Modules = [TensorCategories]
Pages = ["ConvolutionCategory.jl"]
```
