# Groethendieck Ring

TensorCategories can compute the Groethendieck ring of a semisimple multitensor category,
if provided with methods

```
# Return the simple objects
simples(C::YourCategory) ::Vector{YourObject}

# Decompose into simple summands with multiplicity
decompose(X::YourObject) ::Vector{Tuple{YourObject, Int}}

# Check for isomorphy
is_isomorphic(X::YourObject, Y::YourObject) ::Tuple{Bool, Union{YourMorphism, Nothing}}
```

Provided these methods exist you can call the function `grothendieck_ring` on your
category.

```@docs
grothendieck_ring
```

Since this is a generic methods it is by nature not the fastest. Thus if you know more
about your category you should specify `grothendieck_ring` on your type manually.

## Example

Take a look at the convolution category:

```@setup Ex
using TensorCategories, Oscar
```


```@example Ex
G = symmetric_group(2)
X = gset(G, [1,2,3])
F = GF(5)
Conv = ConvolutionCategory(X,F)

R,f = grothendieck_ring(Conv)
```

# The S-Matrix

For a braided fusion category with spherical structure and simple objects ``X_i``
the ``S``-matrix is defined as

```math
\begin{aligned}
(S)_{ij} = Tr(c_{X_iX_j} \circ c_{X_jX_i})
\end{aligned}
```

Hence you see you need to provide a trace function and a braiding. The braiding

```
braiding(X::YourObject, Y::YourObject) ::YourMorphism
```

There are two ways to archive the trace: either you provide the trace by yourself

```
tr(f::YourMorphism) ::YourMorphism
```

or you provide a spherical structure, duals, evaluation and coevaluation

```
spherical(X::YourObject) ::YourMorphism
dual(X::YourObject) ::YourObject
ev(X::YourObject) ::YourMorphism
coev(X::YourObject) ::YourMorphism
```

with which the generic left categorical trace will be computed.

```@docs
tr
```

## Example

We can compute now the ``S``-matrix of for the representations.

```@example Ex
G = symmetric_group(3)
F = GF(5)
Rep = representation_category(F,G)

smatrix(Rep)
```
