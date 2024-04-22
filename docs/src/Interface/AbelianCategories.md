# Interface for Abelian Categories

Abelian categories are all over the place and very important.
Thus we provide an Interface for (pre)additive and abelian categories. First recall the definitions:

A __preadditive category__ is a category $\mathcal C$ such that all Hom-spaces are abelian groups and composition is bilinear. As a consequence all finite products are biproducts, also called direct sums. 

Then $\mathcal C$ is called __additive__ if it is preadditive and all finite products exist.

## Additivity

To implement the preadditive structure you need the following methods

- `direct sum(X::YourObject...)::Tuple{YourObject, Vector{YourMorphism}, Yector{YourMorphism}` returning the direct sum object, the inclusions and projections. You might only implement the binary operation while the package will compile a vector version. This might come with performance issues.
- `+(f::YourMorphism, g::YourMorphism)::YourMorphism` providing the addition on morphisms.
- `zero_morphism(X::YourObject, Y::YourObject)::YourMorphism`

To complete additivity there has to be a zero object.

- `zero(C::YourCategory)::YourObject`
  

## Abelian Categories

A category is called __abelian__ if 

- it is  additive,
- every morphism has a kernel and cokernel,
- every monomorphism and epimorphism is normal.
  
We need the following additional methods:

- `kernel(f::YourMorphism)::Tuple{YourObject, YourMorphism}` providing the kernel tuple $(k,\phi)$ for $f \colon X \to Y$ where $\phi \colon k \hookrightarrow X$ is the embedding.
- `cokernel(f::YourMorphism)::Tuple{YourObject, YourMorphism}` providing the cokernel tuple $(c,\psi)$ for $f \colon X \to Y$ where $\psi \colon Y \twoheadrightarrow c$ is the projection.

## Semisimple Categories

An abelian category is called __semisimple__ if every object decomposes uniquely into a direct sum of simple objects. It might be useful to have the method

- `simples(C::YourCategory)::Vector{YourObject}`


## Categories with fibre functor

Whenever a category ``\mathcal C`` has a fibre functor, i.e. an exact faithful functor ``\mathcal C \to \mathrm{Vec}_k``, we can use matrix calculus to compute technical things we often need to implement certain constructions. Implement an existing fibre functor by providing the a method

```julia
matrix(f::MyMorphism)::MatElem
```


