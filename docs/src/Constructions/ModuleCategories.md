# Internal Module Categories

Let ``\mathcal C`` be a fusion category. Any finite module category over ``\mathcal C``can be realized as an internal module category ``\mathrm{Mod}_A(\mathcal C)`` for an algebra ``A``in ``\mathcal C``. 

## Finding Algebras

There are four kinds of algebras of interest. Let ``(A,m,u)`` be an algebra.

- Algebra objects
- Separable algebra objects, i.e. ``m\colon A\otimes A \to A``splits as a bimodule morphisms

And if ``\mathcal C``admits a braiding ``c_{-,-}``

- Commutative algebras, i.e. ``m ∘ c_{A,A} = m``
- Etale algebras, i.e. separable, commutative algebras.

We can find those structures by setting up a system of quadratic equations. Those systems are often of dimension greater then zero and hence we have to guess some solutions. 

```@docs
algebra_structures
separable_algebra_structures
commutative_algebra_structures
etale_algebra_structures
```

## Internal Module Categories

When obtained an algebra we can set up the left, right and bimodule categories. For compatible modules also the tensor product over ``A`` is available.

```@docs
category_of_left_modules
category_of_right_modules
category_of_bimodules
```

Somtimes it might be handy to construct some free modules by hand:

```@docs
free_right_module
free_left_module
free_bimodule
free_module
```

And also conversions from algebras and bimodules:

```@docs
right_module
left_module
bimodule
``` 
