# Linear Categories

Let ``k`` be any field. A category is called __``k``-linear__ whenever it is enriched over the category of ``k``-vector spaces, i.e. all Hom-spaces are ``k``-vector spaces and composition is ``k``-linear. We need the following method:

- `*(λ, f::YourMorphism)::YourMorphism` returning the multiplication if a scalar λ.

## Rational Forms

In the literature most categories are defined over an algebraically closed field of characteristic 0 or even over ``\mathbb C``. This is technically possible to implement utilizing the implementation of algebraic numbers in [Nemo.jl](http://nemocas.github.io/Nemo.jl/dev/algebraic/).

In general it is very interesting to work with categories not defined over algebraically closed fields. Especially it might be of interest to implement a category that is usually defined over ``\mathbb C`` over a finite extension of ``\mathbb Q``. 

Let ``\mathcal C`` be a ``K``-linear category and ``k \subset K`` a field extension. Then a category ``\tilde \mathcal C``is called a __rational form__ for ``\mathcal C`` if the karoubian envelope of ``\tilde \mathcal C \otimes K``is equivalent to ``\mathcal C``. We call the rational form __complete__ if already ``\tilde \mathcal C \otimes K`` is equivalent to ``\mathcal C``.

Unfortunately the notion of a complete rational form violates the [principle of equivalence](https://ncatlab.org/nlab/show/principle+of+equivalence). For example the center construction does not preserve the completeness of the rational form. 

