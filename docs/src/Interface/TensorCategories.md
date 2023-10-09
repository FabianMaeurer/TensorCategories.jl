# Tensor Categories

Let $k$ be a field. Then a __multitensor category__ is a category $\mathcal C$ which is

- locally finite
- k-linear
- abelian
- rigid monoidal

such that

- ``\otimes`` is bilinear on morphisms.

We call a multitensor category 

- __tensor__ if additionally ``\mathrm{End}(\mathbb 1) = k``,
- __multifusion__ if semisimple and finite and
- __fusion__ if tensor, semisimple and finite.

Thus to implement a one of the above one simply has to implement the interfaces which are part of the definition.

## Multifusion Categories

Semisimple k-linear abelian categories have such a structure that allows one to describe them up to equivalence by matrices. Let ``\{X_i \mid i \in \mathcal I\}`` be the set of non-isomorphic simple objects in ``\mathcal C``. This we can establish an equivalence 

```math
\mathcal C \cong \bigoplus\limits_{i ∈ \mathcal I} \mathrm{Vec}_k
```

of abelian categories via ``X \mapsto \bigoplus \mathrm{Hom(X_i, X)}``. For the detailed construction we refer to [TAMBARA1998692](@cite).

This basically implies that morphisms are merely matrices allowing us to efficiently compute thins like subobjects etc. Thus we want to choose *any* faithful functor 

```math
F \colon \mathcal C \to \mathrm{Vec}_k
```

to provide a method

- `matrix(f::YourMorphism)::MatElem`.

This will open acces to the following operations which are necessary for many computations with fusion categories like for example the computation of the [categorical center](@ref center).

```@docs
eigenvalues
simple_subobjects
express_in_basis
left_inverse
right_inverse
```