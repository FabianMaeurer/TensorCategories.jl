# Basic Categorical Constructions

Let ``\mathcal C``and ``\mathcal D`` be categories.

## Opposite Category

The __opposite category__ ``\mathcal C^{op}`` of ``\mathcal C`` has the same objects and 
morphisms ``f\colon  Y \to X \in \mathrm{Hom}_{\mathcal C^{op}}(X,Y)`` formally switching domain and codomain.

```@docs 
opposite_category
opposite_object
opposite_morphism
```

## Product Categories

Given a family of categories ``\mathcal C_1,...,\mathcal C_n`` we can form the product category ``\mathcal C = \mathcal C_1 \times \cdots \times \mathcal C_n``. Objects and morphism are just families of objects and morphisms. Structures all categories have and are preserved 
by the product will be available.

```@docs
product_category
product_object
product_morphism
```