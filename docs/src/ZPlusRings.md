# ``\mathbb Z_+``-Rings

Especially for the grothendieck ring of a tensor category the notion of a ``\mathbb Z_+``-Ring is very important. 

A ``\mathbb Z_+``-Ring is a ring that is free as a ``\mathbb Z``-module together with a basis ``\{b_i}_{i \in I}`` such that ``b_ib_j = \sum\limits_{i \in I}c_{ij}^k b_k`` with ``c_{ij}^k \in \mathbb Z_+`` and the unit beein a non-negative linear combination of the basis elements. We call such a ring unital if the unit is a basis element.

In TensorCategories this is realized with the structures

```@docs
ZPlusRing
ZPlusRingElem
```

which extend the structures `Oscar.Ring` and `Oscar.RingElem`respectively. 

```@autodocs
Modules = [TensorCategories]
Pages = ["FusionRings.jl"]
````
## Grothendieck Rings

One way to generate ``\mathbb Z_+``-rings is by computing the Grothendieck ring of a multitensor category.

```@docs
grothendieck_ring
```
