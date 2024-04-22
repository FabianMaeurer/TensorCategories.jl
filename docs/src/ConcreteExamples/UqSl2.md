# Representations of ``\mathfrak{sl}_2(\mathbb k)``

The representation category of ``\mathfrak{sl}_2(\mathbb k)`` has countably infinte
simple objects ``\{V_i \mid i \in \mathbb N_0\}`` where the tensor product is given by the 
Clebsch-Gordon rule

```math
V_i \otimes V_j = ⨁\limits_{l = 0}^{\min{(i,j)}} V_{i+j - 2l}.
```

We can construct the category and specify at any value for ``q`` that is either 1 or not 
a rot of unity.

```@docs; canonical = false
sl2_representations
````

## Verlinde type categories

When specifying at a root of unity we arrive at the Verlinde category. 
These categoryies have ``n`` simple objects ``V_0,...,V_{n-1}``and fusion rule 

```math
V_i \otimes V_j = ⨁\limits_{l = 0}^{\min{(i,j)}} V_{i+j - 2l}.