# Skeletal Fusion Categories 

In many applications of fusion categories in mathematical Physics they are assumed to be _skeletal_, i.e. each isomorphism class containes only one object. Then fusion categories reduce to combinatorial objects in the sense that all objects are defined by the multiplicities of their simple subobjects. Morphisms reduce to families of matrices. The associativity constraints are given by the ``F``-symbols.

## ``F``-Symbols 

We use the following conventions for the ``F``-symbols. Let ``\mathcal c`` be a fusion category with simples ``\mathcal O(\mathcal C)``. For all ``a,b,c \in \mathcal O(\mathcal C)`` fix bases for the spaces ``H_{a,b}^c = \mathrm{Hom}(a \otimes b, c)``. Then for all ``a,b,c,d \in \mathcal O(\mathcal C)`` there are canonical bases 

```math 
\begin{align}
    \mathrm{Hom}((a \otimes b) \otimes c, d) = \langle \beta \circ (\alpha \otimes \mathrm{id}_c) \mid \alpha \in H_{a,b}^e, \beta \in H_{e,c}^d, e \in \mathcal O(\mathcal C) \rangle 
\end{align}
```

and 

```math
\begin{align}
    \mathrm{Hom}(a \otimes (b \otimes c), d) = \langle \delta \circ (\mathrm{id}_a \otimes \gamma) \mid \gamma \in H_{b,e}^f, \delta \in H_{a,f}^d, f \in \mathcal O(\mathcal C) \rangle.
\end{align}
```

Now we can express the implied map of the associator ``a_{a,b,c}`` on these spaces in the canonical bases and the resulting matrices 

```math 
\begin{align}
    \left[F_{a,b,c}^d\right]_{(f,\gamma,\delta)}^{(e,\beta,\alpha)}
\end{align}
```
are known as the _``F``-symbols_ of ``\mathcal C``.

## In TensorCategories.jl

We provide a structure for fusion categories defined by fusion rules and ``F``-symbols. The type is `SixJCategory` and can be used with the following constructor.

```@docs; canonical = false 
six_j_category
````

Moreover it is necessary to set some structures like associators, the unit and pivotal structure (if desired).

```@docs
set_associator!
set_one!
set_pivotal!
set_tensor_product!
```

