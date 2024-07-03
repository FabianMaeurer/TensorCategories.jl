
```@setup Vec
using TensorCategories, Oscar
```

# [The Centralizer Construction](@id centralizer)

Let ``\mathcal C`` be a tensor category and ``\mathcal D \subset \mathcal C`` a full topologizing subcategory. Then the centralizer ``\mathcal Z(\mathcal C \colon \mathcal D)`` is given by tuples ``(Z,\gamma)`` where ``\{\gamma_(X)\colon Z\otimes X \to X \otimes Z \mid X \in \mathcal D}`` satisfies the hexagon equations.

## Example 

Let ``G = S_3``. We want to compute the centralizer of ``\langle \delta_{(1,2)}\rangle``.

```@example Vec
G = symmetric_group(3)

Vec = graded_vector_spaces(QQ,G)

C = centralizer(Vec, Vec[2])

simples(C)
```


