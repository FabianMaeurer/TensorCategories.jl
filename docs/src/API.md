# API reference

This appendix is a map of the principal public operations in
TensorCategories.jl. It is deliberately curated: the package uses Julia's
multiple dispatch, so an exported generic function need not have a method for
every category model. The mathematical chapters state the hypotheses,
conventions, and available algorithms for each model.
Some interface operations extend OSCAR generics; the tables qualify a name
when TensorCategories.jl does not export it itself.

Use Julia's help mode and method table for the exact signatures supported by
the installed version:

```julia-repl
help?> center

julia> methods(center)
```

The reference is organized by mathematical role.

| Area | Contents |
|:---|:---|
| [Category interface](API/Framework.md) | Categories, objects, morphisms, Hom spaces, additive and abelian operations |
| [Tensor structure](API/TensorStructure.md) | Tensor products, duality, braidings, coherence checks, and skeletal fusion data |
| [Concrete categories](API/Categories.md) | Realizations by sets, vector spaces, gradings, representations, and sheaves |
| [Fusion data and databases](API/FusionData.md) | Supplied skeletal categories and data loaders |
| [Constructions](API/Constructions.md) | Centers, internal modules, group actions, and functors |
| [Utilities](API/Utilities.md) | Grothendieck rings, numerical conversion, and persistence |

For a linear introduction to the package, begin with [models and the category
interface](Interface/Categories.md). The [catalogue](F-symbols/Examples.md)
records the coefficient fields, conventions, limitations, and references for
the supplied category families.
