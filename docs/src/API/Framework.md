# Category interface

The same generic functions are shared by many category implementations. Which
methods are available depends on the structure provided by the category. See
the [interface chapters](../Interface/Categories.md) for the mathematical
meaning and the [implementation checklist](../Interface/Generic.md) when adding
a new model.

| Role | Principal public names | Manual |
|:---|:---|:---|
| Core types, parents, and underlying data | `Category`, `Object`, `Morphism`, `parent`, `category`, `object`, `morphism`, `==` | [Implementing categories](../Interface/Categories.md) |
| Sources, targets, and composition | `domain`, `codomain`, `id`, `compose`, `∘`, `inv`, `is_invertible`, `left_inverse`, `right_inverse`, `composition_power` | [Objects and morphisms](../Interface/Categories.md) |
| Hom spaces and coordinates | `Hom`, `End`, `HomSpace`, `HomSet`, `basis`, `int_dim`, `zero_morphism`, `is_zero`, `matrix`, `matrices`, `express_in_basis` | [Linear categories and matrix coordinates](../Interface/LinearCategories.md) |
| Direct sums and finite (co)products | `direct_sum`, `⊕`, `product`, `×`, `coproduct`, `∐` | [Linear and abelian categories](../Interface/LinearCategories.md) |
| Simple objects and composition factors | `simples`, `simples_names`, `is_simple`, `simple_subobjects`, `composition_factors` | [Simple objects and finite length](../Interface/SimpleObjects.md) |
| Direct-sum decompositions | `decompose`, `indecomposables`, `is_indecomposable`, `Oscar.is_isomorphic` | [Idempotents and Krull--Schmidt categories](../Interface/KaroubianCategories.md), [finite and semisimple categories](../Interface/SemisimpleCategories.md) |
| Kernels, cokernels, and images | `kernel`, `cokernel`, `image`, `is_monomorphism`, `is_epimorphism`, `is_subobject` | [Linear, additive, and abelian categories](../Interface/LinearCategories.md) |
| Structural predicates | `is_linear`, `is_additive`, `is_abelian`, `is_locally_finite`, `is_finite`, `is_krull_schmidt`, `is_semisimple`, `is_split_semisimple` | [Linear categories](../Interface/LinearCategories.md), [simple objects and finite-length categories](../Interface/SimpleObjects.md), [finite and semisimple categories](../Interface/SemisimpleCategories.md), [splitting](../Interface/SplittingFields.md) |
| Standard category constructions | `opposite_category`, `product_category`, `extension_of_scalars`, `split`, `semisimplify`, `karoubian_envelope` | [Products and scalar extension](../Interface/BasicConstructions.md), [coefficient fields](../Basics/BaseFields.md) |

These lists identify the generic functions to look up; they do not assert that
every row applies to every category. In Julia, `methods(name)` gives the
currently installed methods for a name.
