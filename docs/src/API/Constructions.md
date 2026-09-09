# Centers, modules, actions, and functors

The constructors on this page require mathematical structure from their input
categories. Their precise availability therefore depends on the category
model; use Julia's method table together with the linked manual chapter.

| Construction | Principal public names | Manual |
|:---|:---|:---|
| Drinfeld center and half-braidings | `center`, `CenterCategory`, `CenterObject`, `half_braiding`, `is_central`, `center_simples`, `induction`, `induction_restriction` | [Half-braidings and computation](../Constructions/Center.md) |
| Relative center | `centralizer`, `is_relative_braiding` | [Relative centers](../Constructions/Centralizer.md) |
| Algebra objects | `AlgebraObject`, `algebra_structures`, `commutative_algebra_structures`, `separable_algebra_structures`, `etale_algebra_structures`, `is_algebra`, `is_separable` | [Algebras and internal modules](../Constructions/ModuleCategories.md) |
| Internal module categories | `category_of_left_modules`, `category_of_right_modules`, `category_of_bimodules`, `free_left_module`, `free_right_module`, `free_bimodule`, `internal_hom`, `is_left_module`, `is_right_module`, `is_bimodule` | [Algebras and internal modules](../Constructions/ModuleCategories.md) |
| Group actions and equivariantization | `gtensor_action`, `GTensorAction`, `equivariant_induction`, `equivariantization`, `gcrossed_product`, `action_by_inner_autoequivalences`, `is_tensor_action`, `is_equivariant` | [Group actions](../Constructions/GroupActions.md) |
| Functors and natural transformations | `Functor`, `functor`, `NaturalTransformation`, `MonoidalFunctor`, `monoidal_functor`, `monoidal_structure_candidates`, `autoequivalence_candidates` | [Functors](../Interface/Functors.md), [monoidal functors](../Interface/AdvancedInterface.md) |
