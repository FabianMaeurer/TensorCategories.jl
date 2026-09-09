# Tensor structure and skeletal fusion data

TensorCategories.jl exposes tensor-categorical structure through generic
operations. A category may implement only part of this structure. The manual
explains the required hypotheses before using dimensions, traces, modular
data, or coherence tests.

| Role | Principal public names | Manual |
|:---|:---|:---|
| Tensor product and unit | `tensor_product`, `⊗`, `one`, `associator`, `inv_associator`, `distribute_left`, `distribute_right` | [Tensor products and associators](../Interface/MonoidalCategories.md) |
| Duality | `dual`, `left_dual`, `right_dual`, `ev`, `coev`, `right_ev`, `right_coev` | [Rigidity](../Interface/RigidCategories.md) |
| Traces and dimensions | `left_trace`, `right_trace`, `tr`, `left_dim`, `right_dim`, `dim`, `TensorCategories.squared_norm` | [Pivotal and spherical structures](../Interface/PivotalCategories.md) |
| Pivotal and spherical structure | `pivotal`, `spherical`, `pivotal_structures`, `set_pivotal!`, `set_spherical!`, `is_pivotal`, `is_spherical` | [Pivotal and spherical structures](../Interface/PivotalCategories.md) |
| Braiding and ribbon data | `braiding`, `reverse_braiding`, `twist`, `twists`, `smatrix`, `normalized_smatrix`, `tmatrix` | [Braided categories](../Interface/BraidedCategories.md), [unitarity and modularity](../Interface/UnitaryModularCategories.md), [reversing a braiding](../Interface/BasicConstructions.md) |
| Fusion-category data and predicates | `fusion_coefficient`, `fpdim`, `is_ring`, `is_tensor`, `is_fusion`, `is_multifusion`, `is_unitary`, `is_modular` | [Fusion and multifusion categories](../Interface/TensorCategories.md), [numerical fusion categories](../F-symbols/Numerical.md) |
| Coherence checks | `pentagon_axiom`, `randomized_pentagon_axiom`, `hexagon_axiom`, `monoidal_functor_axiom` | [Monoidal categories](../Interface/MonoidalCategories.md), [braided categories](../Interface/BraidedCategories.md), [skeletal models](../F-symbols/SkeletalFusion.md) |
| Skeletal model and structure setters | `SixJCategory`, `six_j_category`, `skeletonize`, `set_tensor_product!`, `set_one!`, `set_associator!`, `set_braiding!`, `set_pivotal!` | [Skeletal models](../F-symbols/SkeletalFusion.md), [working with fusion data](../F-symbols/WorkedExamples.md) |
| Symbol-dictionary extraction | `F_symbols`, `R_symbols`, `P_symbols` | [Precise conventions](../F-symbols/SkeletalFusion.md#f-conventions), [data exchange](../F-symbols/Data.md) |
| Associator-block generation | `six_j_symbols` | [Skeletal models](../F-symbols/SkeletalFusion.md), [matrix coordinates](../Interface/LinearCategories.md#matrix-realizations) |
| Numerical conversion and symbols | `numeric`, `numeric_F_symbols`, `numeric_R_symbols`, `numeric_P_symbols`, `numeric_smatrix`, `numeric_twists` | [Numerical fusion categories](../F-symbols/Numerical.md), [numerical computations](../Basics/BaseFields.md#numerical-computations) |
| Tensor powers and generated additive closures | `tensor_power`, `tensor_power_category` | [Products and related constructions](../Interface/BasicConstructions.md) |
| Skeletal fusion subcategories | `fusion_subcategory`, `simple_fusion_subcategories` | [Products and related constructions](../Interface/BasicConstructions.md), [skeletal models](../F-symbols/SkeletalFusion.md) |

The structure setters describe a skeletal presentation; a collection of arrays
becomes valid fusion-category data only when the relevant axioms and
nondegeneracy conditions hold. See the linked chapters before constructing or
interpreting $F$- and $R$-symbols.
