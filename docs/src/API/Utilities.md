# Rings, numerical conversion, and persistence

These utilities support the main categorical workflows. The list is selective;
use Julia's help mode and method table for exact signatures.

| Role | Principal public names | Manual |
|:---|:---|:---|
| Grothendieck and based rings | `split_grothendieck_ring`, `ZPlusRing`, `multiplication_table`, `print_multiplication_table`, `fpdim` | [Grothendieck rings](../Interface/GrothendieckRings.md), [computing with fusion rings](../Interface/FusionRings.md) |
| Complex and numerical realization | `complex_embedding`, `complex_embeddings`, `numeric` | [Coefficient fields and numeric computations](../Basics/BaseFields.md), [numerical fusion categories](../F-symbols/Numerical.md) |
| Exact fusion-category archives | `save_fusion_category`, `load_fusion_category` | [Data exchange](../F-symbols/Data.md) |
| Numerical symbol tables | `numeric_symbols_to_csv`, `numeric_symbols_from_csv`, `load_numeric_fusion_category` | [Data exchange](../F-symbols/Data.md) |
| Native object persistence | `save`, `load` | [Data exchange](../F-symbols/Data.md) |

Symbol archives and CSV tables carry convention choices that cannot be inferred
from a bare array. Read the [precise symbol conventions](../F-symbols/SkeletalFusion.md#f-conventions)
and the serialization section of the data-exchange chapter before moving data
between programs.
