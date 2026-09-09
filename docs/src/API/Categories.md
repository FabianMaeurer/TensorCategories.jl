# Concrete categories

Most models on this page retain the objects and morphisms of a mathematical
realization. They are the natural entry point when those objects and their
linear maps matter, rather than only a finite skeletal table of fusion data.
Principal signatures, defaults, and parameter choices are described on the
linked catalogue pages; use `methods(name)` for the complete installed method
table.

| Family | Principal constructors and types | Catalogue entry |
|:---|:---|:---|
| Finite sets and maps | `Sets`, `SetObject`, `SetMorphism`, `SetHomSet` | [Finite sets](../Interface/Categories.md#finite-sets) |
| Finite-dimensional vector spaces | `vector_spaces`, `VectorSpaces`, `VectorSpaceObject`, `VectorSpaceMorphism`, `VSObject`, `VSHomSpace` | [Vector spaces and gradings](../ConcreteExamples/VectorSpaces.md) |
| Group-graded vector spaces | `graded_vector_spaces`, `twisted_graded_vector_spaces`, `GradedVectorSpaces`, `GVSObject`, `GVSMorphism`, `GVSHomSpace`, `Cocycle`, `cyclic_group_3cocycle`, `unitary_cocycle` | [Vector spaces and gradings](../ConcreteExamples/VectorSpaces.md) |
| Finite-group representations | `representation_category`, `rep`, `GroupRepresentationCategory`, `GroupRepresentation`, `GroupRepresentationMorphism`, `GRHomSpace`, `Representation` | [Group representations](../ConcreteExamples/Representations.md) |
| Equivariant coherent sheaves | `coherent_sheaves`, `convolution_category` | [Equivariant sheaves and convolution](../ConcreteExamples/CoherentSheaves.md) |
| Generic quantum $\mathfrak{sl}_2$ model | `sl2_representations` | [$\mathfrak{sl}_2$, Verlinde, and dihedral models](../ConcreteExamples/UqSl2.md) |

The generic quantum $\mathfrak{sl}_2$ entry is the exception: it is a skeletal
recoupling model with infinitely many simple labels, whose objects are finite
direct sums encoded by sparse multiplicity vectors. It does not represent the
action of $U_q(\mathfrak{sl}_2)$ by matrices.

Constructors for supplied skeletal fusion categories are listed under [fusion
data and databases](FusionData.md). General constructions such as scalar
extension, centers, and internal modules are listed separately because their
availability depends on the input category.
