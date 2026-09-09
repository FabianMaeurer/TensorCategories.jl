# Products, scalar extension, and related constructions

This page lists constructions that transform an existing category or select a
category from it. Each operation requires corresponding methods on
the input model; the existence of the generic function does not imply support
for every category type.

## Products, arrows, and opposites

| Julia construction | Mathematical result |
|:--|:--|
| `ArrowCategory(C)` | The category whose objects are arrows of $\mathcal C$ and whose morphisms are commutative squares |
| `opposite_category(C)` | The opposite category $\mathcal C^{\mathrm{op}}$; use `opposite_object` and `opposite_morphism` for its elements |
| `product_category(C,D)` | The categorical product, with pairs of objects and morphisms |
| `C ⊠ D` | The Deligne tensor product of multifusion `SixJCategory` models over coefficient fields admitting a common coercion |

The categorical product and Deligne tensor product are different. For finite
split semisimple categories, simple objects of the Deligne product are pairs of
simple objects, while simple objects of the additive categorical product are
supported in one factor.

The opposite category reverses every arrow. Reversing a braiding is a different
operation: `reverse_braiding(C)` leaves the category and tensor product in
place and replaces
```math
\label{eq:reverse-braiding}
c_{X,Y}\quad\text{by}\quad c^{\mathrm{rev}}_{X,Y}=c_{Y,X}^{-1}.
```
The current `reverse_braiding` method supports `SixJCategory` models and
transports their existing monoidal data rather than constructing an opposite
category.

## Coefficients, splitting, and completions

| Julia construction | Mathematical result |
|:--|:--|
| `extension_of_scalars(C,L; embedding=...)` | A model-specific realization of scalar extension along a specified field homomorphism $k\to L$ |
| `split(Z::CenterCategory; absolute=true)` | Search for a splitting extension of a supported Drinfeld center |
| `split(X; max_degree=64, check=false)`, `split(objects; max_degree=64, check=false)` | Over a finite field, split the indecomposable summands of a specified finite family over one finite extension |
| `karoubian_envelope(Z)` | Add images of idempotents in supported center and relative-center models |
| `semisimplify(C)` | Quotient negligible morphisms in supported pivotal models |

Mathematically, one first forms the Hom-space extension
$\mathcal C\otimes_k^{\mathrm{Hom}} L$, with the same objects as $\mathcal C$
and
```math
\label{eq:scalar-extension-hom}
\operatorname{Hom}_{\mathcal C\otimes_k^{\mathrm{Hom}} L}(X,Y)
=\operatorname{Hom}_{\mathcal C}(X,Y)\otimes_k L.
```
New idempotents need not have images in this category. For a Hom-finite
Krull–Schmidt category, the idempotent-complete additive scalar extension is
the Karoubi completion
```math
\label{eq:scalar-extension-karoubi-envelope}
\mathcal C_L
=\operatorname{Kar}(\mathcal C\otimes_k^{\mathrm{Hom}} L).
```
It records the new indecomposable summands which appear after extending the
endomorphism algebras. For a semisimple category over a perfect field, or more
generally when the relevant endomorphism algebras are separable, this is again
semisimple abelian. For a weak fusion category it agrees with the Deligne product
$\mathcal C\boxtimes_k\operatorname{Vec}_L$; see
[etingof2012descent; §3.1](@cite) and
[lopezfranco2013tensor; Theorem 3 and §5](@cite). The
[splitting chapter](@ref splitting-and-scalars) explains the distinction
between the Hom-space extension and its completion in detail.

For a nonsemisimple abelian category, Karoubi completion alone need not
preserve abelianness. The abelian scalar extension, when it exists, is the
Deligne product $\mathcal C\boxtimes_k\operatorname{Vec}_L$; for
$\mathcal C\simeq A\text{-mod}$ this is $(A\otimes_kL)\text{-mod}$. The two
constructions agree when $L/k$ is finite separable.

The function `extension_of_scalars` realizes this construction according to
the category model. For split skeletal input no simple endomorphism algebra
acquires new idempotents, so the implementation transports the coefficient
arrays directly. For a supported `CenterCategory`, it extends the known simple
central objects and decomposes those whose endomorphism algebras split over the
new field. Scalar extension can therefore change the simple objects and their
endomorphism algebras. Specify an embedding whenever the source field has more
than one embedding into the target.

Scalar extension along a field homomorphism preserves characteristic. When the
same formulas are instead reduced to positive characteristic, semisimplicity
and every denominator in the structure maps must be checked again. The
[coefficient-field chapter](@ref base-fields) explains these distinctions, and
the [Ising-center example](@ref ising-center) exhibits splitting after scalar
extension.

For a chosen object or finite list over a finite field, `split` changes only
that family; it does not enumerate the simple objects of the ambient category
or the tensor closure of the chosen family. Its result records the extension
field, embedding, extended objects, and their decompositions.
The default relative-degree bound is `max_degree=64`. In a nonsemisimple
category the construction makes the resulting indecomposable summands
absolutely indecomposable; this does not by itself make them simple. The degree
construction establishes this result, while `check=true` recomputes the residue
endomorphism algebras afterward.

Idempotent completion and semisimplification solve different problems. The
first adds images of idempotents; the second takes a quotient by negligible
morphisms. Neither operation by itself chooses a splitting field.
The current `semisimplify(C)` call constructs a wrapper. Its Hom spaces and
morphism equality use the radical of the categorical trace pairing, so they
require finite Hom bases and scalar-valued traces. Decomposition and simple
enumeration additionally require the corresponding operations from the input
model. The constructor itself does not verify these hypotheses or independently
prove that the resulting wrapper is semisimple.

## Generated subcategories and skeletal coordinates

| Julia construction | Mathematical result |
|:--|:--|
| `tensor_power(X,k)`, `X ⊗ k` | The chosen recursively bracketed tensor power $X^{\otimes k}$ |
| `fusion_subcategory(X)` | The fusion subcategory generated by a `SixJObject` |
| `simple_fusion_subcategories(C)` | The distinct simple fusion subcategories found from simple generators of a skeletal fusion category |
| `tensor_power_category(X...)` | The additive closure of summands found in tensor words in the generators |
| `six_j_category(C)` | A skeleton with chosen matrix coordinates for a supported split semisimple category whose `is_ring(C)` predicate is true |
| `center(C)` | The Drinfeld center, whose objects carry half-braidings |

For $k\geq0$, `tensor_power(X,k)` returns the unit when $k=0$ and otherwise
forms the power with the category's tensor product; a negative exponent is not
defined. This constructs an object and does not construct a subcategory.

`tensor_power_category(X...)` begins with the tensor unit and enlarges its
stored indecomposable list as tensor words are explored. It does not
automatically close under duals, extensions, or subquotients. The more
specialized `fusion_subcategory(X)` uses the finite skeletal fusion rules of
its input. `simple_fusion_subcategories(C)` considers each simple generator,
keeps those generated fusion rings that have no proper fusion subring, and
removes repetitions. Its current implementation consequently requires
skeletal `SixJCategory` input. It does not return every subcategory generated
by a single simple object.

Skeletonization changes the presentation of a supported split semisimple
category; it does not adjoin new half-braidings. Its associator, braiding, and
pivotal coefficients are expressed using one shared ordered system of simple
representatives and multiplicity-space bases. Thus these data describe one
skeletal gauge; the individual pivotal coefficients need not remain unchanged
under another choice of bases. The
[skeletal-model chapter](@ref skeletal-fusion) describes the resulting
coordinates, while the [Drinfeld-center chapter](@ref center) describes the
different construction performed by `center(C)`.

The [API reference](../API.md) lists the principal public names for these
constructions. Use Julia help mode or `methods(name)` for exact signatures.

Continue with [Algebra objects and internal modules](../Constructions/ModuleCategories.md).
