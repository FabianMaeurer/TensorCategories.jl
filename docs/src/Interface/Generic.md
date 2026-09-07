# [Implementation checklist and generic methods](@id interface-checklist)

Start with the structures your model actually supports. A property predicate is
a declaration used by dispatch and algorithms; it is not a request to construct
missing operations.

| Structure | Operations to provide or obtain from valid fallbacks |
|:--|:--|
| Category | `parent`, `domain`, `codomain`, object/morphism equality, `id`, `compose` |
| Linear category | `base_ring`, morphism addition and scalar multiplication, `zero_morphism`, a finite basis for `Hom` |
| Additive category | `zero(C)`, binary `direct_sum` with inclusions and projections |
| Abelian category | `kernel` with inclusion, `cokernel` with projection |
| Monoidal category | tensor product on objects **and** morphisms, `one(C)`, `associator` and its inverse |
| Rigid category | chosen left duality through `dual`, `ev`, `coev`, and chosen right duality through `right_dual`, `right_ev`, `right_coev`; the generic right-duality methods require `pivotal` |
| Semisimple category | effective `decompose` and simple representatives where enumeration is finite |
| Split fusion category | finite split semisimple tensor structure, simple unit; coherent associators and dualities |
| Pivotal/spherical category | `pivotal` and the required monoidal coherence; equality of left and right traces for spherical structure |
| Braided category | `braiding` and the two hexagon identities |

## Basic methods and generic fallbacks

The default `parent(X::Object)` method reads a field named `parent`.
Similarly, the default `domain(f)` and `codomain(f)` methods read fields of
those names. The fallback `base_ring(C)` first reads a field named `base_ring`;
if none exists, it follows a field named `category` and asks that underlying
category for its coefficient ring. A model with a different representation
must provide the corresponding methods explicitly.

For the linear structure, provide morphism addition, scalar multiplication,
zero maps, `Hom`, and an effective basis. `HomSpace(X,Y,B)` wraps a supplied
basis `B`. A custom Hom-space type can instead subtype `AbstractHomSpace` and
implement `domain`, `codomain`, `basis`, and `base_ring`.

The generic `express_in_basis` method uses `matrix(f)` and linear algebra. A
model without faithful matrix coordinates must provide its own coordinate
method.

An operation can be primitive in one model and derived in another. For example,
`image(f)` can be computed as the kernel of the cokernel, while a representation
model can compute it directly by linear algebra. Generic coordinate routines
typically need faithful matrices and usable Hom bases. Having objects with
numeric fields is not sufficient.

## Return values are part of the interface

`kernel(f)` and `cokernel(f)` return pairs, not just objects.
`direct_sum(X,Y)` returns an object and two lists of structural maps.
When supported, `is_isomorphic(X,Y)` returns a Boolean and a morphism witness;
use its first entry as a condition. A backend may instead throw when it cannot
decide isomorphism over the chosen coefficient field. See the
[built-in abelian examples](@ref built-in-abelian-categories).

Do not implement two fallbacks in terms of each other. Test each primitive
before testing operations that depend on it. Read the implementation of a
fallback before relying on its hypotheses; split semisimple coordinate
algorithms are not general algorithms for non-split abelian categories.

## Mathematical checks

Test identities that are independent of the chosen implementation: identity
and associativity of composition, biproduct equations, rank–nullity and kernel
universal properties, tensor interchange, pentagons, duality triangles, and
hexagons. Include zero objects, rectangular matrices, and repeated simple
summands. A nonsymmetric associator and a fusion multiplicity greater than one
detect errors that an Ising-only test cannot.

`pentagon_axiom(C)` and `hexagon_axiom(C)` exhaust the simple inputs. Over an
exact coefficient field, this is a complete coherence check in supported
finite semisimple additive models, where all objects are finite direct sums of
simples and the structural maps extend additively. Over a numerical ball field,
the same functions exhaust the simple tuples but compare the equations at the
chosen working precision; see [Numerical fusion categories](@ref
numerical-fusion-categories).

For API signatures and source links, use the [API reference](../API.md).

Return to the [built-in abelian examples](@ref built-in-abelian-categories), or
consult the [API reference](../API.md) for signatures and source links.
