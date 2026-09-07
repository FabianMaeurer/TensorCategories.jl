# [Finite and semisimple categories](@id semisimple-categories)

A locally finite $k$-linear abelian category has finite-dimensional Hom
spaces and every object has finite length. It is **finite** if it is equivalent
to the category of finite-dimensional modules over a finite-dimensional
$k$-algebra. Intrinsically, one additionally requires enough projectives and
only finitely many isomorphism classes of simple objects
[EGNO; Definitions 1.8.5--1.8.6](@cite). Here enough projectives means that
every simple object has a projective cover.

The predicate `is_finite(C)` records that the implementation declares this
finiteness property. It is not a test that tries to construct a projective
generator. Some concrete backends currently expose finiteness only through a
stronger fusion-category declaration, so a false result need not prove that the
underlying mathematical category is infinite.

A semisimple category in this manual is an abelian category in which every
object is a finite direct sum of simple objects. This agrees with
[EGNO; Definition 1.5.1](@citet). All semisimple categories considered here
are locally finite over their coefficient field.

The package predicate `is_semisimple(C)` records that this property is known for
the implemented category. When the backend can enumerate the simple objects,
`simples(C)` returns chosen representatives of their isomorphism classes. The
order of this list belongs to the implementation and must not be treated as a
mathematical invariant.

For an object $X$, `decompose(X)` returns pairs `(S,m)` describing an
isomorphism

```math
\label{eq:semisimple-decomposition}
X\cong\bigoplus_{S\in\operatorname{Irr}(\mathcal C)}S^{\oplus m_S}.
```

If $D_S=\operatorname{End}_{\mathcal C}(S)$, then
$\operatorname{Hom}_{\mathcal C}(S,X)$ is a right $D_S$-module by
precomposition and

```math
\label{eq:nonsplit-hom-multiplicity}
m_S=\dim_{D_S}\operatorname{Hom}_{\mathcal C}(S,X).
```

The generic implementation of `decompose(X,S)` uses equation
\eqref{eq:nonsplit-hom-multiplicity}; it does not silently replace the
$D_S$-dimension by the dimension over $k$.

Schur's lemma says that $D_S$ is a division algebra
[EGNO; Lemma 1.5.2](@cite). In a semisimple category, simple and indecomposable
objects coincide. Outside a semisimple category, the converse in Schur's lemma
fails: an object with division endomorphism algebra need not be simple, and an
indecomposable object need not be simple. The generic `is_simple(X)` therefore
uses indecomposability only when the parent category is semisimple; otherwise a
category-specific test is required.

## The representation example

For a finite group $G$, Maschke's theorem determines whether the implemented
representation category is semisimple:

```math
\label{eq:maschke-condition}
\operatorname{Rep}_k(G)\text{ is semisimple}
\quad\Longleftrightarrow\quad
\operatorname{char}(k)\nmid |G|.
```

```@example representation_semisimplicity
using TensorCategories, Oscar
G = cyclic_group(5)
C0 = representation_category(QQ, G)
C5 = representation_category(GF(5), G)
@assert is_semisimple(C0)
@assert !is_semisimple(C5)
(is_semisimple(C0), is_semisimple(C5))
```

Equation \eqref{eq:maschke-condition} concerns semisimplicity only. Even when it
holds, the simple representations need not remain simple after extending the
coefficient field.

When a semisimple category has only finitely many simple isomorphism classes,
it is finite: every object is projective, so the direct sum of representatives
of the simples is a projective generator. This is the finiteness condition used
later in the definitions of fusion and multifusion categories.

Continue with [splitting and scalar extension](@ref splitting-and-scalars).
