# [Scalar extension and splitting](@id splitting-and-scalars)

## Scalar extension

Let $\mathcal C$ be a $k$-linear category and let
$\iota\colon k\hookrightarrow L$ be a field embedding. The **Hom-space
extension** $\mathcal C\otimes_k^{\mathrm{Hom}}L$ has the same objects as
$\mathcal C$ and morphism spaces

```math
\label{eq:hom-space-scalar-extension}
\operatorname{Hom}_{\mathcal C\otimes_k^{\mathrm{Hom}}L}(X,Y)
=
\operatorname{Hom}_{\mathcal C}(X,Y)\otimes_{k,\iota}L.
```

Composition is extended $L$-bilinearly. There is a canonical $k$-linear
functor

```math
\label{eq:hom-space-extension-functor}
T^{\mathrm{Hom}}_{L/k}\colon\mathcal C
\longrightarrow\mathcal C\otimes_k^{\mathrm{Hom}}L,
\qquad
X\longmapsto X,
\qquad
f\longmapsto f\otimes1.
```

If $\mathcal C$ is additive, then its Hom-space extension is additive and
$L$-linear. However, when $\mathcal C$ is Karoubian, its Hom-space extension
need not be Karoubian: an algebra
$\operatorname{End}_{\mathcal C}(X)\otimes_kL$ may contain idempotents whose
images are not represented by the old objects. We therefore define the
**scalar extension** of an idempotent-complete additive category by

```math
\label{eq:completed-scalar-extension}
\mathcal C_L
=
\operatorname{Kar}\!\left(\mathcal C\otimes_k^{\mathrm{Hom}}L\right).
```

An object of $\mathcal C_L$ is a pair $(X,e)$ with
$e^2=e\in\operatorname{End}_{\mathcal C}(X)\otimes_kL$, and

```math
\label{eq:completed-scalar-extension-hom}
\operatorname{Hom}_{\mathcal C_L}\bigl((X,e),(Y,d)\bigr)
=
d\bigl(\operatorname{Hom}_{\mathcal C}(X,Y)\otimes_kL\bigr)e.
```

The scalar-extension functor sends $X$ to
$(X,\operatorname{id}_X)$. Some sources write
$\mathcal C\otimes_kL$ for this completed category; others use that notation
first for the Hom-space extension and mention the Karoubi envelope separately
[morrison2012noncyclotomic; §2.1](@cite).

If $\mathcal C$ is a Hom-finite Krull–Schmidt category, then $\mathcal C_L$ is
again Hom-finite and idempotent complete, hence Krull–Schmidt
[krause2015krull; Corollary 4.4](@cite). If $\mathcal C$ is finite
semisimple and $k$ is perfect, then $\mathcal C_L$ is again finite
semisimple. Indeed, over a perfect field every finite-dimensional division
algebra is separable over $k$, and a separable algebra remains semisimple after
any field extension. More generally, it is enough that the endomorphism
algebras of the simple objects are separable over $k$; this is the
**absolutely semisimple** hypothesis in
[etingof2012descent; §3.1, footnote 1](@citet).

Perfectness cannot simply be omitted. Let $k$ have characteristic $p>0$, let
$a\in k$ not be a $p$th power, and put
$E=k[u]/(u^p-a)$. The category $\operatorname{Vec}_E$, regarded as a
$k$-linear category, is finite semisimple. After extending from $k$ to $E$,
the endomorphism algebra of its simple object becomes

```math
\label{eq:inseparable-scalar-extension-counterexample}
E\otimes_kE
\cong
E[\varepsilon]/(\varepsilon^p),
```

which is not semisimple. Consequently the scalar-extended category is not
semisimple.

### Scalar-extension interface

| Operation | Meaning |
|:---|:---|
| `extension_of_scalars(C,L; embedding=iota)` | construct a scalar-extended category along $\iota\colon k\to L$ |
| `extension_of_scalars(X,L,D; embedding=iota)` | map an object into a chosen scalar-extended category $\mathcal D$ |
| `extension_of_scalars(f,L,D; embedding=iota)` | map a morphism into $\mathcal D$ |

A concrete category model must implement scalar extension for its categories,
objects, and morphisms. The category-level method constructs a suitable target
model rather than a universal category of formal pairs $(X,e)$. The
object-level method returns the scalar extension of $X$ in that target; it may
be decomposable. For `VectorSpaces` and `GroupRepresentationCategory`, the
target is the full corresponding category over $L$, such as
$\operatorname{Rep}_L(G)$.

Independently constructed isomorphic fields do not necessarily identify their
chosen generators. Pass `embedding=iota` whenever the required embedding is not
the canonical one.

## Splitting

Let $\mathcal C$ now be Hom-finite and Krull–Schmidt. If $X$ is
indecomposable, then $E_X=\operatorname{End}_{\mathcal C}(X)$ is local. Its
residue algebra

```math
\label{eq:residue-endomorphism-division-algebra}
D_X=E_X/\operatorname{rad}(E_X)
```

is a division algebra. After extending scalars,
$E_X\otimes_kL$ can acquire nontrivial idempotents, and
$X_L$ can therefore decompose. Over a separable extension, this decomposition
is governed by the primitive idempotents of $D_X\otimes_kL$, which lift to
$E_X\otimes_kL$.

An indecomposable object is **absolutely indecomposable** if it remains
indecomposable after extension to an algebraic closure. To **split** an object
or a finite family means to choose an extension over which all resulting
indecomposable summands are absolutely indecomposable.

For a finite semisimple category, the indecomposable objects are precisely the
simple objects. A simple object is **absolutely simple** if it remains simple
after every field extension, equivalently after extension to an algebraic
closure. By Schur's lemma,

```math
\label{eq:simple-endomorphism-division-algebra}
D_S=\operatorname{End}_{\mathcal C}(S)
```

is a finite-dimensional division algebra for every simple $S$. The category
is **split** over $k$ if the canonical map

```math
\label{eq:split-simple-condition}
k\longrightarrow D_S
```

is an isomorphism for every simple $S$. This is equivalent to every simple
object being absolutely simple: after scalar extension,
$\operatorname{End}(S_L)=D_S\otimes_kL$. No perfectness hypothesis is needed
for this equivalence. A field $L$ is a **splitting field** for $\mathcal C$ if
$\mathcal C_L$ is split; compare [EGNO; §4.16](@cite).

If $\mathcal C$ is split semisimple, then its Hom-space extension is already
idempotent complete, so the Karoubi envelope introduces no new objects. When
scalar extension remains semisimple, extending a simple object $S$ amounts to
finding primitive idempotents in $D_S\otimes_kL$ and taking their images.
Every simple object of $\mathcal C_L$ arises in this way from a simple object
of $\mathcal C$.

### Splitting interface

| Operation | Meaning |
|:---|:---|
| `split(X)` | choose an extension and split the indecomposable summands of $X$ |
| `split(objects)` | split a finite family over one common extension |
| `split(C)` | construct a split scalar extension when the category model supplies such an algorithm |
| `is_split_semisimple(C)` | test the implemented split-semisimplicity condition |

Computing a splitting field is a category-dependent problem. It can require
endomorphism algebras, their radicals and primitive idempotents, images of
idempotents, and decompositions of the resulting objects. A user implementing
a new category must provide the necessary operations, or a specialized `split`
method, when such algorithms are available.

Over finite fields, the generic methods `split(X)` and `split(objects)` compute
the required extension degree from the residue endomorphism algebras, choose
one extension for the specified object or family, and return the extended
objects together with their decompositions. They work for
$\operatorname{Rep}_k(G)$ and for any other model which provides the required
endomorphism, decomposition, and scalar-extension operations. The resulting
summands are absolutely indecomposable; in a nonsemisimple category they need
not be simple. A category-level `split(C)` remains model-dependent.

The predicate `is_split_semisimple(C)` first tests semisimplicity and then
checks that every enumerated simple object has one-dimensional endomorphism
algebra over the base field. It therefore depends on the implementation's
ability to enumerate simple objects.

The two-dimensional simple representation of $C_3$ over $\mathbb F_2$
provides a small example. Its endomorphism field is $\mathbb F_4$. After
extension to $\mathbb F_4$,

```math
\mathbb F_4\otimes_{\mathbb F_2}\mathbb F_4
\cong
\mathbb F_4\times\mathbb F_4.
```

The two primitive idempotents give the two nontrivial one-dimensional
characters.

```@example finite_field_extension
using TensorCategories, Oscar
G = cyclic_group(3)
C = representation_category(GF(2), G)
V = only(filter(X -> int_dim(X) == 2, simples(C)))

result = split(V)
summands = only(result.decompositions)

@assert result.extension_degree == 2
@assert order(result.field) == 4
@assert length(summands) == 2
@assert all(m == 1 && int_dim(X) == 1 for (X, m) in summands)
summands
```

Extending the coefficient field cannot restore semisimplicity. For example,
no extension of $\mathbb F_5$ makes
$\operatorname{Rep}_{\mathbb F_5}(C_5)$ semisimple.

## Field embeddings, Galois conjugacy, and forms

Let $\iota_1,\iota_2\colon k\hookrightarrow L$. If
$\iota_2=\tau\circ\iota_1$ for an automorphism $\tau$ of $L$, then

```math
\label{eq:semilinear-conjugate-scalar-extensions}
f\otimes\lambda\longmapsto f\otimes\tau(\lambda)
```

defines a $\tau$-semilinear equivalence between the two scalar extensions.
It is an equivalence of ordinary categories, but it is not generally
$L$-linear, and an $L$-linear equivalence need not exist. This is the usual
Galois twist [etingof2012descent; §2.2 and §3](@cite).

Under this hypothesis, if both scalar extensions are split finite semisimple
categories, their underlying $L$-linear categories are nevertheless
noncanonically equivalent:
each is equivalent to $\operatorname{Vec}_L^{\oplus r}$, where $r$ is the
number of simple isomorphism classes.

Conversely, let $\mathcal D$ be an $L$-linear category. A **$k$-form** of
$\mathcal D$ is a $k$-linear category $\mathcal C$ together with an
$L$-linear equivalence

```math
\label{eq:form-of-linear-category}
\mathcal C_L\simeq\mathcal D.
```

Constructing a form is also called **descent**. Thus forms provide the reverse
viewpoint on scalar extension [etingof2012descent; §2.1 and §3](@cite).

Choosing an embedding and extending scalars are different operations. An
embedding selects how the old coefficients map into the target field, whereas
scalar extension can create new idempotents and new direct-sum decompositions.

Continue with [functors and natural transformations](@ref linear-functors).
