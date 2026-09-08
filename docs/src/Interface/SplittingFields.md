# [Splitting and scalar extension](@id splitting-and-scalars)

Scalar extension is available for linear categories independently of
semisimplicity. Its effect on direct-sum decompositions is already visible in
Hom-finite Krull–Schmidt categories. For semisimple categories, the same
construction leads to the notions of split categories and splitting fields
used later for fusion categories.

## Scalar extension

Let $\iota\colon k\hookrightarrow L$ be a field embedding. There are two
closely related constructions which are both denoted by scalar extension in
the literature. The first one changes the morphism spaces but keeps the old
objects. We call it the **Hom-space extension** and denote it temporarily by
$\mathcal C\otimes_k^{\mathrm{Hom}}L$. It is defined by

```math
\label{eq:hom-space-scalar-extension}
\operatorname{Ob}(\mathcal C\otimes_k^{\mathrm{Hom}}L)
=\operatorname{Ob}(\mathcal C),
\qquad
\operatorname{Hom}_{\mathcal C\otimes_k^{\mathrm{Hom}}L}(X,Y)
=\operatorname{Hom}_{\mathcal C}(X,Y)\otimes_{k,\iota}L.
```

Composition is extended $L$-bilinearly, and there is a canonical $k$-linear
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

The Hom-space extension is additive and $L$-linear, but it need not be
idempotent complete. Indeed, an algebra
$\operatorname{End}_{\mathcal C}(X)\otimes_kL$ may contain idempotents which
are not split by objects already present in
$\mathcal C\otimes_k^{\mathrm{Hom}}L$. Since every idempotent splits in an
abelian category, the Hom-space extension then cannot be abelian.

This phenomenon is not restricted to semisimple categories. Suppose that
$\mathcal C$ is a Hom-finite Krull–Schmidt $k$-linear category. An
indecomposable object $X$ can become decomposable after extending scalars,
because the local algebra $\operatorname{End}_{\mathcal C}(X)$ may acquire
nontrivial idempotents after tensoring with $L$. The natural scalar extension
among idempotent-complete additive categories therefore includes a second
step:

```math
\label{eq:completed-scalar-extension}
\mathcal C_L
=\operatorname{Kar}\!\left(\mathcal C\otimes_k^{\mathrm{Hom}}L\right).
```

An object of $\mathcal C_L$ is a pair $(X,e)$ with
$e^2=e\in\operatorname{End}_{\mathcal C}(X)\otimes_kL$, and

```math
\label{eq:completed-scalar-extension-hom}
\operatorname{Hom}_{\mathcal C_L}\bigl((X,e),(Y,d)\bigr)
=d\bigl(\operatorname{Hom}_{\mathcal C}(X,Y)\otimes_kL\bigr)e.
```

The scalar-extension functor sends $X$ to $(X,\operatorname{id}_X)$. The
category $\mathcal C_L$ is again Hom-finite and idempotent complete, hence
Krull–Schmidt [krause2015krull; Corollary 4.4](@cite). Decomposing the new
idempotents in $\operatorname{End}_{\mathcal C}(X)\otimes_kL$ gives the
indecomposable summands of the extended object.

For an indecomposable $X$, the algebra
$E_X=\operatorname{End}_{\mathcal C}(X)$ is local and

```math
\label{eq:residue-endomorphism-division-algebra}
D_X=E_X/\operatorname{rad}(E_X)
```

is a division algebra. Over the separable extensions relevant below, the
splitting of $X$ is governed by the primitive idempotents of
$D_X\otimes_kL$, which lift to $E_X\otimes_kL$. Thus the Krull–Schmidt
analogue of a split simple object is an **absolutely indecomposable** object.
The finite-field methods `split(X)` and `split(objects)` use precisely these
residue endomorphism algebras; their resulting summands need not be simple in
a nonsemisimple category.

Some sources write $\mathcal C\otimes_kL$ for the completed category, while
others first use this notation for the Hom-space extension and mention its
idempotent completion separately; see also
[morrison2012noncyclotomic; §2.1](@citet). In this
manual, **scalar extension** means the completed construction in equation
\eqref{eq:completed-scalar-extension}, unless the Hom-space extension is named
explicitly.

## The interface

| Operation | Meaning |
|:---|:---|
| `extension_of_scalars(C,L; embedding=iota)` | construct the scalar-extended category along $\iota\colon k\to L$ |
| `extension_of_scalars(X,L,D; embedding=iota)` | map an object into a chosen scalar-extended category $D$ |
| `split(X)`, `split(objects)` | choose a common extension and split the indecomposable summands of an object or finite family |
| `split(C)` | construct a split form of a category when the category model supplies this operation |

These operations require algorithms supplied by the concrete category model.
In particular, scalar extension must be implemented for its categories,
objects, and morphisms. Finding a splitting field may additionally require
endomorphism algebras, their radicals and primitive idempotents, images of
idempotents, and decomposition of the resulting objects. There is no general
algorithm over an arbitrary coefficient field.

The public function `extension_of_scalars` realizes scalar extension through a
concrete target model rather than by constructing a universal category of
formal pairs $(X,e)$. For `VectorSpaces` and `GroupRepresentationCategory`, it
returns the full category over $L$, such as $\operatorname{Rep}_L(G)$. For a
nonsemisimple representation category this is the abelian scalar extension,
which can in general contain more objects than the additive completion in
equation \eqref{eq:completed-scalar-extension}. For a finite separable
extension $L/k$, the two agree: every module over the extended algebra is a
direct summand of the scalar extension of its restriction to $k$. A
`SixJCategory` is split, so its
coefficient arrays can be transported directly and no new summands occur. For
supported center and relative-center categories, the implementation extends
the known simple objects and explicitly computes the new simple summands. In
the semisimple setting these different Julia representations model the same
completed scalar extension.

The object-level call returns the image of
$X$ under the scalar-extension functor; this image may be decomposable, and
`decompose` finds its summands. Independently constructed isomorphic fields do
not necessarily identify their chosen generators, so a noncanonical embedding
should be passed explicitly.

Over finite fields, the generic methods `split(X)` and
`split(objects)` compute the required degree from
$\operatorname{End}(X)/\operatorname{rad}\operatorname{End}(X)$, choose one
extension for a specified object or finite family, and return the extended
objects together with their decompositions. They work for
$\operatorname{Rep}_k(G)$ and for any other model providing the operations
listed above. The category-level method `split(C)` chooses a common splitting
field for a supported center category. Separate
`karoubian_envelope` methods are available for center and relative-center
models; scalar extension of a center already performs the required splitting
of its known simple objects.

!!! note "Beyond the semisimple setting"
    Equation \eqref{eq:completed-scalar-extension} is the appropriate additive
    scalar extension of a Hom-finite Krull–Schmidt category and records all
    new direct-sum decompositions. If $\mathcal C$ is nonsemisimple abelian,
    however, this Karoubi envelope need not be abelian: it adds images of
    idempotents, but not arbitrary missing kernels, cokernels, or extensions.
    The abelian scalar extension is instead the Deligne tensor product
    $\mathcal C\boxtimes_k\operatorname{Vec}_L$ when it exists. For
    $\mathcal C\simeq A\text{-mod}$ it is
    $(A\otimes_kL)\text{-mod}$; compare
    [lopezfranco2013tensor; Theorem 3 and Example 11](@citet).
    When $L/k$ is finite separable, every
    $(A\otimes_kL)$-module is a direct summand of the scalar extension of its
    restriction to $k$, so the additive and abelian constructions agree in
    this case.

## Split semisimple categories

Now let $\mathcal C$ be a finite semisimple $k$-linear category. By Schur's
lemma, the endomorphism algebra

```math
\label{eq:simple-endomorphism-division-algebra}
D_S=\operatorname{End}_{\mathcal C}(S)
```

of a simple object $S$ is a finite-dimensional division algebra over $k$.
The category is **split** over $k$ if the canonical map

```math
\label{eq:split-simple-condition}
k\longrightarrow D_S
```

is an isomorphism for every simple $S$. This is automatic when $k$ is
algebraically closed, but it need not hold over a number field. A simple
object is **absolutely simple** if it remains simple after extension to an
algebraic closure. Over a perfect field this agrees with the
scalar-endomorphism condition in equation
\eqref{eq:split-simple-condition}. Over an imperfect field, inseparability
requires additional care; compare [EGNO; §4.16](@cite).

If $k$ is perfect, or more generally if the endomorphism algebras of the
simple objects are separable over $k$, then $\mathcal C_L$ is again a finite
semisimple $L$-linear category. In this case it agrees with the Deligne
product $\mathcal C\boxtimes_k\operatorname{Vec}_L$; see
[etingof2012descent; §3.1](@cite) and
[lopezfranco2013tensor; Theorem 3 and §5](@cite). A weak fusion category
therefore remains weak fusion under this hypothesis.

The field $L$ is a **splitting field** for $\mathcal C$ if $\mathcal C_L$
is split. Extending a simple $S$ to $L$ amounts to finding primitive
idempotents in $D_S\otimes_kL$ and taking their images. Every simple object
of $\mathcal C_L$ arises in this way, and scalar extensions of two
nonisomorphic simple objects have no common simple summand.

The predicate `is_split_semisimple(C)` checks that the category is
semisimple and that every enumerated simple object has one-dimensional
endomorphism algebra over the base field. It depends on the backend's
ability to enumerate simple objects and does not separately certify absolute
simplicity over an imperfect field.

The two-dimensional simple representation of $C_3$ over $\mathbb F_2$
provides a small example. Its endomorphism field is $\mathbb F_4$. After
extension to $\mathbb F_4$, this endomorphism algebra becomes
$\mathbb F_4\otimes_{\mathbb F_2}\mathbb F_4\cong\mathbb F_4\times\mathbb F_4$.
The two primitive idempotents do not split in the Hom-space extension, whereas
the completed scalar extension contains their images: the two nontrivial
one-dimensional characters.

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

## Field embeddings and Galois conjugacy

As explained in the [coefficient-field section](@ref base-fields), an abstract
number field has no preferred embedding into $\mathbb C$. The embedding
$\iota\colon k\hookrightarrow L$ is therefore part of a scalar extension, not
an implicit identification of the two fields. Applying a field automorphism to
all coefficients of an algebraic category model gives its **Galois
conjugate**; the defining algebraic equations are preserved.

More precisely, let $\iota_1,\iota_2\colon k\hookrightarrow L$. If
$\iota_2=\tau\circ\iota_1$ for an automorphism $\tau$ of $L$, then

```math
\label{eq:semilinear-conjugate-scalar-extensions}
f\otimes\lambda\longmapsto f\otimes\tau(\lambda)
```

defines a $\tau$-semilinear equivalence between the two scalar extensions.
It is an equivalence of ordinary categories, but it is not generally
$L$-linear, and an $L$-linear equivalence need not exist. This is the usual
Galois twist; compare [etingof2012descent; §2.2 and §3](@citet).

When such a $\tau$ exists and both extensions are split finite semisimple
categories, their underlying $L$-linear categories are nevertheless
noncanonically equivalent: each is equivalent to
$\operatorname{Vec}_L^{\oplus r}$, where $r$ is the number of simple
isomorphism classes.

Conversely, let $\mathcal D$ be an $L$-linear category. A **$k$-form** of
$\mathcal D$ is a $k$-linear category $\mathcal C$ together with an
$L$-linear equivalence

```math
\label{eq:form-of-linear-category}
\mathcal C_L\simeq\mathcal D.
```

Constructing such a form is also called **descent**. Thus forms provide the
reverse viewpoint on scalar extension; see
[etingof2012descent; §2.1 and §3](@citet).

Changing the embedding and extending scalars are nevertheless different
operations. An embedding selects a realization of the existing coefficients,
whereas scalar extension can create new idempotents and new direct-sum
decompositions.

Positive characteristic introduces a second, independent issue. Extending the
field can make simple objects split, but it does not turn a nonsemisimple
category into a semisimple one. For instance, extending scalars in
$\operatorname{Rep}_{\mathbb F_5}(C_5)$ does not restore Maschke's theorem.

Continue with [functors and natural transformations](@ref linear-functors).
