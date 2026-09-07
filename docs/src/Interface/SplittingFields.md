# [Splitting and scalar extension](@id splitting-and-scalars)

Let $\mathcal C$ be a finite semisimple $k$-linear category. By Schur's
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
algebraically closed, but it need not hold over a number field. A simple object
is **absolutely simple** if it remains simple after extension to an algebraic
closure. Over the perfect fields used in the examples below, this agrees with
the scalar-endomorphism condition in equation
\eqref{eq:split-simple-condition}. Over an imperfect field, inseparability
requires additional care; compare [EGNO; §4.16](@cite).

The predicate `is_split_semisimple(C)` tests the implemented version of this
condition. It first requires semisimplicity and then checks that every
enumerated simple object has one-dimensional endomorphism algebra over the base
field. It does not separately certify absolute simplicity over an imperfect
field, and it depends on the backend's ability to enumerate simples.

## Scalar extension

Given an embedding $\iota\colon k\hookrightarrow K$, scalar extension forms
a $K$-linear category $\mathcal C_K$ and a $k$-linear functor

```math
\label{eq:categorical-scalar-extension}
K\otimes_k-\colon\mathcal C\longrightarrow\mathcal C_K.
```

For a concrete matrix model, TensorCategories.jl applies $\iota$ to all
matrix entries and reconstructs the object or morphism in its new parent. The
corresponding calls are `extension_of_scalars(C,K; embedding=iota)` and
`extension_of_scalars(X,K,D; embedding=iota)`, with methods depending on the
category model. Independently constructed isomorphic fields do not necessarily
identify their chosen generators, so a noncanonical embedding should be passed
explicitly.

Scalar extension can split a simple object. Algebraically, one extends
$D_S$ to $K$, finds primitive idempotents in $K\otimes_kD_S$, and takes
their images. The same procedure applied to endomorphism algebras decomposes
arbitrary extended objects. The package uses this principle in its splitting
algorithms. The one-argument call `split(C)` is currently implemented for
center categories; generic object splitting over finite fields is available
through `split(X)` and `split(objects)` for the supported object types.

The two-dimensional simple representation of $C_3$ over $\mathbb F_2$
provides a small example. Its endomorphism field is $\mathbb F_4$, and after
extension to $\mathbb F_4$ it decomposes as the sum of the two nontrivial
one-dimensional characters:

```@example finite_field_extension
using TensorCategories, Oscar
G = cyclic_group(3)
C = representation_category(GF(2), G)
V = only(filter(X -> int_dim(X) == 2, simples(C)))

K = GF(2, 2)
D = extension_of_scalars(C, K)
V_K = extension_of_scalars(V, K, D)
summands = decompose(V_K)

@assert length(summands) == 2
@assert all(m == 1 && int_dim(X) == 1 for (X, m) in summands)
summands
```

## Complex embeddings and Galois conjugacy

An abstract number field does not come with a preferred embedding into
$\mathbb C$. If $K=\mathbb Q(s)$ with $s^2=2$, its two complex
embeddings send $s$ to the two square roots of $2$. OSCAR represents their
values by certified complex balls:

```@example complex_embeddings
using Oscar
K, s = quadratic_field(2)
embeddings = complex_embeddings(K)
@assert length(embeddings) == 2
@assert overlaps(embeddings[1](s), -embeddings[2](s))
length(embeddings)
```

Applying an embedding to every coefficient of a system of polynomial
structural equations preserves those equations. The resulting realization is
a **Galois conjugate**. Its decomposition rules can remain unchanged while
signs, phases, dimensions, and positivity or unitarity properties change.
Different embeddings are restrictions of automorphisms after passing to a
normal closure; the original coefficient field need not itself be Galois.

Choosing a complex embedding and extending the coefficient field are distinct
operations. An embedding selects a realization of existing exact coefficients.
A scalar extension can additionally create new idempotents and hence new
direct-sum decompositions.

Positive characteristic introduces a second, independent issue. Extending the
field can make simple objects split, but it does not turn a nonsemisimple
category into a semisimple one. For instance, extending scalars in
$\operatorname{Rep}_{\mathbb F_5}(C_5)$ does not restore Maschke's theorem.

Continue with [functors and natural transformations](@ref linear-functors).
