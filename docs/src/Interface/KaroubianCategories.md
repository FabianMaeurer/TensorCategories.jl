# [Idempotents and Krull–Schmidt categories](@id karoubian-categories)

Let $\mathcal C$ be an additive category. An endomorphism
$e\colon X\to X$ is **idempotent** if $e^2=e$. It **splits** if there are an
object $Y$ and morphisms

```math
Y\xrightarrow{i}X\xrightarrow{p}Y
```

such that

```math
\label{eq:idempotent-splitting}
p\circ i=\operatorname{id}_Y,
\qquad
i\circ p=e.
```

An additive category is **Karoubian**, or **idempotent complete**, if every
idempotent splits. In an abelian category this is automatic: one may take
$Y=\operatorname{im}(e)$. For an arbitrary additive category, its Karoubi
envelope formally adjoins objects $(X,e)$ for idempotents
$e\in\operatorname{End}(X)$; a morphism $(X,e)\to(Y,d)$ is a morphism
$f\colon X\to Y$ satisfying $f=d\circ f\circ e$.

Idempotents encode direct-sum decompositions. If

```math
\label{eq:direct-sum-components}
X\cong X_1\oplus\cdots\oplus X_n
```

has inclusions $i_r\colon X_r\to X$ and projections
$p_r\colon X\to X_r$, then

```math
\label{eq:orthogonal-idempotents}
p_r\circ i_r=\operatorname{id}_{X_r},
\qquad
p_r\circ i_s=0\quad(r\ne s),
\qquad
\operatorname{id}_X=\sum_{r=1}^n i_r\circ p_r.
```

Consequently, $e_r=i_r\circ p_r$ are pairwise orthogonal idempotents whose
sum is $\operatorname{id}_X$. Conversely, splitting any such family recovers
the direct sum in equation \eqref{eq:direct-sum-components}. Decomposing an
object into indecomposable summands therefore amounts to decomposing its
identity into primitive pairwise orthogonal idempotents.

An additive category is a **Krull–Schmidt category** if every object is a
finite direct sum of objects with local endomorphism rings. These summands are
indecomposable, and the decomposition is unique up to permutation and
isomorphism. In this setting an object is indecomposable precisely when its
endomorphism ring is local. An additive category with split idempotents and
semiperfect endomorphism rings is Krull–Schmidt
[krause2015krull; Corollary 4.4](@cite). In particular, a Hom-finite
$k$-linear additive category is Krull–Schmidt precisely when it is idempotent
complete. Every locally finite abelian category is therefore Krull–Schmidt;
compare [EGNO; Definition 1.8.1 and the paragraph following it](@cite).

The number of indecomposable isomorphism classes is a separate finiteness
question. A finite abelian category has only finitely many simple isomorphism
classes, but it may have infinitely many indecomposable ones. For an
algebraically closed field $k$ of characteristic $p$, the category
$\operatorname{Rep}_k(G)$ has finite representation type precisely when the
Sylow $p$-subgroups of $G$ are cyclic [higman1954indecomposable](@cite).
Thus finite cyclic groups have only finitely many indecomposable
representations, whereas the Klein four group in characteristic $2$ already
has infinitely many.

## The interface

As with simple-object computations, determining indecomposability and finding
a decomposition are generally difficult, category-specific problems. An
implementation may provide:

| Operation | Meaning |
|:---|:---|
| `is_indecomposable(X)` | test whether $X$ is indecomposable |
| `decompose(X)` | return pairs `(Y,m)` of indecomposable summands and their multiplicities |
| `TensorCategories.is_krull_schmidt(C)` | record that the implementation treats $\mathcal C$ as Krull–Schmidt |

Over a finite field, the generic decomposition backend forms
$A=\operatorname{End}_{\mathcal C}(X)$ and decomposes the regular right
$A$-module. Projecting $1_A$ onto its indecomposable summands gives primitive
idempotents in $A$, and their images give the indecomposable summands of $X$.
The group-representation backend instead applies GAP's
`MTX.Indecomposition` directly to the representation. Hecke.jl's `ModAlgAss`
MeatAxe routines provide related irreducibility and composition-series
algorithms, as discussed in the preceding section.

## Example: Modular group representations

Let $G=C_5$ and $k=\mathbb F_5$. This category has one simple isomorphism
class, while the unipotent Jordan blocks $J_r$, for $1\leq r\leq5$, are its
five indecomposable isomorphism classes. The two-dimensional block $J_2$ is
not simple, and
$\operatorname{End}(J_2)\cong k[t]/(t^2)$ is local.

```@example krull_schmidt_representations
using TensorCategories, Oscar

function jordan_representation(C, n)
    J = identity_matrix(base_ring(C), n)
    for i in 1:n-1
        J[i, i+1] = 1
    end
    Representation(C, gens(base_group(C)), [J])
end

C = representation_category(GF(5), cyclic_group(5))
J1 = jordan_representation(C, 1)
J2 = jordan_representation(C, 2)
@assert is_indecomposable(J2) && !is_simple(J2)

X, inclusions, projections = direct_sum(J1, J2)
e1 = inclusions[1] ∘ projections[1]
e2 = inclusions[2] ∘ projections[2]
@assert e1 ∘ e1 == e1 && e2 ∘ e2 == e2
@assert e1 ∘ e2 == zero_morphism(X, X)
@assert e2 ∘ e1 == zero_morphism(X, X)
@assert e1 + e2 == id(X)

sort([(int_dim(Y), m) for (Y, m) in decompose(X)])
```

The two primitive idempotents split $X$ into summands of dimensions $1$ and
$2$, and `decompose(X)` returns `[(1,1),(2,1)]` after recording each summand by
its dimension and multiplicity. This is a direct-sum decomposition into
indecomposables; it is different from the composition series of $J_2$
discussed in the preceding section.

Continue with [finite and semisimple categories](@ref semisimple-categories).
