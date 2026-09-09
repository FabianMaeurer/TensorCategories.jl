# [Tensor products and associators](@id monoidal-categories)

A monoidal category has a tensor product bifunctor, a tensor unit, and coherent
associativity and unit isomorphisms. We use the conventions of
[EGNO; §2.2](@cite), while recording the package's stricter presentation of the
unit.

## Tensor products and the unit

The tensor product acts on both objects and morphisms. If $f\colon X\to X'$
and $g\colon Y\to Y'$, then

```math
\label{eq:tensor-product-morphism}
f\otimes g\colon X\otimes Y\longrightarrow X'\otimes Y'.
```

The implementation must satisfy the interchange law

```math
\label{eq:tensor-bifunctoriality}
(f'\circ f)\otimes(g'\circ g)
=(f'\otimes g')\circ(f\otimes g).
```

A general monoidal category has left and right unit constraints.
TensorCategories.jl uses a unit-strict presentation: tensoring a represented
object with $\mathbb 1$ returns that object, and the unit constraints are
identities. The associator need not be an identity.

## The associator

The package fixes the direction

```math
\label{eq:monoidal-associator}
a_{X,Y,Z}\colon (X\otimes Y)\otimes Z
\longrightarrow X\otimes(Y\otimes Z).
```

This direction agrees with [EGNO; Definition 2.2.8](@citet). Parentheses should
be written explicitly even when the two bracketings happen to be equal as
represented objects.

The associator satisfies Mac Lane's pentagon equation. With the direction in
equation \eqref{eq:monoidal-associator}, it is

```math
\label{eq:monoidal-pentagon}
(\operatorname{id}_X\otimes a_{Y,Z,W})
\circ a_{X,Y\otimes Z,W}
\circ(a_{X,Y,Z}\otimes\operatorname{id}_W)
=
a_{X,Y,Z\otimes W}\circ a_{X\otimes Y,Z,W}.
```

Both sides map $((X\otimes Y)\otimes Z)\otimes W$ to
$X\otimes(Y\otimes(Z\otimes W))$. Writing the sources and targets is a
useful check when translating associator formulas from another convention.

## The interface

| Operation | Meaning |
|:---|:---|
| `tensor_product(X,Y)`, `X ⊗ Y` | the tensor product $X\otimes Y$ |
| `tensor_product(f,g)`, `f ⊗ g` | the tensor product of morphisms |
| `one(C)` | the tensor unit $\mathbb 1$ |
| `associator(X,Y,Z)` | the associator in equation \eqref{eq:monoidal-associator} |
| `inv_associator(X,Y,Z)` | the inverse associator |
| `pentagon_axiom(C)` | exhaustively check the pentagon on the listed simple objects when supported |
| `randomized_pentagon_axiom(C,n)` | check the pentagon on $n$ sampled quadruples |

An implementation must provide the tensor product on objects and morphisms,
the tensor unit, and the associator, and it is responsible for the
bifunctoriality and coherence axioms. Merely providing an `associator` method
does not establish the pentagon. The randomized check is a diagnostic rather
than an exhaustive verification.

## Example: Vector spaces and representations

For the implemented categories $\operatorname{Vec}_k$ and
$\operatorname{Rep}_k(G)$, tensor-product bases are ordered so that the
coordinate from the right tensor factor varies fastest. Consequently,

```math
\label{eq:concrete-tensor-kronecker}
M_{f\otimes g}=M_f\mathbin{\operatorname{\otimes}_{\mathrm{Kr}}}M_g,
```

where the right-hand side is the Kronecker product in OSCAR. The canonical
rebracketing of these bases is represented by an identity matrix.

```@example concrete_monoidal
using TensorCategories, Oscar
C = vector_spaces(QQ)
X = VectorSpaceObject(C, 2)
Y = VectorSpaceObject(C, 3)
f = morphism(X, X, matrix(QQ, [1 1; 0 1]))
g = id(Y)
@assert matrix(f ⊗ g) == kronecker_product(matrix(f), matrix(g))
@assert associator(X,Y,X) == id((X⊗Y)⊗X)
int_dim(X ⊗ Y)
```

These formulas belong to these concrete models. A general monoidal category
can have a nontrivial associator and need not provide matrices at all.

Continue with [implementing a monoidal category](@ref implementing-monoidal).
