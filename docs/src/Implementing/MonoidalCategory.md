# [Implementing a monoidal category](@id implementing-monoidal)

To add a monoidal structure to a category model, implement the tensor product
on both objects and morphisms, a tensor unit, and the associator described in
the preceding chapter:

| Required method | Meaning |
|:---|:---|
| `tensor_product(X,Y)` | the object $X\otimes Y$ |
| `tensor_product(f,g)` | the morphism $f\otimes g$ |
| `one(C)` | the tensor unit $\mathbb 1$ |
| `associator(X,Y,Z)` | $a_{X,Y,Z}\colon (X\otimes Y)\otimes Z\to X\otimes(Y\otimes Z)$ |

The object and morphism methods for `tensor_product` must define the same
bifunctor. In particular, they must preserve identities and composition. The
represented tensor unit is strict: tensoring a represented object with
$\mathbb 1$ must return that object. The associator must be natural and satisfy
the pentagon equation. These requirements are mathematical obligations of the
implementation; defining methods with the correct signatures does not verify
them.

For a category with a matrix realization, the implementation must also specify
the ordered basis of $U(X\otimes Y)$. In the built-in model of vector spaces,
the coordinate from the right tensor factor varies fastest. The matrix of a
tensor product is therefore the Kronecker product

```math
\label{eq:implemented-tensor-kronecker}
M_{f\otimes g}=M_f\mathbin{\operatorname{\otimes}_{\mathrm{Kr}}}M_g.
```

The induced bases on $(X\otimes Y)\otimes Z$ and
$X\otimes(Y\otimes Z)$ have the same order, so the associator is represented
by an identity matrix:

```@example monoidal_implementation
using TensorCategories, Oscar
C = vector_spaces(QQ)
X = VectorSpaceObject(C, 2)
Y = VectorSpaceObject(C, 3)
f = morphism(X, Y, matrix(QQ, [1 0 0; 0 0 0]))

@assert int_dim(X ⊗ Y) == 6
@assert matrix(f ⊗ id(X)) ==
    kronecker_product(matrix(f), matrix(id(X)))
@assert associator(X,Y,X) == id((X⊗Y)⊗X)
is_monoidal(C)
```

This identity associator is a feature of the chosen coordinate model. A
different ordering of tensor-product bases requires the corresponding change
of basis, and a general monoidal category may have a nontrivial associator or
no matrix realization at all.

An implementation should declare `is_monoidal(C) == true` only after supplying
this structure and verifying bifunctoriality, naturality, the strict unit
conditions, and the pentagon. The declaration does not perform those checks
automatically.

Continue with [rigidity](@ref rigid-categories).
