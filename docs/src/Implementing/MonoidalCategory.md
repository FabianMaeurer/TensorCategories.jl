# [Implementing a monoidal category](@id implementing-monoidal)

We extend the matrix category from the preceding chapter by giving it the
standard tensor product of vector spaces. The complete version is available as
[matrix_category.jl](matrix_category.jl); its linear and abelian methods agree
with the earlier tutorial.

On objects and morphisms the additional methods are:

```julia
function TensorCategories.tensor_product(X::MatObject, Y::MatObject)
    parent(X) == parent(Y) || throw(ArgumentError("different categories"))
    MatObject(parent(X), X.n*Y.n)
end

TensorCategories.tensor_product(f::MatMorphism, g::MatMorphism) =
    morphism(domain(f)⊗domain(g), codomain(f)⊗codomain(g),
             kronecker_product(matrix(f), matrix(g)))

Base.one(C::MatCategory) = MatObject(C, 1)

TensorCategories.associator(X::MatObject, Y::MatObject, Z::MatObject) =
    id((X⊗Y)⊗Z)
```

The ordered tensor basis has the coordinate from the right factor varying
fastest. This makes the morphism formula the usual Kronecker product and makes
the canonical rebracketing matrix an identity. Other basis orders would require
corresponding permutation matrices; the identity associator is a feature of
this representation, not part of the definition of a monoidal category.

The implementation declares `is_monoidal(C) == true` only after providing the
tensor product on both objects and morphisms, the unit, and the associator. The
declaration does not check bifunctoriality or the pentagon automatically.

```@example monoidal_implementation
using TensorCategories, Oscar
include("matrix_category.jl")
using .MatrixCategoryTutorial
C = MatCategory(QQ)
X, Y = MatObject(C, 2), MatObject(C, 3)
f = morphism(X, Y, matrix(QQ, [1 0 0; 0 0 0]))
@assert int_dim(X ⊗ Y) == 6
@assert matrix(f ⊗ id(X)) ==
    kronecker_product(matrix(f), matrix(id(X)))
@assert associator(X,Y,X) == id((X⊗Y)⊗X)
is_monoidal(C)
```

This tutorial has not implemented duality, so it does not declare rigidity or
any of the stronger structures introduced next.

Continue with [rigidity](@ref rigid-categories).
