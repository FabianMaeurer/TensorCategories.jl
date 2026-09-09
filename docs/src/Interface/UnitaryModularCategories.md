# [Unitarity and modularity](@id unitary-modular-categories)

## [Dagger structures and unitarity](@id unitary-categories)

For a category realized over $\mathbb C$, a dagger sends a morphism
$f\colon X\to Y$ to $f^\dagger\colon Y\to X$, is conjugate-linear and
involutive, reverses composition, and is compatible with tensor products. A
unitary category has a compatible positive dagger structure.

### The interface

The method `dagger(f)` is category-specific. In concrete matrix models it is
implemented by conjugate transpose when the chosen bases carry the intended
Hermitian structures. The predicate `is_unitary(C)` concerns the structure in
those chosen coordinates; it does not search for a change of basis. Exact
coefficients in an abstract number field do not by themselves choose complex
conjugation or positivity, so a complex embedding or complex ball realization
is needed for analytic unitarity questions.

## [Premodular and modular categories](@id premodular-categories)

A premodular category is a braided spherical fusion category. For ordered
simple representatives $X_1,\ldots,X_r$, TensorCategories.jl uses the
unnormalized matrix

```math
\label{eq:categorical-s-matrix}
S_{ij}=\operatorname{Tr}\!\left(
c_{X_i,X_j}\circ c_{X_j,X_i}
\right).
```

The displayed endomorphism acts on $X_j\otimes X_i$. Cyclicity of the
spherical trace identifies this with the usual trace of the double braiding on
$X_i\otimes X_j$. A premodular category is **modular** when this matrix is
invertible [EGNO; Definitions 8.13.1, 8.13.2, and 8.13.4](@cite).

### The interface

The public methods use the following conventions:

- `smatrix(C)` returns the unnormalized matrix in equation
  \eqref{eq:categorical-s-matrix}, ordered by `simples(C)`;
- `normalized_smatrix(C)` returns
  $S/\sqrt{\dim(\mathcal C)}$, including the required choice of square root;
- `tmatrix(C)` is diagonal with entries `twist_scalar(X)` in the same order;
- the generic `is_modular(C)` requires `is_fusion(C)`, `is_braided(C)`, and
  `is_spherical(C)` and then tests whether `smatrix(C)` is invertible.

Thus the generic predicate presently imposes the split fusion condition. A
non-split weak fusion category can satisfy an analogous nondegeneracy condition,
but `is_modular(C)` does not report it as modular unless the category-specific
implementation supplies a different method.

For `ArbField` and `AcbField`, invertibility is decided by checking that the
determinant ball excludes zero. A successful result is rigorous for the input
enclosures at the selected precision; failure to exclude zero is inconclusive
and can require a higher-precision realization.

Continue with [monoidal functors](@ref monoidal-functors).
