# TensorCategories.jl


[![Citation](https://img.shields.io/badge/Citation-arXiv%3A2406.13438-B31B1B?logo=arxiv)](https://arxiv.org/abs/2406.13438)

TensorCategories.jl is an open-source software package for computations with tensor categories, especially [fusion categories](@ref tensor-conventions). Built on the [Julia](https://julialang.org/) programming language and the [OSCAR](https://www.oscar-system.org/) computer algebra system, it follows the standard mathematical framework described by [EGNO](@citet): objects, morphisms, tensor products, associators, and other categorical structures are represented as such, while concrete combinatorial descriptions, such as [$F$-symbols](@ref skeletal-fusion), are also supported. The package supports exact symbolic computations over a range of [coefficient fields](@ref base-fields), including number fields and finite fields of positive characteristic, as well as [numerical computations](@ref numerical-computations) intended for applications in mathematical physics such as anyon models and conformal field theory.

Current highlights include:

* A general, extensible [category framework](@ref category-interface) together with additive, linear, abelian, monoidal, tensor, and fusion structures.

* Support for [skeletal fusion categories](@ref skeletal-fusion), including exact and numerical access to $F$-symbols, $R$-symbols, pivotal data, and related invariants.

* Integration of fusion-category data from [AnyonWiki](F-symbols/AnyonWiki.md), providing access to a large collection of fusion categories.

* A generic algorithm for computing [Drinfeld centers](@ref center), producing explicit central objects with half-braidings rather than only abstract equivalence classes [maurer2024computing](@cite).

* Computation of the Drinfeld centers, including $F$-symbols and $R$-symbols, of the 76 monoidal presentations underlying the 279 multiplicity-free AnyonWiki entries of rank at most 5; the results are stored in our [TensorCategoriesDatabase](https://github.com/TensorCategories/TensorCategoriesDatabase). Some of these centers yield new [modular categories](@ref premodular-categories) with high rank and nontrivial fusion multiplicities, for example rank 21 and fusion multiplicity 2 [maurer2024computing](@cite).

* Explicit computation of $F$-symbols, $R$-symbols, and pivotal coefficients for the [Drinfeld center of the Haagerup subfactor](F-symbols/Haagerup.md) [maurer2026haagerup](@cite).


## Showcase

Here is a showcase example computing the [center](@ref center)
$\mathcal{Z}(\mathcal{C})$ of the Ising fusion category $\mathcal{C}$ over the
field $\mathbb{Q}(\sqrt{2})$. The computation shows that
$\mathcal{Z}(\mathcal{C})$ is [not split](@ref splitting-and-scalars) over
$\mathbb{Q}(\sqrt{2})$, i.e. some simple objects will decompose after scalar
extension to $\mathbb{C}$. We then compute the multiplication table of its
[Grothendieck ring](@ref grothendieck-rings) and the
[$S$-matrix](@ref premodular-categories) of this non-split modular category.
Simple-object enumeration uses randomized algebra algorithms, so the order
below is one possible output.

```julia-repl
julia> using TensorCategories, Oscar

julia> K,r2 = quadratic_field(2)
(Real quadratic field defined by x^2 - 2, sqrt(2))

julia> C = ising_category(K,r2)

julia> simples(C)
3-element Vector{SixJObject}:
 𝟙
 χ
 X

julia> Z = center(C)
Drinfeld center of Ising fusion category

julia> S = simples(Z)
5-element Vector{CenterObject}:
 Central object: 𝟙
 Central object: 𝟙
 Central object: 𝟙 ⊕ χ
 Central object: 2⋅χ
 Central object: 4⋅X

julia> T = only([T for T in S if int_dim(End(T)) == 2]);

julia> End(T)
Vector space of dimension 2 over Real quadratic field defined by x^2 - 2.

julia> print_multiplication_table(S, ["X$i" for i in eachindex(S)])
5×5 Matrix{String}:
 "X1"  "X2"  "X3"            "X4"           "X5"
 "X2"  "X1"  "X3"            "X4"           "X5"
 "X3"  "X3"  "X1 ⊕ X2 ⊕ X4"  "2⋅X3"         "2⋅X5"
 "X4"  "X4"  "2⋅X3"          "2⋅X1 ⊕ 2⋅X2"  "2⋅X5"
 "X5"  "X5"  "2⋅X5"          "2⋅X5"         "4⋅X1 ⊕ 4⋅X2 ⊕ 8⋅X3 ⊕ 4⋅X4"

julia> smatrix(Z)
[        1            1    2    2    4*sqrt(2)]
[        1            1    2    2   -4*sqrt(2)]
[        2            2    0   -4            0]
[        2            2   -4    4            0]
[4*sqrt(2)   -4*sqrt(2)    0    0            0]
```

The two-dimensional endomorphism algebra of $T$ shows
why this center is not split over $\mathbb Q(\sqrt2)$: a split simple would
have endomorphism algebra equal to the coefficient field. The five displayed
objects become nine simple objects after extension to a splitting field; the
[center chapter](@ref ising-center) carries out that computation.

## Installation

You need to have [Julia](https://julialang.org/downloads/) installed. To install
TensorCategories.jl, run:

```julia-repl
julia> import Pkg
julia> Pkg.add("TensorCategories")
```

This also installs all dependencies, including [OSCAR](https://www.oscar-system.org/).


## How to cite

If TensorCategories.jl contributes to your research, please cite the paper that introduced the software:

```bibtex
@misc{MaeurerThiel2024ComputingCenter,
  author        = {M{\"a}urer, Fabian and Thiel, Ulrich},
  title         = {Computing the center of a fusion category},
  year          = {2024},
  eprint        = {2406.13438},
  archivePrefix = {arXiv},
  primaryClass  = {math.RT},
  doi           = {10.48550/arXiv.2406.13438}
}
```

The software itself is archived on Zenodo and can be cited as follows:

```bibtex
@software{Maeurer2026TensorCategories,
  author    = {M{\"a}urer, Fabian},
  title     = {{TensorCategories.jl}},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.18760250},
  url       = {https://doi.org/10.5281/zenodo.18760250}
}
```

## License

The TensorCategories.jl package is licensed under the GNU General Public License v3.0 or later.

Copyright (c) 2021 Fabian Mäurer and contributors.

See [`LICENSE`](LICENSE) for the full license text.
See [`COPYRIGHT`](COPYRIGHT) for copyright information.

## Acknowledgements

TensorCategories.jl was initiated by [**Ulrich Thiel**](https://agag-thiel.math.rptu.de/math/) (RPTU University Kaiserslautern-Landau) within his project A20 "Towards unipotent character sheaves associated to Coxeter groups" (2020–2024) of the SFB-TRR 195 ["Symbolic Tools in Mathematics and their Application"](https://www.computeralgebra.de/sfb/), funded by the German Research Foundation (DFG). The package was created and developed by **Fabian Mäurer** as part of his Master's and PhD work under Thiel's supervision (2021–2026); see [maeurer2026thesis](@citet). Its development is currently supported by Thiel's project A20 "Categorical representation theory" (2024–2028) in the SFB-TRR 195. Additional support is provided by the Forschungsinitiative "SymbTools" of the state of Rheinland-Pfalz, in which Thiel is one of the project leaders.

[**Gert Vercleyen**](https://gert-vercleyen.github.io/) contributed to the integration of the data from his [AnyonWiki](https://anyonwiki.github.io/).

Since version 0.7, AI assistance has been used in the development of TensorCategories.jl.
