# TensorCategories.jl

[![][docs-dev-img]][docs-dev-url][![][ga-img]][ga-url] [![][codecov_img]][codecov_url] [![DOI](https://zenodo.org/badge/413334326.svg)](https://doi.org/10.5281/zenodo.15756399)

TensorCategories is a package under development with the intention to provide a framework as well a examples for computations in the realm of categories.

## Installation

You need to have Julia installed. To use TensorCategories
do the following:

```julia-repl
julia> import Pkg
julia> Pkg.add("TensorCategories")
```

## Usage

TensorCategories relies on the algebraic structures from the [OSCAR-System](https://github.com/oscar-system/Oscar.jl). Here a minimal usage Example.

```@repl
using TensorCategories, Oscar;
C = graded_vector_spaces(QQ, symmetric_goup(3))
Z = center(C)
simples(Z)
smatrix(Z)
```

## Features

TensorCategories provides a vast framework for constructions with finite tensor categories and especially fusion categories. 

# The Center of a fusion category

The current pinnacle feature is the computation of the center of a fusion category in explicit form. The theoretical ground work for 
this approach is layed in [https://arxiv.org/abs/2406.13438](https://arxiv.org/abs/2406.13438)  

## Acknowledgements

This project was started under supervision of [Prof. Ulrich Thiel](https://ulthiel.com/math/)  (University of Kaiserslautern). This work is a
contribution to the SFB-TRR 195 'Symbolic Tools in Mathematics and their
Application' of the German Research Foundation (DFG).


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://fabianmaeurer.github.io/TensorCategories.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://fabianmaeurer.github.io/TensorCategories.jl/dev/

[build-status-img]: https://ci.appveyor.com/api/projects/status/egtv4niuustg4kpc?svg=true
[build-status-url]: https://ci.appveyor.com/project/FabianMaeurer/TensorCategories-jl

[codecov_img]: https://codecov.io/gh/FabianMaeurer/TensorCategories.jl/branch/master/graph/badge.svg?token=axGHAcozx5
[codecov_url]: https://codecov.io/gh/FabianMaeurer/TensorCategories.jl

[ga-img]: https://github.com/FabianMaeurer/TensorCategories.jl/actions/workflows/runtests.yml/badge.svg
[ga-url]: https://github.com/FabianMaeurer/TensorCategories.jl/actions/workflows/runtests.yml
