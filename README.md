# TensorCategories.jl

[![][docs-stable-img]][docs-stable-url][![][ga-img]][ga-url] [![][codecov_img]][codecov_url]

TensorCategories is a package under development with the intention to provide a framework as well a examples for computations in the realm of categories.

## Installation

You need to have Julia installed. For reliable results Julia version at least 1.6 is required. To use TensorCategories
do the following:

```julia-repl
julia> import Pkg
julia> Pkg.add(url = "https://github.com/FabianMaeurer/TensorCategories.jl")
```

## Usage

To use TensorCategories the structures from the [OSCAR-System](https://github.com/oscar-system/Oscar.jl) are required. Here a minimal usage Example.

```@repl
using TensorCategories, Oscar;
F,_ = FiniteField(5)
G = symmetric_group(2)
X = gset(G,[1,2,3])
C = ConvolutionCategory(X,F)
simples(C)
grothendieck_ring(C)
```

## Acknowledgements

This project was started under supervision of [Prof. Ulrich Thiel](https://ulthiel.com/math/)  (University of Kaiserslautern). This work is a
contribution to the SFB-TRR 195 'Symbolic Tools in Mathematics and their
Application' of the German Research Foundation (DFG).


[docs-stable-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-url]: https://fabianmaeurer.github.io/TensorCategories.jl/

[build-status-img]: https://ci.appveyor.com/api/projects/status/egtv4niuustg4kpc?svg=true
[build-status-url]: https://ci.appveyor.com/project/FabianMaeurer/TensorCategories-jl

[codecov_img]: https://codecov.io/gh/FabianMaeurer/TensorCategories.jl/branch/master/graph/badge.svg?token=axGHAcozx5
[codecov_url]: https://codecov.io/gh/FabianMaeurer/TensorCategories.jl

[ga-img]: https://github.com/FabianMaeurer/TensorCategories.jl/actions/workflows/runtests.yml/badge.svg
[ga-url]: https://github.com/FabianMaeurer/TensorCategories.jl/actions/workflows/runtests.yml
