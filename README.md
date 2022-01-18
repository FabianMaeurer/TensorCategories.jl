# JuCat.jl

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url]|                                  | [![Build status](https://ci.appveyor.com/api/projects/status/egtv4niuustg4kpc?svg=true)](https://ci.appveyor.com/project/FabianMaeurer/jucat-jl) [![codecov](https://codecov.io/gh/FabianMaeurer/JuCat.jl/branch/master/graph/badge.svg?token=axGHAcozx5)](https://codecov.io/gh/FabianMaeurer/JuCat.jl)|

# JuCat.jl

JuCat is a package under development with the intention to provide a framework as well a examples for computations in the realm of categories.

## Installation

You need to have Julia installed. For reliable results Julia version at least 1.6 is required. To use JuCat
do the following:

```julia-repl
julia> import Pkg
julia> Pkg.add("https://github.com/FabianMaeurer/JuCat.jl")
```

## Usage

To use Jucat the structures from the [OSCAR-System](https://github.com/oscar-system/Oscar.jl) are required. Here a minimal usage Example.

```@repl
using JuCat, Oscar;
F = FiniteField(5)
G = symmetric_group(2)
X = gset(G,[1,2,3])
C = ConvolutionCategory(X,F)
simples(C)
groethendieck_ring(C)
```

## Acknowledgements

This project was started under supervision of [Prof. Ulrich Thiel](https://ulthiel.com/math/)  (University of Kaiserslautern). This work is a
contribution to the SFB-TRR 195 'Symbolic Tools in Mathematics and their
Application' of the German Research Foundation (DFG).


[docs-stable-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-url]: https://fabianmaeurer.github.io/JuCat.jl/

[ga-img]: https://github.com/fabianmaeurer/JuCat.jl/workflows/Run%20tests/badge.svg
[ga-url]: https://github.com/fabianmaeurer/JuCat.jl/actions?query=workflow%3A%22Run+tests%22
