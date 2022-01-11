# JuCat.jl

JuCat is a package under development with the intention to provide a framework as well a examples for computations in the realm of categories.

## Installation

You need to have Julia installed. For reliable results Julia version at least 1.6 is required. To use JuCat
do the following:

```julia
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
