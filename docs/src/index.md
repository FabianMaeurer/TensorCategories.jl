# TensorCategories.jl

TensorCategories.jl is a software package based on the programming language [Julia](https://julialang.org) and the open-source computer algebra system [Oscar.jl](https://github.com/oscar-system/Oscar.jl) for computations with tensor categories. 

## Installation

You need to have Julia installed. For reliable results Julia version at least 1.10 is required. To use TensorCategories
do the following:

```julia-repl
julia> import Pkg
julia> Pkg.add("TensorCategories.jl")
```

## Usage

To use TensorCategories the structures from the [OSCAR-System](https://github.com/oscar-system/Oscar.jl) are required. Here a minimal usage example.

```jldoctest 
using TensorCategories;
I = ising_category()
C = center(I)
S = smatrix(C)

# output 
1
```

## Acknowledgements

This project was started under supervision of [Prof. Ulrich Thiel](https://ulthiel.com/math/)  (University of Kaiserslautern). This work is a
contribution to the SFB-TRR 195 'Symbolic Tools in Mathematics and their
Application' of the German Research Foundation (DFG).
