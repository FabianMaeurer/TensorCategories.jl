# JuCat.jl

JuCat is a package under development with the intention to provide a framework as well a examples for computations in the realm of categories.

##Installation

You need to have Julia installed. For reliable results Julia version at least 1.6 is required

```julia
julia> import Pkg
julia> Pkg.add("https://github.com/FabianMaeurer/JuCat.jl")
```

Then to use JuCat you can just type

```julia
julia> using JuCat, Oscar
```

to utilize JuCat and the required algebraic datatypes from the [OSCAR-System][2]


## Acknowledgements

This project was started under supervision of [Prof. Ulrich Thiel][1]  (University of Kaiserslautern). This work is a
contribution to the SFB-TRR 195 'Symbolic Tools in Mathematics and their
Application' of the German Research Foundation (DFG).


[1] = https://ulthiel.com/math/ 'Ulrich Thiel'
[2] = https://github.com/oscar-system/Oscar.jl'OSCAR'
