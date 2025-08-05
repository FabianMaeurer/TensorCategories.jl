# Fusion categories from the AnyonWiki 

In [vercleyen2024lowrankmultiplicityfreefusioncategories](@cite) they were able to compute all multiplicity free fusion catogories up to rank 7. 
In cooperation with the author we were able to include the full database into the package. 

## The naming 

Please visit their website [AnyonWiki](https://anyonwiki.github.io) for the details on the naming scheme. The categories are given by 7 indices

```math 
(i,j,k,l,m,n,o)
```

where 

 - ``i`` is the rank, i.e the number of simple objects
 - ``j`` is the multiplicity
 - ``k`` is the number of _not_ self-dual simples 
 - ``l`` is the index of the fusion rules 
 - ``m`` is the index of the associator 
 - ``n`` is the index of the braiding - 0 if not braided 
 - ``o`` is the index of the pivotal structure 

## Access

The categories of the AnyonWiki can be accessed as instances of `SixJCategory` via the function 

```@docs
anyonwiki
```