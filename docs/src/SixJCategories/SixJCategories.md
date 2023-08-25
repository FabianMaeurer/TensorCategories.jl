# [Fusion Categories from `6j`-Symbols] (@id 6j_categories)

In most literature concerning fusion categories they are characterized by 
so called ``6j``-symbols. Often only those and the fusion rules are provided or of interest. Thus we provide a datatype that allows to 
get a workable fusion category from provided ``6j`-symbols.

## ``6j``-Symbols

Let ``\mathcal C`` be a locally finite semisimple multitensor category.  Then, if ``\{X_i\mid i \in \mathcal I\}`` is a collection of the non-isomorphic simple objects, there is an equivalence of abelian categories 

```math
F \colon \mathcal C \cong \bigoplus\limits_{i \in \mathcal I} \mathrm{Vec}_k
```

given by 

```math
X \mapsto \mathrm{Hom}(X_i,X).
```

We define ``H_{ij}^k := /mathrm{Hom}(X_k, X_i\otimes X_j)`` to be the multiplicity spaces. Now considering the image of a tensor product $X_i \otimes X_j$ of two simple objects we obtain 

```math
X_i \otimes X_j \mapsto \bigoplus\limits_{k \in \mathcal I} H_{ij}^k
```

After fixing a natural isomorphism 

```math
(X_i \otimes X_j) \otimes X_k \cong X_i \otimes (X_j \otimes X_k)
```

we obtain morphisms 

```math
\bigoplus\limits_{k ∈ I} H_{ij}^k \xrightarrow
```


```@autodocs
Modules = [TensorCategories]
Pages = ["FusionCategory.jl"]
```