# [Fiber functors and matrix realizations](@id fiber-functors)

A fiber functor on a ring category $\mathcal C$ over $k$ is a tensor
functor

```math
\label{eq:fiber-functor}
F\colon\mathcal C\longrightarrow\operatorname{Vec}_k,
```

hence an exact faithful $k$-linear functor together with coherent tensor and
unit isomorphisms [EGNO; Definition 5.1.1](@cite). The monoidal data are
essential. A faithful realization of morphisms by matrices need not be a fiber
functor.

For $\operatorname{Rep}_k(G)$, forgetting the group action gives the standard
fiber functor. The package stores its underlying vector spaces and matrices
inside each representation and intertwiner, so computations use this
realization even though it is not necessarily constructed as a separate Julia
functor value.

By contrast, suppose $\mathcal C$ is a finite split semisimple $k$-linear
category with chosen simple representatives $S_1,\ldots,S_r$. The assignment

```math
\label{eq:semisimple-linear-realization}
U(X)=\bigoplus_{i=1}^r\operatorname{Hom}_{\mathcal C}(S_i,X)
```

is a faithful exact linear realization after bases are chosen. It gives block
matrices for morphisms, but it has no automatic coherent identification of
$U(X)\otimes U(Y)$ with $U(X\otimes Y)$. It is therefore not generally a
fiber functor. Its vector-space dimension is the total multiplicity of simple
summands of $X$, rather than the categorical or Frobenius–Perron dimension.

In a non-split semisimple category,
$\operatorname{Hom}(S_i,X)$ is naturally a module over the division algebra
$\operatorname{End}(S_i)$. Forgetting that module structure can still give
vector-space coordinates, but it must not identify one copy of $S_i$ with a
single scalar coordinate.

TensorCategories.jl consequently treats `matrix(f)`, `tensor_product`, and the
associator as separate parts of the interface. A category can support all three
without possessing a fiber functor to ordinary vector spaces.

Continue with [skeletal fusion categories](@ref skeletal-fusion).
