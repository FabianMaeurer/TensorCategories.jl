```@setup ising_category
using TensorCategories
```

# [The Center Construction](@id center)

Let ``\mathcal C`` be a monoidal category. The __Drinfeld center__ of ``\mathcal C`` is
a category whose objects are tuples ``(X,\gamma)`` such that ``X \in \mathcal C`` and
``\{\gamma_Z \colon X \otimes Z \to Z \otimes X \mid Z \in \mathcal C\}`` is a natural isomorphism such that

```@raw html
<img src="https://i.upmath.me/svg/%5Cbegin%7Btikzcd%7D%0A%09%26%20(Y%20%5Cotimes%20X)%20%5Cotimes%20Z%20%5Car%5Br%2C%20%22a_%7BY%2CX%2CZ%7D%22%5D%20%26%20Y%20%5Cotimes%20(X%20%5Cotimes%20Z)%20%5Car%5Bdr%2C%20%22%5Cid_Y%20%5Cotimes%20%5Cgamma_X(Z)%22%5D%20%26%20%5C%5C%0A%20%20%20%20%20%20%20%20(X%20%5Cotimes%20Y)%20%5Cotimes%20Z%20%5Car%5Bur%2C%20%22%5Cgamma_X(Y)%20%5Cotimes%20%5Cid_Z%22%5D%20%5Car%5Bdr%2C%20%22a_%7BX%2CY%2CZ%7D%22%5D%20%26%20%26%20%26%20Y%20%5Cotimes%20(Z%20%5Cotimes%20X)%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%26%20X%20%5Cotimes%20(Y%20%5Cotimes%20Z)%20%5Car%5Br%2C%20%22%5Cgamma_X(Y%20%5Cotimes%20Z)%22%5D%20%26%20(Y%20%5Cotimes%20Z)%20%5Cotimes%20X%20%5Car%5Bur%2C%20%22a_%7BY%2CZ%2CX%7D%22%5D%20%20%26%0A%5Cend%7Btikzcd%7D" alt="\begin{tikzcd}
	&amp; (Y \otimes X) \otimes Z \ar[r, &quot;a_{Y,X,Z}&quot;] &amp; Y \otimes (X \otimes Z) \ar[dr, &quot;\id_Y \otimes \gamma_X(Z)&quot;] &amp; \\
        (X \otimes Y) \otimes Z \ar[ur, &quot;\gamma_X(Y) \otimes \id_Z&quot;] \ar[dr, &quot;a_{X,Y,Z}&quot;] &amp; &amp; &amp; Y \otimes (Z \otimes X) \\
            &amp; X \otimes (Y \otimes Z) \ar[r, &quot;\gamma_X(Y \otimes Z)&quot;] &amp; (Y \otimes Z) \otimes X \ar[ur, &quot;a_{Y,Z,X}&quot;]  &amp;
\end{tikzcd}" /></p>
```

commutes for all ``Y,Z \in \mathcal C`` and ``\gamma_{\mathbb 1} = \mathrm{id}_X``.

# Computing the Center

The Drinfeld center can be computed explicitely for reasonably small
fusion categories using the algorithm described in [maurer2024computing](@cite). Any fusion category implementing the corresponding 
interface is supported. 


## Example

```@example ising_category
I = ising_category()
C = center(I)
simples(C)
```

# Methods

```@autodocs
Modules = [TensorCategories]
Pages = ["Center.jl", "CenterChecks.jl", "Induction.jl"]
```