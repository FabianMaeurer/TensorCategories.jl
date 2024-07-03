
# Monoidal Categories

A monoidal category is a quintuplet ``(\mathcal C, \otimes, \mathbb 1, a, \iota)`` where 

- ``\mathcal C`` is a category
- ``\otimes\colon \mathcal C \times \mathcal C \to \mathcal C``is a
  bifunctor
- ``a_{X,Y,Z} \colon (X \otimes Y) \otimes Z \to X \otimes (Y \otimes Z)`` is a natural transformation
- ``\iota \colon \mathbb 1 \otimes \mathbb 1 \to \mathbb 1`` is an isomorphism

such that 

```@raw html
<img src="https://i.upmath.me/svg/%5Cbegin%7Btikzcd%7D%0A%09%26%20%7B((W%20%5Cotimes%20X)%20%5Cotimes%20Y)%5Cotimes%20Z%7D%20%5C%5C%0A%09%7B(W%20%5Cotimes%20(X%20%5Cotimes%20Y))%20%5Cotimes%20Z%7D%20%26%26%20%7B(W%20%5Cotimes%20X)%20%5Cotimes%20(Y%5Cotimes%20Z)%7D%20%5C%5C%0A%09%7BW%20%5Cotimes%20((X%5Cotimes%20Y)%20%5Cotimes%20Z)%7D%20%26%26%20%7BW%20%5Cotimes%20(X%20%5Cotimes%20(Y%5Cotimes%20Z))%7D%0A%09%5Carrow%5B%22%7Ba_%7BW%2CX%2CY%7D%20%5Cotimes%20%5Cmathrm%7Bid%7D_Z%7D%22'%7Bpos%3D1%7D%2C%20shorten%20%3C%3D4pt%2C%20from%3D1-2%2C%20to%3D2-1%5D%0A%09%5Carrow%5B%22%7Ba_%7BW%2CX%5Cotimes%20Y%2CZ%7D%7D%22'%2C%20shift%20right%3D5%2C%20draw%3Dnone%2C%20from%3D2-1%2C%20to%3D3-1%5D%0A%09%5Carrow%5Bfrom%3D2-1%2C%20to%3D3-1%5D%0A%09%5Carrow%5B%22%7B%5Cmathrm%7Bid%7D%5Cotimes%20a_%7BX%2CY%2CZ%7D%7D%22'%2C%20from%3D3-1%2C%20to%3D3-3%5D%0A%09%5Carrow%5B%22%7Ba_%7BW%5Cotimes%20X%2C%20Y%2CZ%7D%7D%22%7Bpos%3D0.8%7D%2C%20from%3D1-2%2C%20to%3D2-3%5D%0A%09%5Carrow%5B%22%7Ba_%7BW%2CX%2CY%5Cotimes%20Z%7D%7D%22%2C%20shift%20left%3D5%2C%20draw%3Dnone%2C%20from%3D2-3%2C%20to%3D3-3%5D%0A%09%5Carrow%5Bfrom%3D2-3%2C%20to%3D3-3%5D%0A%5Cend%7Btikzcd%7D" alt="\begin{tikzcd}
	&amp; {((W \otimes X) \otimes Y)\otimes Z} \\
	{(W \otimes (X \otimes Y)) \otimes Z} &amp;&amp; {(W \otimes X) \otimes (Y\otimes Z)} \\
	{W \otimes ((X\otimes Y) \otimes Z)} &amp;&amp; {W \otimes (X \otimes (Y\otimes Z))}
	\arrow[&quot;{a_{W,X,Y} \otimes \mathrm{id}_Z}&quot;'{pos=1}, shorten &lt;=4pt, from=1-2, to=2-1]
	\arrow[&quot;{a_{W,X\otimes Y,Z}}&quot;', shift right=5, draw=none, from=2-1, to=3-1]
	\arrow[from=2-1, to=3-1]
	\arrow[&quot;{\mathrm{id}\otimes a_{X,Y,Z}}&quot;', from=3-1, to=3-3]
	\arrow[&quot;{a_{W\otimes X, Y,Z}}&quot;{pos=0.8}, from=1-2, to=2-3]
	\arrow[&quot;{a_{W,X,Y\otimes Z}}&quot;, shift left=5, draw=none, from=2-3, to=3-3]
	\arrow[from=2-3, to=3-3]
\end{tikzcd}" />
```


commutes for all objects ``W,X,Y,Z`` in ``\mathcal C`` and

```math
\begin{align*}
	L_{\mathbb 1} \colon & X \to \mathbb 1 \otimes X \\
	R_{\mathbb 1} \colon & X \to X \otimes \mathbb 1
\end{align*}
```

are autoequivalences.

## Conventions and Restrictions

At the current state all monoidal categories are assumed to satisfy ``X \otimes \mathbb 1 \cong X \cong \mathbb 1 \otimes X`` and ``\iota = \mathrm{id}_{\mathbb 1}``.

But building non-strict monoidal categories is explicitly encouraged, as this support is a strength of our Package. 

## Monoidal Categories

Following the definition we need the following methods.

- `tensor_product(X::YourObject...)::YourObject` returning the monoidal product. You can access your method by invoking the infix operate `⊗`.
- `tensor_product(f::YourMorphism...)::YourMorphism` returning the the monoidal product of morphisms. You can access your method by invoking the infix operate `⊗`.
- `one(C::YourCategory)::YourObject` returning the monoidal unit.
- `associator(X::YourObject, Y::YourObject, Z::YourObject)`.

## Rigidity

Whenever there are objects which admit duals it is feasible to acces them.

- `left_dual(X::YourObject)::YourObject` return the left dual ``X^\ast``.
- `right_dual(X::YourObject)::YourObject` return the right dual ``{}^\ast X``.
- `ev(X::YourObject)::YourMorphism` return the evaluation morphism ``\mathrm{ev}_X\colon X^\ast \otimes X \to \mathbb 1``.
- `coev(X::YourObject)::YourMorphism` return the coevaluation morphism ``\mathrm{coev}_X\colon \mathbb 1 \to X\otimes X^\ast``. 

This allows to generically compute 

```@docs 
left_dual(::Morphism)
```

Note that `dual` will always call `left_dual`.

## Checks

To verify for oneself the pentagon and hexagon axioms can be checked.

```@docs
pentagon_axiom
hexagon_axiom
```





