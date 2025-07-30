# Tambara Yamagami Categories

Tambara and Yamagami classified a class of near-group categories of multiplicity one over algebraically closed fields in [TAMBARA1998692](@cite). Let ``A`` be an abelian group. Choose a square root ``\tau \in k`` of ``\vert A\vert`` and a non-degenerate symmetric bilinear form ``\chi \colon A \times A \to k^\times``. Then the Tambara-Yamagimi category ``TY(A,\chi,\tau)`` has objects ``A \cup \{m\}`` with fusion rules 

```math
\begin{align*}
    a ⊗ b = a+b,~~a \otimes m = m \otimes a = m,~~ m \otimes m = \sum\limits_{a \in G} a
\end{align*}
```
for ``a,b \in G``. The non-trivial associativity constraints are given by 

```math
\begin{align*}
    a_{a,m,b} = \chi(a,b)\mathrm{id}_m & a_{m,a,m} = \bigoplus\limits_{b\in A}\chi(a,b)\mathrm{id}_b & a_{m,m,m} = \bigoplus\limits_{a,b \in A} \frac{1}{\tau\chi(a,b)}\mathrm{id}_m.
\end{align*}
```

Those categories can be constructed with a generic symmetric bilinear form or with a custom bilinear form and over arbitrary fields.

```@docs
tambara_yamagami
```

Tambara-Yamagami categories are implemented as an instance of `SixJCategory` and hence all functionality follows from there.

## Ising Category

The Ising fusion category is a special example of a Tambara-Yamagami category with ``A = \mathbb Z_2``.  

```@docs
ising_category()
ising_category(::Ring)
```

# The Haagerup Subfactor

The fusion categories stemming from the Haagerup subfactor are a well known and important example of a fusion category. Details can be found in 
[wolf2021microscopic](@cite).

In the Morita equivalence class of the Haagerup subfactor lie three categories. We call them ``\mathcal H_1,\mathcal H_2`` and ``\mathcal H_3``. The third has multiplicity 1 and is also known as the Haagerup-Izumi category for ``\mathbb Z_3``. It has six simple objects and teh same fusion rules as ``\mathcal H_2``:

```math
\begin{array}{c||c|c|c|c|c|c}
        & \mathbb 1 & \alpha & \alpha^\ast & \rho & {}_{\alpha}\rho & {}_{\alpha^\ast}\rho \\ \hline \hline
        \mathbb 1 & \mathbb 1 & \alpha & \alpha^\ast & \rho & {}_{\alpha}\rho & {}_{\alpha^\ast}\rho \\ \hline
        \alpha & \alpha & \alpha^\ast & \mathbb 1 & {}_\alpha\rho & {}_{\alpha^\ast}\rho & \rho \\ \hline
        \alpha^\ast & \alpha^\ast & \mathbb 1 & \alpha & {}_{\alpha^\ast}\rho & \rho & {}{\alpha}\rho \\ \hline
        \rho & \rho & {}_{\alpha^\ast}\rho & {}_\alpha\rho & \mathbb 1 \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho & \alpha \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho & \alpha^\ast \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho \\ \hline 
        {}_\alpha\rho & {}_\alpha\rho & {}_{\alpha^\ast}\rho & \rho & \alpha \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho & \mathbb 1 \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho & \alpha^\ast  \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho \\ \hline
        {}_{\alpha^\ast}\rho & {}_{\alpha^\ast}\rho & \rho & {}_\alpha\rho & \alpha^\ast \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho & \alpha \oplus \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho & \mathbb 1 \oplus  \rho \oplus {}_\alpha\rho \oplus {}_{\alpha^\ast}\rho 
    \end{array}
```

All three can be accessed via

```@docs
haagerup_H1
haagerup_H2
haagerup_H3
```

# Fusion Categories From Truncated Hecke Categories

TODO: Explanation

```@docs
I2
I2subcategory
```

# Various Other Categories Given by ``6J``-Symbols

Here are some more examples to play around. 

## Categorifications by Vercleyen and Singerland

In [vercleyen2023low](@cite) they found a huge number of fusion rings and some explicit categorifications that are neither Tambara-Yamagami nor Haagerup-Izumi categories. 

### FR``{}_2^{82}``

They provide 97 different (but maybe equivalent) associators one can access.

```@docs
cat_fr_8122
```

### FR``{}_3^{94}``

We have a singel associator for this ring.

```@docs
cat_fr_9143
```