# [Computing with fusion rings](@id computing-fusion-rings)

This page applies the preceding [Grothendieck-ring conventions](@ref grothendieck-rings).
It first constructs the ring of a category, then enters fusion rules without a
categorification, and finally compares a non-split category with its splitting
field.

## A split representation ring

Over $\mathbb F_3$, the cyclic group $C_2$ has the trivial and sign
representations. Their classes satisfy $[\varepsilon]^2=[\mathbb 1]$.

```@example grothsplit
using TensorCategories, Oscar
G = cyclic_group(2)
C = representation_category(GF(3), G)
S = simples(C)
R = split_grothendieck_ring(C)
unit_index = only(findall(X -> is_isomorphic(X, one(C))[1], S))
sign_index = only(setdiff(eachindex(S), [unit_index]))
u, e = R[unit_index], R[sign_index]
@assert u == one(R) && e*e == u
@assert involution(e) == e
e*e
show(stdout, MIME"text/plain"(), e*e); println() # hide
```

The basis of `R` follows the simple-object order of `C`. To obtain the class
of an object, pass its integer multiplicities to `R`:

```@example grothsplit
Y = S[sign_index] ⊕ S[sign_index]
y = R(ZZ.(coefficients(Y,S)))
@assert y == 2*e
coefficients(y)
```

Here `ZZ.(...)` converts each multiplicity to an OSCAR integer. A virtual
class such as `e-u` is also an element of `R`:

```@example grothsplit
@assert base_ring(R) == ZZ
@assert fpdim(e) == 1
@assert fpdim(R) == 2
fpdim.(basis(R))
```

Both Frobenius–Perron dimensions are $1$, and their squared sum is $2$.
The main ring operations are:

| Operation | Result |
|:---|:---|
| `basis(R)`, `R[i]` | Distinguished basis elements |
| `rank(R)` | Number of basis elements |
| `one(R)`, `zero(R)` | Ring identity and zero |
| `R(ZZ.(coeffs))` | Element with coefficient vector $(a_1,\ldots,a_r)$ supplied as `coeffs` |
| `coefficients(r)` | Coefficient vector of a ring element |
| `multiplication_table(R)` | Integer array `N[i,j,l]` |
| `involution(r)` | Dual class, when the involution is stored |
| `fpdim(r)` | Additive Frobenius--Perron dimension |

For semisimple rigid input, `split_grothendieck_ring` stores the involution
obtained from the duality permutation.

## Entering a ring without a category

`ZPlusRing` constructs a ring directly from its basis names, multiplication
table, and unit coefficient vector. Its aliases are `ℤ₊Ring` and `ℕRing`.
For the Fibonacci rule $t^2=1+t$:

```@example grothfibonacci
using TensorCategories, Oscar
N = zeros(Int,2,2,2)
N[1,1,1] = N[1,2,2] = N[2,1,2] = 1
N[2,2,1] = N[2,2,2] = 1
R = ZPlusRing(["1","t"], N, [1,0])
t = R[2]
@assert t*t == one(R)+t
d = fpdim(t)
@assert d > 0 && d^2 == 1+d
d
show(stdout, MIME"text/plain"(), d); println() # hide
```

Thus $\operatorname{FPdim}(t)=(1+\sqrt5)/2$. No associator or coefficient
field for a categorification was needed. Constructing a category with this
Grothendieck ring requires substantially more data.
The constructor converts the supplied table and unit coordinates to integers;
it does not itself certify nonnegativity, associativity, the unit equations, or
the based-ring identities. These are assumptions on directly entered data.

## A non-split representation ring

Over $\mathbb F_2$, the group $C_3$ has two irreducible representations: the
trivial representation and a two-dimensional representation $V$ with
endomorphism field $\mathbb F_4$. The category is semisimple, and

```math
\label{eq:rank-two-fusion-rule}
[V]^2=2[\mathbb 1]+[V].
```

Indeed, over $\mathbb F_4$ the representation $V$ splits as
$\chi\oplus\chi^{-1}$, so its square is
$2\cdot\mathbb 1\oplus\chi\oplus\chi^{-1}$.

```@example grothnonsplit
using TensorCategories, Oscar
G = cyclic_group(3)
C = representation_category(GF(2),G)
S = simples(C)
R = split_grothendieck_ring(C)
i = only(findall(Y -> int_dim(Y) == 2,S))
V, v = S[i], R[i]
@assert int_dim(End(V)) == 2
@assert v*v == 2*one(R)+v
@assert fpdim(v) == 2
@assert fpdim(C) == 3
@assert fpdim(R) == 5
v*v
show(stdout, MIME"text/plain"(), v*v); println() # hide
```

The coefficient $2$ is an integer multiplicity; it does not vanish in the
Grothendieck ring. The category's Frobenius--Perron dimension is
$1^2+2^2/2=3$, whereas the ring method's unweighted sum is $1^2+2^2=5$.

Over the splitting field, there are three one-dimensional simples:

```@example grothnonsplit
Cs = representation_category(GF(2,2),G)
Rs = split_grothendieck_ring(Cs)
@assert rank(Rs) == 3
@assert all(b -> fpdim(b) == 1,basis(Rs))
@assert fpdim(Rs) == 3
rank(Rs)
```

The field extension changes the simple basis and hence the Grothendieck ring;
its multiplication coefficients are integers over both fields.

Continue with [pivotal and spherical structures](@ref pivotal-braided).
