```@setup FC
using TensorCategories, Oscar
```

# Semisimple Ring Categories from 6j-Symbols

We provide a structure `SixJCategory` for finite semisimple skeletal (multi)ring categories. 

They are constructed by setting the multiplication table for the monoidal product and associativity for the simple objects. Objects in finite semisimple ring categories are just vectors with multiplicities of the simple objects in the decomposition. Thus morphism are vectors of matrices of according size.

As an example we will construct the Ising fusion category. This is a fusion category with three simple objects ``1``, ``\chi`` and ``X``. The multiplication is given by ``\chi \otimes \chi = 1``, ``\chi \otimes X = X \otimes \chi = X`` and ``X \otimes X = 1 \oplus \chi``. There are precisely 3 non-trivial 6j-symbols ``\Phi_X^{\chi\,X\,\chi} = -1``, ``\Phi_{\chi}^{X\,\chi\,X} = -1`` und ``\Phi_X^{X\,X\,X} = \frac{1}{\sqrt 2}\begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix}``.

```@example FC
F = QQBarField()
I = six_j_category(F,["𝟙", "χ", "X"])

M = zeros(Int,3,3,3)

M[1,1,:] = [1,0,0]
M[1,2,:] = [0,1,0]
M[1,3,:] = [0,0,1]
M[2,1,:] = [0,1,0]
M[2,2,:] = [1,0,0]
M[2,3,:] = [0,0,1]
M[3,1,:] = [0,0,1]
M[3,2,:] = [0,0,1]
M[3,3,:] = [1,1,0]

set_tensor_product!(I,M)

set_associator!(I, 2,3,2,3, [-1])
set_associator!(I, 3,2,3,2, [-1])
set_associator!(I, 3,3,3,3, inv(sqrt(F(2))).*[1 1; 1 -1])
```

At the current state of development the unit object has to be set manually. Other structure like spherical structures are set as identities, which might or might not be correct. If you know your category to be spherical you can set the canonical spherical structure such that ``\dim X_i = \mathrm{fpdim}\,X_i`` by applying ``set_canonical_spherical(::SixJCategory)`` to it. It is also possible to set a braiding.

```@example FC
set_one!(I, [1,0,0])

#set one of four possible braidings 
α = root_of_unity(F,8)

braid = Array{MatElem,3}(undef, 3,3,3)

braid[1,1,:] = matrices(id(I[1]))
braid[1,2,:] = braid[2,1,:] = matrices(id(I[2]))
braid[2,2,:] = -1 .* matrices(id(I[1]))

braid[1,3,:] = braid[3,1,:] = matrices(id(I[3]))
braid[2,3,:] = braid[3,2,:] = α^2 .* matrices(id(I[3]))

braid[3,3,:] = matrices((α * id(I[1]) ⊕ (α^7 * id(I[2]))))

set_braiding!(I,braid)

# Check the pentagon equations for all simple objects in I
@show pentagon_axiom(I)
```


