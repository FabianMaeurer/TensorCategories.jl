```@setup FC
using TensorCategories, Oscar
```

# Semisimple Ring Categories from 6j-Symbols

We provide a structure `RingCategory` for finite semisimple skeletal (multi)ring categories. 

They are constructed by setting the multiplication table for the monoidal product and associativity for the simple objects. Objects in finite semisimple ring categories are just vectors with multiplicities of the simple objects in the decomposition. Thus morphism are vectors of matrices of according size.

As an example we will construct the Ising fusion category. This is a fusion category with three simple objects ``1``, ``\chi`` and ``X``. The multiplication is given by ``\chi \otimes \chi = 1``, ``\chi \otimes X = X \otimes \chi`` and ``X \otimes X = 1 \oplus \chi``. There are precisely 3 non-trivial 6j-symbols ``\Phi_X^{\chi\,X\,\chi} = -1``, ``\Phi_{\chi}^{X\,\chi\,X} = -1`` und ``\Phi_X^{X\,X\,X} = \frac{1}{\sqrt 2}\begin{pmatrix} 1 & 1 \\ 1 & -1\end{pmatrix}``.

```@example FC
F = QQBar
I = RingCategory(F,["𝟙", "χ", "X"])

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
set_associator!(I, )