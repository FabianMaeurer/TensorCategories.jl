using TensorCategories, Oscar
#=------------------------------------------------
    Ising
------------------------------------------------=#
Ising_time = @elapsed begin

I = Ising()

𝟙, χ, X = simples(I)

S₁ = simple_subobjects(induction(𝟙))
S₂ = simple_subobjects(induction(χ))
S₃ = simple_subobjects(induction(X))

C = Center(I)

add_simple!(C, [S₁; S₂; S₃])

end
#=------------------------------------------------
    Subcategory of I2
------------------------------------------------=#
I26_time = @elapsed begin
    
B = I2subcategory(6)

Bs, Bsts, Bststs = simples(B)

S₁ = simple_subobjects(induction(Bs))
S₂ = simple_subobjects(induction(Bsts))
S₃ = simple_subobjects(induction(Bststs))

D = Center(B)

add_simple!(D, [S₁; S₂; S₃])

end