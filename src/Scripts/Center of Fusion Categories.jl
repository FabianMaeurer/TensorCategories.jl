using TensorCategories, Oscar
#=------------------------------------------------
    ising_category
------------------------------------------------=#
ising_category_time = @elapsed begin

I = ising_category()

𝟙, χ, X = simples(I)

S₁ = simple_subobjects(induction(𝟙))
S₂ = simple_subobjects(induction(χ))
S₃ = simple_subobjects(induction(X))

C = center(I)

add_simple!(C, [S₁; S₂; S₃])

end
#=------------------------------------------------
    Subcategory of I2
------------------------------------------------=#
I26_time = @elapsed begin
    
B = I2subcategory(5)

Bs, Bsts = simples(B)

S₁ = simple_subobjects(induction(Bs))
S₂ = simple_subobjects(induction(Bsts))


D = center(B)

add_simple!(D, [S₁; S₂; S₃])

end

I27_time = @elapsed begin
    
    B7 = I2subcategory(7)
    
    Bs, Bsts, Bststs = simples(B7)
    
    S₁ = simple_subobjects(induction(Bs))
    S₂ = simple_subobjects(induction(Bsts))
    S₃ = simple_subobjects(induction(Bststs))
    
    D2 = center(B)
    
    add_simple!(D2, [S₁; S₂; S₃])

    #s_matrix_7 = smatrix(D2)
    
end
    

#=------------------------------------------------
    Tambara-Yamagami Categories
------------------------------------------------=#

TY_time = @elapsed begin
    
    # the group ZZ/2×ZZ/2
    A = abelian_group(PcGroup,[2,2])

    TY = tambara_yamagami(A)
    
    S = vcat([simple_subobjects(induction(s)) for s ∈ simples(TY)]...)

    ZTY = center(TY)
    add_simple!(ZTY, S)
end

TY2_times = @elapsed begin
    
    # The group ZZ/2×ZZ/4
    A = abelian_group(PcGroup, [2,4])

    TY2 = tambara_yamagami(A)

    S = vcat([simple_subobjects(induction(s)) for s ∈ simples(TY2)]...)

    ZTY2 = center(TY2)
    add_simples(ZTY2, S)
end