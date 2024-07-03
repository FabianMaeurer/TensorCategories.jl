#=-------------------------------------------------
    Test sets for the induction functor
        𝒞 → 𝒵(𝒞:\scr) 
-------------------------------------------------=#

G    = symmetric_group(3)
H    = cyclic_group(3) 
F,ξ  = cyclotomic_field(3,"ξ")
c    = cyclic_group_3cocycle(H,F,ξ)  

VecG = graded_vector_spaces(F,G)
VecH = graded_vector_spaces(F,H,c)

@testset "Graded Vector Spaces" begin
    Z = centralizer(VecG, VecG[2])

    simps = simples(Z)

    for s ∈ simps
        @test is_central(s)
    end
end


@testset "Graded Vector Spaces" begin
    Z = centralizer(VecH, VecH[2])

    simps = simples(Z)

    for s ∈ simps
        @test is_central(s) 
    end
end


F = GF(23)
RepG = representation_category(F,G)

# @testset "Group Representation Category" begin
#     S = simples(RepG)
#     induction_S = induction.(S)
#     for X ∈ induction_S
#         @test is_half_brading(object(X), half_braiding(X))
#     end
# end