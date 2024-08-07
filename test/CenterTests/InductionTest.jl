#=-------------------------------------------------
    Test sets for the induction functor
        𝒞 → 𝒵(𝒞) 
-------------------------------------------------=#

G    = symmetric_group(3)
H    = cyclic_group(3) 
F,ξ  = cyclotomic_field(3,"ξ")
c    = cyclic_group_3cocycle(H,F,ξ)  

VecG = graded_vector_spaces(F,G)
VecH = graded_vector_spaces(F,H,c)

@testset "Graded Vector Spaces" begin
    S = simples(VecG)
    induction_S = induction.(S)
    for X ∈ induction_S
        @test is_half_braiding(object(X), half_braiding(X))
    end
end

@testset "Twisted Graded Vector Spaces" begin
    S = simples(VecH)
    induction_S = induction.(S)
    for X ∈ induction_S
        @test is_half_braiding(object(X), half_braiding(X))
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