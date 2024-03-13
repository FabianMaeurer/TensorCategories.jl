#=-------------------------------------------------
    Test sets for the induction functor
        𝒞 → 𝒵(𝒞) 
-------------------------------------------------=#

G    = symmetric_group(3)
H    = cyclic_group(3) 
F,ξ  = cyclotomic_field(3,"ξ")
c    = cyclic_group_3cocycle(H,F,ξ)  

VecG = GradedVectorSpaces(F,G)
VecH = GradedVectorSpaces(F,H,c)

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
RepG = RepresentationCategory(G,F)

# @testset "Group Representation Category" begin
#     S = simples(RepG)
#     induction_S = induction.(S)
#     for X ∈ induction_S
#         @test is_half_brading(object(X), half_braiding(X))
#     end
# end