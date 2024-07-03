#=----------------------------------------------------------
    Test UqSL2Representations 
----------------------------------------------------------=#

# Uq(𝔰𝔩₂) reps at q = √2
C = sl2_representations(quadratic_field(2)...)

@testset "Uq(sl2) representations" begin
    # test some Associators
    @test pentagon_axiom([C[1],C[2],C[3]])
    
    @test C[2] ⊗ C[3] == C[1,3,5]
end