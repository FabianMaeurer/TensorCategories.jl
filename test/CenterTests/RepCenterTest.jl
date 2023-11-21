#=----------------------------------------------------------
    Test the Center for The Example Rep(G) 
----------------------------------------------------------=#

G = alternating_group(4)
F = GF(13,2)

Rep = RepresentationCategory(G,F)

D = Center(Rep)
S = simples(D)

@testset "Simples are central" begin
    for s in S
        @test is_central(s)
    end
end

@testset "Hom spaces are Central" begin
    H = End(S[3]⊗S[4])
    for f ∈ H
        @test is_central(f)
    end
end