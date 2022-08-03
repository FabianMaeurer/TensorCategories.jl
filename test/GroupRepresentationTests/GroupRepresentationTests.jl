G = symmetric_group(3)
F,a = FiniteField(23)

RepG = RepresentationCategory(G,F)
simple_objects = simples(RepG)
𝟙,σ,τ = simple_objects

@testset "Simple objects of Rep(S₃)" begin
    @test length(simple_objects) == 3
    @test dim.(simple_objects) == F.([1,1,2])
    @test dual.(simple_objects) == simple_objects
end
