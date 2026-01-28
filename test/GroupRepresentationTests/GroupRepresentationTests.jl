G = symmetric_group(3)
F = GF(23)

RepG = representation_category(F,G)
simple_objects = simples(RepG)
𝟙,σ,τ = simple_objects

@testset "Simple objects of Rep(S₃)" begin
    @test length(simple_objects) == 3
    @test dim.(simple_objects) == [1,1,2]
   # @test dual.(simple_objects) == simple_objects
end

@testset "Rep center" begin 
    S = simples(center(RepG))
    @test length(S) == 8
end
