#=----------------------------------------------------------
    Test center for graded Vector spaces
----------------------------------------------------------=#

V = graded_vector_spaces(QQ, symmetric_group(2))
V2 = graded_vector_spaces(QQ, symmetric_group(3))
Z = center(V)
Z2 = center(V2)

@testset "Compute Center" begin
    
    @test length(simples(Z)) == 4
    @test length(simples(Z2)) == 7

end

@testset "split Center" begin
    Z3,_ = split(Z2)
    @test length(simples(Z3)) == 8 
    
    @test all(is_central.(simples(Z3)))
    @test all([int_dim(End(s)) == 1 for s in simples(Z3)])
end