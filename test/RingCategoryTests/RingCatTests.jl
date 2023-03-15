#=------------------------------------------------
    Ising Category
------------------------------------------------=#

I = Ising()
a,b,c = simples(I)

@testset "Ising: CategoryObjects" begin
    X = a^2 ⊕ b ⊕ c^2
    Y = a ⊕ c^2
    @test X == RingCategoryObject(I, [2,1,2])
    @test X⊗Y == RingCategoryObject(I, [6,5,8])
end

@testset "Associator" begin
    @test pentagon_axiom(I)
    @test pentagon_axiom(c,c,c^2,c)
    @test pentagon_axiom(c^2,a⊕b,c⊕b,c)
end

#=------------------------------------------------
    I2 
------------------------------------------------=#

B = I2(5)

@testset "I2" begin
    @test pentagon_axiom(B)
    @test pentagon_axiom(B[2]⊕B[3], B[3]⊕B[4]⊕B[7], B[2]⊕B[1]⊕B[2], B[6]^2)
end

