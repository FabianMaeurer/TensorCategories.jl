# Test File for RingCategory with the Ising category

I = Ising()
a,b,c = simples(I)

@testset "Ising: Objects" begin
    X = a^2 ⊕ b ⊕ c^2
    Y = a ⊕ c^2
    @test X == RingCatObject(I, [2,1,2])
    @test X⊗Y == RingCatObject(I, [6,5,8])
end

@testset "Ising: Associativity" begin
    f = id(a)⊕ (-id(b))
    @test compose((id(b)⊗f)⊗id(c), associator(b,a⊕b,c)) == compose(associator(b,a⊕b,c), id(b)⊗(f⊗id(c)))
    @test pentagon_axiom(I)
    @test pentagon_axiom([simples(I).^2])
end

# @testset "Ising: Duals, Ev, Coev" begin
#     X = a^2 ⊕ b^2 ⊕ c
#     @test dual(X) == X
#     @test (id(X)⊗ev(X))∘associator(X,X,X)∘(coev(X)⊗id(X)) == id(X)
# end
