#=----------------------------------------------------------
    Test sets for internal Module Categories 
----------------------------------------------------------=#

I = Ising()

𝟙,χ,X = simples(I)

A = separable_algebra_structures(𝟙 ⊕ χ)[1]

M1 = category_of_right_modules(A)

Funcs = category_of_bimodules(A)

@testset "Modules in Ising" begin
    @test length(simples(M1)) == 3
    @test pentagon_axiom(Funcs)
end

