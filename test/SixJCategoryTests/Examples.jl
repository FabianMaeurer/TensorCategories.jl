#=----------------------------------------------------------
    Test 6j Examples
----------------------------------------------------------=#


@testset "6j-Categories" begin
    # test Ising
    @test pentagon_axiom(Ising())
    @test hexagon_axiom(Ising(cyclotomic_field(16)[1]))

    # test Tambara Yamagami for A = [2,2]
    @test pentagon_axiom(tambara_yamagami(2,2))

    # test HaagerupH3
    @test pentagon_axiom(haagerup_H3())

    # test Verline Categories
    @test pentagon_axiom(verlinde_category(5))

    # test I2 and I2subcategory
    @test pentagon_axiom(I2(4))
    @test pentagon_axiom(I2subcategory(5))
end