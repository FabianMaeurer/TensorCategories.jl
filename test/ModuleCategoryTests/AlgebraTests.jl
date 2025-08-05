#=----------------------------------------------------------
    Test Algebras 
----------------------------------------------------------=#

V = graded_vector_spaces(QQ, symmetric_group(3))

# @testset "Algebras in VecG" begin
#     algs = algebra_structures(V[1,2])
#     algs2 = algebra_structures(V[1,4,5])
#     @test all(is_algebra.(algs))
#     @test all(is_algebra.(algs2))
# end

I = ising_category()

@testset "Algebras in Ising" begin 
    algs = algebra_structures(I[1,2])
    @test all(is_algebra.(algs))

    algs2 = algebra_structures(I[1,3])
    @test all(is_algebra.(algs2))
    @test !any(is_separable.(algs2))
end