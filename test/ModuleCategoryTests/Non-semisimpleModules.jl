#=----------------------------------------------------------
    Test structures for non semisimple module categories 
----------------------------------------------------------=#

C = graded_vector_spaces(QQ, symmetric_group(3))

# take a non-separable_algebra
A = filter(!is_separable, algebra_structures(C[1]⊕C[2])) |> first

@testset "Non-semisimple modules" begin
    @test !is_seperable(A)
    @test !is_semisimple(endomorphism_ring(free_right_module(C[1,2],A)))    
end