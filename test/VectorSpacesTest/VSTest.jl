
m = matrix(QQ, [1 -1 1; -1 1 -1; 1 -1 1])
f = Morphism(m)

@testset "Objects" begin
    V = VectorSpaceObject(QQ, 2)
    W = VectorSpaceObject(QQ, ["v", "w", "x"])
    @test dim(V ⊕ W) == 5
    @test dim(V⊗W) == 6
end

@testset "(Co)Kernel" begin
    K,k = kernel(f)
    C,c = cokernel(f)
    @test iszero(matrix(f∘k))
    @test iszero(matrix(c∘f))
end
