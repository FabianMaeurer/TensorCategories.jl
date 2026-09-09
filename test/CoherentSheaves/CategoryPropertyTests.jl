@testset "Finite equivariant-sheaf categories" begin
    K = GF(2)
    G = cyclic_group(2)
    X = gset(G, (x, g) -> x, [1])

    S = coherent_sheaves(K, X)
    C = convolution_category(K, X)

    # The one-point trivial action identifies the stabilizer category with
    # Rep_K(G), which is finite but not semisimple in this characteristic.
    @test !is_semisimple(S) && !is_semisimple(C)
    @test is_finite(S) && is_locally_finite(S)
    @test is_finite(C) && is_locally_finite(C)
end
