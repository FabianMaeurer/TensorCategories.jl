#=------------------------------------------------
    Test S-Matrix of Center of I₂(5) by comparison
    with Lusztigs and Malles
------------------------------------------------=#

# To speed up computation compute the equivalent centre category category
I = I2subcategory(5)
C = Center(I)
simples(C)

#Sort the simples in 𝒵(I₂(5)) to match Lusztigs
C.simples = C.simples[[]]

S = normalized_smatrix(C)

# Lusztigs S-matrix according to 
# https://projecteuclid.org/journals/duke-mathematical-journal/volume-73/issue-1/Exotic-Fourier-transform/10.1215/S0012-7094-94-07309-2.short

K = base_ring(I)
ξ = gen(K)
λ = -ξ^3 + ξ^2 + 1

Lusztigs_S = inv(sqrt(K(5))) * matrix(K,[  λ-1 λ   1    1;
                        λ   λ-1 -1   -1;
                        1   -1  λ    -λ+1;
                        1   -1  -λ+1 λ])

# Malle S-matrix according to 
# https://www.sciencedirect.com/science/article/pii/S0021869302006312

tuples = [(i,j) for i ∈ 0:5, j ∈ 0:5 if 0 < i < j < i+j < 5 || 0 == i < j < 5/2]

Malles_S = matrix(K, [inv(K(5)) * (ξ^(i*l+j*k) + ξ^(-i*l-j*k) - ξ^(i*k + j*l) - ξ^(-i*k-j*l)) for (i,j) ∈ tuples, (k,l) ∈ tuples])

@testset "𝒵(I₂(5)) S-matrix" begin
    @test "S == "
end