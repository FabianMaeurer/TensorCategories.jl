using TensorCategories, Oscar

G = small_group(30,3)
𝒞 = RepresentationCategory(G,GF(3,2))

# Construct a non-zero dimensional non-simple,
# indecomposable representation
x = decompose(regular_representation(𝒞))[1][1]
y = collect(values(eigenvalues(basis(End(x))[2])))[1]

# Tensor Powercategory genereated by y
𝒯 = TensorPowerCategory(y)

𝒮 = Semisimplification(𝒯)

S = simples(𝒮)

