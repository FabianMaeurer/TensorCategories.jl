#=----------------------------------------------------------
    Test 6j Examples
----------------------------------------------------------=#

# These constructors encode established exact F-symbol examples.  A monoidal
# category requires the pentagon identity, and a proposed braiding additionally
# requires the two hexagon identities.  See P. Etingof, S. Gelaki,
# D. Nikshych, and V. Ostrik, Tensor Categories, AMS (2015), Definitions 2.2.8
# and 8.1.1.  These checks validate the supplied symbols; they do not classify
# all categorifications of the corresponding fusion rings.
@testset "6j-Categories" begin
    # test Ising
    @test pentagon_axiom(ising_category())
    @test hexagon_axiom(ising_category(cyclotomic_field(16)[1]))

    # test Tambara Yamagami for A = [2,2]
    @test pentagon_axiom(tambara_yamagami(2,2))

    # test HaagerupH3
    @test pentagon_axiom(haagerup_H3())

    # test Verlinde Categories
    @test pentagon_axiom(verlinde_category(5))

    # test I2 and I2subcategory
    @test pentagon_axiom(I2(4))
    @test pentagon_axiom(I2subcategory(5))
end

# Ardonne--Slingerland, arXiv:1004.5456v2, Appendix B, gives the F- and
# R-symbols for this rank-four SU(3)_3 subcategory, including multiplicity two.
@testset "SU(3)_3 rank-four subcategory" begin
    C = TensorCategories.su_3_3_subcategory()
    S = simples(C)
    X = S[2]
    K = base_ring(C)
    @test X ⊗ X == S[1] ⊕ X ⊕ X ⊕ S[3] ⊕ S[4]
    @test int_dim(Hom(X ⊗ X,X)) == 2
    @test dim.(S) == K.([1,3,1,1])
    @test dual(S[3]) == S[4]
    @test all(pentagon_axiom(W,X,Y,Z) for W in S,X in S,Y in S,Z in S)
    @test hexagon_axiom(C)
    @test is_pivotal(C;check=true)
    # Relabeling must transport the multiplicity-two intermediate channel in
    # 8⊗8. Outer-index permutation alone changes this F-matrix and destroys
    # coherence; Appendix B supplies the exact F/R oracle used here.
    perm = [2,4,1,3]
    oldF = copy(C.ass)
    olddims = dim.(S)
    TensorCategories.sort_simples!(C,perm)
    @test pentagon_axiom(C) && hexagon_axiom(C)
    @test dim.(simples(C)) == olddims[perm]
    TensorCategories.sort_simples!(C,invperm(perm))
    @test C.ass == oldF
    # Complete ten-index F dictionaries must retain the multiplicity-two
    # channel under a numerical CSV round trip, with unit label 3 after the
    # chosen permutation.
    TensorCategories.sort_simples!(C,perm)
    embedding = first(complex_embeddings(base_ring(C)))
    mktempdir() do dir
        f = joinpath(dir,"su3-F.csv")
        r = joinpath(dir,"su3-R.csv")
        numeric_symbols_to_csv(f,Dict(k => embedding(v,128)
                                      for (k,v) in TensorCategories.F_symbols(C)))
        numeric_symbols_to_csv(r,Dict(k => embedding(v,128)
                                      for (k,v) in R_symbols(C)))
        D = load_numeric_fusion_category(f,r,AcbField(64);
                                         pivotal=embedding.(C.pivotal))
        @test multiplication_table(D) == multiplication_table(C)
        @test D.one == C.one
        @test pentagon_axiom(D) && hexagon_axiom(D)
    end
end

# Rowell--Stong--Wang, arXiv:0712.1377v4, pp. 3--4 and Section 5.3.2.
@testset "Fibonacci and its Galois conjugate" begin
    K,_ = cyclotomic_field(5)
    C,D = fibonacci_category(K,1),fibonacci_category(K,2)
    d,e = dim(C[2]),dim(D[2])
    @test d^2 == d+1 && e^2 == e+1
    @test d != e && d*e == -1
    @test pentagon_axiom(C) && pentagon_axiom(D)
    @test is_pivotal(C;check=true) && is_pivotal(D;check=true)
end

# EGNO Exercise 4.7.16: pivotal structures on Vec_C3 are the three
# characters of C3, and only the trivial character is spherical.
@testset "Pivotal enumeration over a nonreal splitting field" begin
    K,z = cyclotomic_field(3)
    N = zeros(Int,3,3,3)
    for i in 0:2,j in 0:2
        N[i+1,j+1,mod(i+j,3)+1] = 1
    end
    C = six_j_category(K,N)
    set_one!(C,1)
    pivs = pivotal_structures(C)
    @test length(pivs) == 3
    @test all(p in pivs for p in [[K(1),z^i,z^(2i)] for i in 0:2])
    @test TensorCategories.spherical_structures(C) == [K.([1,1,1])]

    # The two embeddings of Q(zeta_3) conjugate the nontrivial pivotal
    # characters. Every component must be evaluated by the same exact field
    # homomorphism used for the rest of the category.
    set_pivotal!(C,K.([1,z,z^2]))
    images = Any[]
    for embedding in complex_embeddings(K)
        E = extension_of_scalars(C,QQBarField(),embedding)
        exact = TensorCategories._qqbar_embedding(embedding)
        @test E.pivotal == exact.(C.pivotal)
        @test is_pivotal(E;check=true)
        push!(images,E.pivotal[2])
    end
    @test length(unique(images)) == 2

    P = [2,1,3]
    D = six_j_category(K,N[P,P,P])
    set_one!(D,2)
    @test Set(Tuple(p[P]) for p in pivs) ==
          Set(Tuple(p) for p in pivotal_structures(D))

    U = six_j_category(QQ,ones(Int,1,1,1))
    set_one!(U,1)
    set_pivotal!(U,[QQ(2)])
    # Installing structural data is a declaration; the expensive coherence
    # equations are checked only on request.  On the rank-one fusion ring,
    # monoidality forces the sole pivotal component p to satisfy p^2=p, hence
    # the nonzero solution is p=1 rather than the deliberately supplied p=2.
    @test is_pivotal(U)
    @test !is_pivotal(U;check=true)
    @test pivotal_structures(U) == [[QQ(1)]]
    @test U.pivotal == [QQ(2)]
end

# EGNO Proposition 4.8.4 makes the pivotal trace a nonzero coordinate on the
# double-dual Hom space of a split simple. Skeletonization uses that coordinate
# in the same fusion bases as its associator and braiding.
@testset "Pivotal transport under skeletonization" begin
    G = cyclic_group(3)

    K,z = cyclotomic_field(3)
    C = graded_vector_spaces(K,G)
    g = first(gens(G))
    C.spherical = Dict(g^i => z^i for i in 0:2)
    D = skeletonize(C;check=true)
    @test D.pivotal == [C.spherical[x] for x in elements(G)]
    @test is_pivotal(D;check=true) && !is_spherical(D;check=true)

    F = GF(7)
    E = graded_vector_spaces(F,G)
    E.spherical = Dict(g^i => F(2)^i for i in 0:2)
    H = skeletonize(E;check=true)
    @test H.pivotal == [E.spherical[x] for x in elements(G)]
    @test all(!iszero,dim.(simples(H)))
    @test is_pivotal(H;check=true) && !is_spherical(H;check=true)
end

# EGNO Section 4.6: the Deligne square of pointed Vec_C2 has componentwise
# strict associator. Lazy symbols must survive construction and checking.
@testset "Lazy associators in a pointed Deligne product" begin
    N = zeros(Int,2,2,2)
    N[1,1,1] = N[1,2,2] = N[2,1,2] = N[2,2,1] = 1
    C = six_j_category(QQ,N)
    set_one!(C,1)
    D = tensor_product(C,C,String[],String[],true)
    @test rank(D) == 4
    @test pentagon_axiom(D)
end
