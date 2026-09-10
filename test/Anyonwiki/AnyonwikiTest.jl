# Deterministic full-suite integration tests for the AnyonWiki loaders and
# center pipeline. The database is a fixture, not a classification oracle.
# Fusion-category coherence is the complete pentagon identity (EGNO (2015),
# Section 2.2), and splitness means End(S)=K for every simple S; see
# Mäurer--Thiel, arXiv:2406.13438v2, Section 2.1.

@testset "AnyonWiki nonsplit center and skeletonization" begin
    # This rank-three entry has the representation ring of S3:
    #   ε²=1, εX=X, X²=1⊕ε⊕X,
    # but uses the second stored associator and no stored braiding.
    code = (3,1,0,2,2,0,1)
    C = anyonwiki(code...)

    N = zeros(Int,3,3,3)
    N[1,1,1] = N[1,2,2] = N[1,3,3] = 1
    N[2,1,2] = N[2,2,1] = N[2,3,3] = 1
    N[3,1,3] = N[3,2,3] = 1
    N[3,3,1] = N[3,3,2] = N[3,3,3] = 1
    @test multiplication_table(C) == N
    @test fpdim(C[3]) == QQBarField()(2)
    @test pentagon_axiom(C)

    # The center is semisimple but not split over Q(ζ3): five simples have
    # scalar endomorphisms and one has a three-dimensional division algebra.
    # This fixed profile prevents the test from silently bypassing splitting.
    Z = center(C)
    S = simples(Z)
    @test length(S) == 6
    @test sort(int_dim.(End.(S))) == [1,1,1,1,1,3]
    @test dim(Z) == dim(C)^2

    Zsplit, = split(Z)
    T = simples(Zsplit)
    @test length(T) == 8
    @test all(S -> int_dim(End(S)) == 1,T)
    @test base_ring(Zsplit) != base_ring(Z)

    # F-symbol computation must transport the split center to a skeletal
    # category with the same fusion rules. Check all 8^4 pentagons, not a
    # random sample.
    Zskel = six_j_category(Zsplit)
    @test multiplication_table(Zskel) == multiplication_table(Zsplit)
    @test pentagon_axiom(Zskel)
    @test is_pivotal(Zskel;check=true)
    @test is_spherical(Zskel;check=true)

    # A second associator for the same rank-three fusion ring must use the
    # identical simple order for its fusion bases and pivotal coordinates.
    C2 = anyonwiki(3,1,0,2,3,0,1)
    Z2split, = split(center(C2))
    Z2skel = six_j_category(Z2split)
    @test is_pivotal(Z2skel;check=true)
    @test is_spherical(Z2skel;check=true)
end

# Scalar conversion preserves the exact polynomial pentagon relations. QQBar
# and GF(17) are exact; AcbField uses enclosure overlap, so its result is
# numerical compatibility rather than a proof.
@testset "AnyonWiki scalar conversion" begin
    code = (3,1,0,1,1,1,1)
    @test pentagon_axiom(anyonwiki(QQBarField(),code...))
    @test pentagon_axiom(anyonwiki(GF(17),code...))
    @test pentagon_axiom(anyonwiki(AcbField(64),code...))
end

# Save/load tests use the same fixed Ising fixture and check every pentagon
# after the round trip. The scalar CSV subtest also fixes serialization order.
@testset "AnyonWiki saving and loading" begin
    code = (3,1,0,1,1,1,1)

    @testset "Numeric" begin
        mktempdir() do path
            K = AcbField(64)
            scalars = Dict([2,1] => K(3//2,5//4),
                           [1,2] => K(-2),
                           [1,1] => K(0,-1))
            scalar_file = joinpath(path,"scalar-symbols.csv")
            numeric_symbols_to_csv(scalar_file,scalars)
            labels = [parse.(Int,split(line,", ")[1:2])
                      for line in readlines(scalar_file)]
            @test labels == [[1,1],[1,2],[2,1]]
            @test numeric_symbols_from_csv(scalar_file,K) == scalars

            C = anyonwiki(code...)
            file = joinpath(path,"ising-numeric")
            numeric_symbols_to_csv(file,numeric_F_symbols(C;precision=64))
            D = load_numeric_fusion_category(file,AcbField(32))
            @test pentagon_axiom(D)
        end
    end

    @testset "Symbolic" begin
        mktempdir() do path
            C = anyonwiki(code...)
            save_fusion_category(C,path,"ising-symbolic")
            D = load_fusion_category(joinpath(path,"ising-symbolic"))
            @test multiplication_table(D) == multiplication_table(C)
            @test pentagon_axiom(D)
        end
    end
end

# These counts are a deterministic data-index contract, not a theorem about
# completeness of the AnyonWiki classification.
@testset "AnyonWiki index" begin
    @test length(anyonwiki_keys(5)) == 279
    @test length(anyonwiki_keys(5,"unitary")) == 56
end
