# Generic interface regressions. References and counterexamples are recorded
# with the testsets that use them.

struct UndeclaredAuditCategory <: Category end
@attributes mutable struct DeclaredAuditCategory <: Category end
struct MultifusionShortcutAuditCategory <: Category end
TensorCategories.is_multifusion(::MultifusionShortcutAuditCategory) = true
TensorCategories.is_fusion(::MultifusionShortcutAuditCategory) =
    error("the fusion predicate must not be needed here")
struct SpuriousSMatrixAuditCategory <: Category end
TensorCategories.smatrix(::SpuriousSMatrixAuditCategory) = identity_matrix(QQ,1)
struct BrokenModularAuditCategory <: Category end
TensorCategories.is_fusion(::BrokenModularAuditCategory) = true
TensorCategories.is_braided(::BrokenModularAuditCategory) = true
TensorCategories.is_spherical(::BrokenModularAuditCategory) = true
TensorCategories.smatrix(::BrokenModularAuditCategory) = error("broken S-matrix")

# The center-specific SixJ implementation reuses precomputed central Hom bases.
# Keyword dispatch must retain that specialization; falling back to `Category`
# would reconstruct half-braidings inside every basis composition.
@testset "Center SixJ keyword dispatch" begin
    kwtypes = Tuple{
        NamedTuple{(:homs,),Tuple{Dict{NTuple{3,Int},HomSpace}}},
        typeof(TensorCategories.six_j_symbols),
        CenterCategory,
        Vector{CenterObject},
    }
    @test which(Core.kwcall,kwtypes).sig.parameters[4] === CenterCategory
end

# Structural predicates record axioms supplied by an implementation. The mere
# presence of a generic method does not establish an axiom. The fusion/weak
# fusion convention is Mäurer--Thiel, arXiv:2406.13438v2, Section 2.1; weak
# multifusion and multifusion are the corresponding variants without a simple
# unit object.
@testset "Declared categorical structures" begin
    C = UndeclaredAuditCategory()
    # Previously these two queries called each other indefinitely.
    @test !is_additive(C) && !is_abelian(C) && !is_linear(C) && !is_monoidal(C)
    @test !is_finite(C) && !is_locally_finite(C)

    # EGNO (2015), Section 8.13 requires a braided spherical fusion category:
    # an arbitrary nonsingular matrix does not establish modularity.
    @test !is_modular(SpuriousSMatrixAuditCategory())
    # Once all hypotheses are asserted, an implementation failure in S must
    # remain visible rather than being converted into the mathematical answer.
    @test_throws ErrorException is_modular(BrokenModularAuditCategory())

    D = DeclaredAuditCategory(Dict{Symbol, Any}(
        :additive => true,
        :rigid => true,
        :spherical => false,
        :krull_schmidt => true,
    ))
    # Additivity does not imply k-linearity, and every property reads its own
    # declaration rather than the unrelated spherical flag.
    @test is_additive(D) && !is_linear(D)
    @test TensorCategories.is_rigid(D) && is_krull_schmidt(D)

    K = DeclaredAuditCategory(Dict{Symbol, Any}(:krull_schmidt => true))
    @test is_krull_schmidt(K) && is_additive(K)
    @test !is_abelian(K) && !is_linear(K)

    S = DeclaredAuditCategory(Dict{Symbol, Any}(:spherical => true))
    # Spherical structure is rigid monoidal structure, but does not by itself
    # assert the finiteness hypotheses used for Krull--Schmidt decomposition.
    @test TensorCategories.is_rigid(S) && is_monoidal(S)
    @test !is_krull_schmidt(S)

    SS = DeclaredAuditCategory(Dict{Symbol, Any}(:semisimple => true))
    @test is_semisimple(SS) && is_abelian(SS) && is_additive(SS)
    @test !is_linear(SS) && !is_locally_finite(SS)
    @test !is_krull_schmidt(SS)

    L = DeclaredAuditCategory(Dict{Symbol, Any}(:locally_finite => true))
    @test is_locally_finite(L) && is_linear(L) && is_abelian(L) && is_additive(L)
    @test is_krull_schmidt(L)
    @test !is_finite(L) && !is_semisimple(L)

    Fin = DeclaredAuditCategory(Dict{Symbol, Any}(:finite => true))
    @test is_finite(Fin) && is_locally_finite(Fin)
    @test is_linear(Fin) && is_abelian(Fin) && is_additive(Fin)
    @test is_krull_schmidt(Fin) && !is_semisimple(Fin)

    MR = DeclaredAuditCategory(Dict{Symbol, Any}(:multiring => true))
    @test is_multiring(MR) && is_locally_finite(MR)
    @test is_linear(MR) && is_abelian(MR) && is_additive(MR)
    @test is_krull_schmidt(MR) && !is_finite(MR)

    W = DeclaredAuditCategory(Dict{Symbol, Any}(:weak_fusion => true))
    @test is_weak_fusion(W) && is_weak_multifusion(W) && is_semisimple(W)
    # Weak fusion categories retain the simple unit of a tensor category; only
    # splitness is dropped. Hence weak fusion implies tensor and ring.
    @test is_tensor(W) && is_multitensor(W) && is_ring(W) && is_multiring(W)
    @test !is_fusion(W) && !is_multifusion(W)

    WM = DeclaredAuditCategory(Dict{Symbol, Any}(:weak_multifusion => true))
    # A weak multifusion category may have a non-simple unit, so it implies the
    # multi-variants but not the tensor/ring variants.
    @test is_weak_multifusion(WM) && is_multitensor(WM) && is_multiring(WM)
    @test !is_weak_fusion(WM) && !is_tensor(WM) && !is_ring(WM)

    M = DeclaredAuditCategory(Dict{Symbol, Any}(:multifusion => true))
    @test is_multifusion(M) && is_weak_multifusion(M) && is_semisimple(M)
    @test !is_fusion(M) && !is_weak_fusion(M) && is_split_semisimple(M)

    # Multifusion already implies weak multifusion and semisimplicity. Checking
    # the stronger fusion condition caused recursive centralizer decomposition.
    M0 = MultifusionShortcutAuditCategory()
    @test is_weak_multifusion(M0) && is_semisimple(M0)

    F = DeclaredAuditCategory(Dict{Symbol, Any}(:fusion => true))
    @test is_fusion(F) && is_multifusion(F)
    @test is_weak_fusion(F) && is_weak_multifusion(F)
    @test is_split_semisimple(F)

    V = vector_spaces(QQ)
    @test is_braided(V) && is_spherical(V)
    @test is_additive(id(V)) && is_linear(id(V))

    R = representation_category(GF(5), cyclic_group(5))
    @test is_finite(R) && is_locally_finite(R)
    @test is_linear(R) && is_abelian(R) && is_additive(R)
    @test is_krull_schmidt(R) && !is_semisimple(R)

    A = ArrowCategory(V)
    @test is_finite(A) && is_locally_finite(A)
    @test is_linear(A) && is_abelian(A) && is_krull_schmidt(A)

    P = product_category(R, R)
    @test is_finite(P) && is_locally_finite(P)
    @test is_linear(P) && is_abelian(P) && is_krull_schmidt(P)
    @test is_finite(op(R)) && is_locally_finite(op(R))

    # A product over unrelated coefficient fields is still abelian, but this
    # model does not equip it with a common linear structure.
    mixed = product_category(vector_spaces(QQ), vector_spaces(GF(5)))
    @test is_abelian(mixed) && !is_linear(mixed) && !is_locally_finite(mixed)
    @test !is_multiring(mixed) && !is_multitensor(mixed)
end

# The categorical product has objects and arrows componentwise, with
# Hom((X_i),(Y_i)) = direct_sum_i Hom(X_i,Y_i).  The tensor product is the
# bifunctor C x C -> C on both objects and arrows; see EGNO (2015),
# Sections 1.1 and 2.1.
@testset "Product category and tensor bifunctor" begin
    V = VectorSpaces(QQ)
    X, Y = VectorSpaceObject(V,2), VectorSpaceObject(V,3)
    f = morphism(X, X, matrix(QQ, [1 2; 0 1]))
    g = morphism(Y, Y, matrix(QQ, [1 0 0; 0 2 0; 0 0 3]))

    P = product_category([V,V])
    @test P == product_category(V,V) == V × V

    XY = product_object(X,Y)
    fg = product_morphism(f,g)
    @test parent(XY) == P
    @test domain(fg) == XY && codomain(fg) == XY
    @test int_dim(Hom(XY,XY)) == 4 + 9

    T = TensorCategories.TensorProductFunctor(V)
    @test domain(T) == P && codomain(T) == V
    @test T(XY) == X ⊗ Y
    @test T(fg) == f ⊗ g
    @test T(id(XY)) == id(T(XY))
end

# In a skeletal semisimple category, a component vector records multiplicities
# relative to the chosen simple objects of one parent; see EGNO (2015),
# Sections 1.1 and 4.9. Equal coordinates in different categories do not define
# equal objects or mixed-category morphisms.
@testset "Parent-sensitive skeletal operations" begin
    F = GF(11)
    M = zeros(Int, 2, 2, 2)
    for i in 1:2, j in 1:2
        M[i, j, mod(i + j - 2, 2) + 1] = 1
    end
    C = six_j_category(F, M, ["1", "g"])
    set_one!(C, [1, 0])
    @test is_finite(C) && is_locally_finite(C)
    D = fibonacci_category(F)
    g, X = simples(C)[2], simples(D)[2]

    # Both have coordinates [0,1], but g⊗g=1 whereas X⊗X=1⊕X.
    @test g != X && !is_isomorphic(g, X)[1]
    @test_throws ArgumentError Hom(g, X)
    @test_throws ArgumentError zero_morphism(g, X)
    @test_throws ArgumentError morphism(g, X, matrices(id(g)))
    @test_throws ArgumentError g ⊗ X
    @test_throws ArgumentError direct_sum(g, X)
    @test_throws ArgumentError g ⊕ X
    @test_throws ArgumentError express_in_basis(id(g), End(X))

    g2 = simples(C)[2]
    ok, f = is_isomorphic(g, g2)
    @test ok && domain(f) === g && codomain(f) === g2
    @test inv(f) ∘ f == id(g)
    @test Dict(g => 7)[g2] == 7

    C2 = six_j_category(F, M, ["1", "g"])
    set_one!(C2, [1, 0])
    g3 = simples(C2)[2]
    # Structurally identical presentations still require explicit transport.
    @test g != g3
    @test parent(extension_of_scalars(g, F, C2)) === C2

    G = cyclic_group(5)
    R11, R13 = representation_category(F, G), representation_category(GF(13), G)
    # Zero representations retain their group and coefficient field.
    @test zero(R11) != zero(R13)
    @test !is_isomorphic(zero(R11), zero(R13))[1]
end

# Finite coproducts in an additive category are biproducts; see EGNO (2015),
# Section 1.2. The symbolic operator returns the object, while `coproduct`
# also returns the universal injections.
@testset "Coproduct operator" begin
    F = GF(5)
    C = six_j_category(F, ones(Int, 1, 1, 1), ["1"])
    set_one!(C, [1])
    U = one(C)
    @test ∐(U, U) == U ⊕ U

    V, inc = coproduct(U, U)
    h = morphism(V, U, [matrix(F, 2, 1, [2, 3])])
    # The copairing [2*id,3*id] is characterized by its restrictions.
    @test h ∘ inc[1] == F(2) * id(U)
    @test h ∘ inc[2] == F(3) * id(U)
end

# A scalar multiple is an equality M=cN over the common matrix base ring.
# The coefficient must be selected using a nonzero entry of N, and conversion
# of categorical scalars must land in the ring requested by the caller.
@testset "Scalar-multiple fallbacks" begin
    F = GF(5)
    N = matrix(F, [0 1; 0 0])
    I = identity_matrix(F, 2)
    # E12 is not scalar; its nonzero off-diagonal entry faces a zero entry of I.
    @test TensorCategories.is_scalar_multiple(N, I) == (false, nothing)
    @test TensorCategories.is_scalar_multiple(zero_matrix(F, 2, 2), I) ==
          (true, F(0))
    @test TensorCategories.is_scalar_multiple(F(3) * N, N) == (true, F(3))
    # Matrix equality, hence scalar-multiple equality, requires equal shapes.
    @test TensorCategories.is_scalar_multiple(I, identity_matrix(F, 1)) ==
          (false, nothing)

    C = six_j_category(F, ones(Int, 1, 1, 1), ["1"])
    set_one!(C, [1])
    L = GF(5, 2)
    # The previous fallback returned an F5 element despite being called as L(f).
    c = L(F(3) * id(one(C)))
    @test c == L(3) && parent(c) === L
end

# Scalar extension applies one chosen field embedding to objects and every
# Hom-space basis morphism. For finite fields these embeddings are controlled
# by Frobenius; see K. Conrad, Finite Fields, Section 5,
# https://kconrad.math.uconn.edu/blurbs/galoistheory/finitefields.pdf.
@testset "Finite-field Hom scalar extension" begin
    F = GF(5)
    C = six_j_category(F, ones(Int, 1, 1, 1), ["1"])
    set_one!(C, [1])
    TensorCategories.set_twist!(C, [F(2)])
    U = one(C)
    L = GF(5, 2)
    CL = extension_of_scalars(C, L)
    # Optional metadata is transported only when present, along the same map.
    @test getfield(CL, :twist) == [L(2)] && !isdefined(CL, :name)
    HL = extension_of_scalars(Hom(U, U), L, CL)
    @test int_dim(HL) == 1
    @test only(basis(HL)) == id(domain(HL))
    @test parent(domain(HL)) === CL

    F9, L81 = GF(3, 2), GF(3, 4)
    D = six_j_category(F9, ones(Int, 1, 1, 1), ["1"])
    set_one!(D, [1])
    e = Oscar.embed(F9, L81)
    DL = extension_of_scalars(D, L81; embedding=e)
    a = gen(F9)
    H = HomSpace(one(D), one(D), [a * id(one(D))])
    H1 = extension_of_scalars(H, L81, DL; embedding=e)
    H2 = extension_of_scalars(H, L81, DL; embedding=x -> e(x)^3)
    # Frobenius x↦x³ is the nontrivial F3-automorphism of F9.
    @test matrix(only(basis(H1)))[1,1] == e(a)
    @test matrix(only(basis(H2)))[1,1] == e(a)^3
end

# Vec_G is the direct sum of its homogeneous components, and tensor degrees
# multiply in the order gh; see EGNO (2015), Section 2.3.
@testset "Graded multiplicities and tensor order" begin
    G = cyclic_group(2)
    e, g = one(G), only(gens(G))
    C = graded_vector_spaces(QQ, G)
    X, Y, Z = C[e,e,g], C[e,g,g], C[g,e,e]
    # Equal support and total dimension do not imply equal graded dimensions.
    @test is_isomorphic(X, Y) == (false, nothing)
    ok, f = is_isomorphic(X, Z)
    @test ok && inv(f) ∘ f == id(X) && f ∘ inv(f) == id(Z)

    # A graded map has no nonzero component between unequal degrees.
    @test_throws ErrorException morphism(C[e], C[g], matrix(QQ, 1, 1, [1]))
    @test_throws ArgumentError morphism(C[e], C[e], identity_matrix(GF(5), 1))
    @test_throws ArgumentError morphism(C[e], C[e], identity_matrix(QQ, 2))

    H = symmetric_group(3)
    D = graded_vector_spaces(QQ, H)
    a, b = gens(H)[1:2]
    @test a*b != b*a
    @test TensorCategories.grading(D[a] ⊗ D[b]) == [a*b]
    @test_throws ArgumentError C[e] ⊗ D[a]
end

# Hom spaces are vector spaces (EGNO (2015), Section 1.2): zero has the unique
# empty coordinate vector in the zero subspace, while a nonzero vector is not
# in that span. A one-factor tensor product is the factor itself.
@testset "Zero coordinates and unary tensor products" begin
    V = VectorSpaces(QQ)
    U = one(V)
    empty_basis = VectorSpaceMorphism[]
    @test express_in_basis(zero_morphism(U,U), empty_basis) == QQFieldElem[]
    @test_throws ArgumentError express_in_basis(id(U), empty_basis)
    @test tensor_product(U) === U
    @test_throws ArgumentError tensor_product()

    X, Y = VectorSpaceObject(V,2), VectorSpaceObject(V,3)
    @test isempty(basis([zero_morphism(X,Y)]))
    @test_throws ArgumentError express_in_basis(
        zero_morphism(X,Y), [zero_morphism(Y,X)])
end

# Row reduction in a Hom space must preserve the original source and target.
# Rectangular matrices detect both accidental transposition and flattening in
# the wrong order; compare the vector-space structure in EGNO (2015), §1.2.
@testset "Rectangular Hom-span coordinates" begin
    V = VectorSpaces(QQ)
    X, Y = VectorSpaceObject(V,2), VectorSpaceObject(V,3)
    f = morphism(X,Y,matrix(QQ,2,3,[1,2,0,3,0,4]))
    g = morphism(X,Y,matrix(QQ,2,3,[0,1,2,0,3,0]))
    B = basis([f,g,f+g])
    @test length(B) == 2
    @test all(domain(b) == X && codomain(b) == Y for b in B)
    # The independent output spans each original generator exactly.
    for h in (f,g)
        c = express_in_basis(h,B)
        @test sum((a*b for (a,b) in zip(c,B));
                  init=zero_morphism(X,Y)) == h
    end
end

# A linear map Hom(X,Y)→Hom(X,Z) acts on row-coordinate vectors. Its
# represented matrix must preserve the Hom endpoints and coefficient field.
@testset "Linear maps between Hom spaces" begin
    V = VectorSpaces(QQ)
    X, Y, Z = (VectorSpaceObject(V,n) for n in 1:3)
    H, L = Hom(X,Y), Hom(X,Z)
    A = matrix(QQ,2,3,[1,2,3,4,5,6])
    T = morphism(H,L,A)
    for i in 1:2
        expected = sum((A[i,j]*b for (j,b) in enumerate(basis(L)));
                       init=zero_morphism(X,Z))
        @test T(basis(H)[i]) == expected
    end
    @test_throws ArgumentError T(only(basis(Hom(Y,X))))
    @test_throws ArgumentError morphism(H,L,matrix(GF(5),2,3,[1,2,3,4,0,1]))
    ordinary = morphism(VectorSpaceObject(V,1),VectorSpaceObject(V,1),
                        identity_matrix(QQ,1))
    @test_throws ArgumentError ordinary(id(X))
end

# The opposite category reverses arrows and composition, interchanges products
# with coproducts, and uses the inverse associator; see EGNO (2015), §1.1.
@testset "Opposite-category variance" begin
    V = VectorSpaces(QQ)
    U, A, B = [VectorSpaceObject(V,n) for n in 1:3]
    f = morphism(U,A,matrix(QQ,1,2,[1,2]))
    g = morphism(A,B,matrix(QQ,2,3,[1,0,2,0,1,3]))
    O = opposite_category(V)
    of, og = O(f), O(g)

    @test all(domain(h) == O(A) && codomain(h) == O(B)
              for h in basis(Hom(O(A),O(B))))
    @test of ∘ og == O(g ∘ f)
    @test opposite_morphism(of) == f
    @test opposite_object(O(A)) == A

    Z, inc, proj = direct_sum(O(U),O(A))
    @test proj[1] ∘ inc[1] == id(O(U))
    @test proj[2] ∘ inc[2] == id(O(A))
    @test inc[1] ∘ proj[1] + inc[2] ∘ proj[2] == id(Z)
    P, ps = product(O(U),O(A))
    Q, js = coproduct(O(U),O(A))
    @test domain(ps[1]) == P && codomain(ps[1]) == O(U)
    @test domain(js[2]) == O(A) && codomain(js[2]) == Q

    a = associator(O(U),O(A),O(B))
    @test domain(a) == (O(U)⊗O(A))⊗O(B)
    @test codomain(a) == O(U)⊗(O(A)⊗O(B))
    @test O(f) ⊗ O(g) == O(f ⊗ g)

    K, k = kernel(of)
    C, c = cokernel(of)
    @test of ∘ k == zero_morphism(K,codomain(of))
    @test c ∘ of == zero_morphism(domain(of),C)
end

# Hom(X,-) is covariant by postcomposition, while Hom(-,X) is a functor on
# C^op by precomposition; see EGNO (2015), §1.1.
@testset "Hom-functor variance" begin
    V = VectorSpaces(QQ)
    U, A, B = [VectorSpaceObject(V,n) for n in 1:3]
    f = morphism(U,A,matrix(QQ,1,2,[1,2]))
    g = morphism(A,B,matrix(QQ,2,3,[1,0,2,0,1,3]))
    H = Hom(U,:)
    @test H(g)(f) == g ∘ f
    @test H(g ∘ f) == H(g) ∘ H(f)
    @test H(id(A)) == id(H(A))

    O = opposite_category(V)
    of, og = O(f), O(g)
    K = Hom(:,B)
    @test K(of)(g) == g ∘ f
    @test K(of ∘ og) == K(of) ∘ K(og)
    @test domain(K) == O && is_additive(K) && is_linear(K)
end

function audit_pointed_c3()
    K,z = cyclotomic_field(3)
    N = zeros(Int,3,3,3)
    for i in 0:2,j in 0:2
        N[i+1,j+1,mod(i+j,3)+1] = 1
    end
    C = six_j_category(K,N,["1","g","g^2"])
    set_one!(C,1)
    R = [N[i,j,k] == 1 ? matrix(K,1,1,[z^((i-1)*(j-1))]) :
         zero_matrix(K,0,0) for i in 1:3,j in 1:3,k in 1:3]
    set_braiding!(C,R)
    C,z
end

function audit_pointed_c2()
    N = zeros(Int,2,2,2)
    N[1,1,1] = N[1,2,2] = N[2,1,2] = N[2,2,1] = 1
    C = six_j_category(QQ,N,["1","g"])
    set_one!(C,1)
    set_braiding!(C,[identity_matrix(QQ,N[i,j,k])
                     for i in 1:2,j in 1:2,k in 1:2])
    C
end

# A skeletal morphism is a block matrix in the actual Hom space. In
# particular, Schur's lemma gives Hom(1,g)=0 in split Vec_C2 (EGNO Lemma
# 1.5.2); malformed blocks must not manufacture an isomorphism between them.
@testset "Skeletal morphisms respect their Hom spaces" begin
    C = audit_pointed_c2()
    @test int_dim(Hom(C[1],C[2])) == 0
    @test_throws ArgumentError morphism(C[1],C[2],
        [identity_matrix(QQ,1),zero_matrix(QQ,0,0)])
    @test_throws ArgumentError morphism(C[1],C[1],[identity_matrix(QQ,1)])
    L = GF(5)
    @test_throws ArgumentError morphism(C[1],C[1],
        [identity_matrix(L,1),zero_matrix(L,0,0)])

    # Matrices act on row vectors, so Hom(1⊕1,1) has a 2-by-1 unit block.
    f = morphism(C[1]⊕C[1],C[1],
        [matrix(QQ,2,1,[1,2]),zero_matrix(QQ,0,0)])
    @test id(C[1]) ∘ f == f
    @test int_dim(Hom(domain(f),codomain(f))) == 2
    @test iszero(morphism(C[1],C[2],
        [zero_matrix(QQ,1,0),zero_matrix(QQ,0,1)]))
end

# Morphisms between biproducts are matrices of component maps (EGNO (2015),
# Section 1.2). Skeletal Hom bases must use the same matrix-unit order as their
# coordinate vectors, including rectangular source and target multiplicities.
@testset "Skeletal Hom matrix-unit coordinates" begin
    C = six_j_category(GF(5),ones(Int,1,1,1),["1"])
    set_one!(C,[1])
    U = one(C)
    H = Hom(U ⊕ U,U ⊕ U ⊕ U)
    for f in basis(H)
        @test sum(express_in_basis(f,H) .* basis(H)) == f
    end
end

# EGNO Proposition 2.6.1 identifies tensor structures on the identity of a
# pointed category with group 2-cocycles modulo coboundaries.  Consequently a
# finite algebraic search must not present its output as a complete list unless
# completeness has independently been proved.
@testset "Verified monoidal candidates versus classification" begin
    C,_ = audit_pointed_c3()
    @test_throws ArgumentError monoidal_structures(id(C))

    # In rank one the normalized tensorator is forced by the unit constraint.
    V = VectorSpaces(QQ)
    structures = monoidal_structures(id(V))
    @test length(structures) == 1
    @test monoidal_functor_axiom(only(structures))
end


# EGNO Sections 2.10, 4.7, and 8.13 distinguish left and right pivotal
# traces and define T from twist eigenvalues. Rowell--Stong--Wang,
# arXiv:0712.1377v4, Section 5.3.3 supplies the pointed C3 data.
@testset "Pivotal and modular data" begin
    C,z = audit_pointed_c3()
    K = base_ring(C)
    @test pentagon_axiom(C) && hexagon_axiom(C)
    @test tmatrix(C) == diagonal_matrix(K.([1,z,z]))

    # A nontrivial character of C3 is pivotal but not spherical. Left and
    # right dimensions use inverse character values.
    set_pivotal!(C,[K(1),z,z^2])
    X = C[2]
    @test is_pivotal(C) # constant-time declaration installed by the setter
    @test is_pivotal(C;check=true) && !is_spherical(C;check=true)
    @test left_dim(X) == z && right_dim(X) == inv(z)
    @test (right_ev(X)⊗id(X)) ∘ inv_associator(X,right_dual(X),X) ∘
          (id(X)⊗right_coev(X)) == id(X)
    # The trivial C3 character is spherical (EGNO, Exercise 4.7.16). Once
    # declared, unrelated braiding installation and categorical relabeling
    # must preserve that structural flag.
    set_spherical!(C,K.([1,1,1]))
    @test is_spherical(C)
    set_braiding!(C,copy(C.braiding))
    @test is_pivotal(C) && is_spherical(C)
    TensorCategories.sort_simples!(C,[2,1,3])
    @test is_pivotal(C) && is_spherical(C)

    # Unit-containing pentagons are part of coherence. The skeletal API fixes
    # unit associators to identities and rejects conflicting supplied data.
    U = six_j_category(QQ,ones(Int,1,1,1))
    set_one!(U,1)
    @test_throws ArgumentError set_associator!(U,1,1,1,1,
                                                matrix(QQ,1,1,[2]);check=true)
    @test pentagon_axiom(U)
    @test TensorCategories.randomized_pentagon_axiom(U,1)
end

# The negligible ideal is the radical of the opposite-Hom trace pairing;
# see Etingof--Ostrik, arXiv:1801.04409v4, Definition 2.1. On a skeletal
# biproduct, EGNO (2015), Proposition 4.7.3 gives the dimension-weighted block
# trace used by the optimized implementation.
@testset "Trace pairings and weighted skeletal traces" begin
    F = GF(5)
    C = six_j_category(F,ones(Int,1,1,1),["1"])
    set_one!(C,[1])
    X = one(C) ⊕ one(C) ⊕ one(C)
    H = End(X)
    @test all(tr(f) == TensorCategories.left_trace(f) for f in basis(H))
    @test rank(trace_pairing(H,H)) == 9
    @test quotient_hom_dimension(X) == 9
    @test size(trace_pairing(zero(C),X)) == (0,0)
    @test_throws ArgumentError trace_pairing(Hom(one(C),X),Hom(one(C),X))

    N = zeros(Int,2,2,2)
    N[1,1,1] = N[1,2,2] = N[2,1,2] = N[2,2,1] = 1
    D = six_j_category(F,N,["1","g"])
    set_one!(D,[1,0])
    g = simples(D)[2]
    Y = g ⊕ g ⊕ g
    @test F(tr(id(Y))) == F(3)
    set_pivotal!(D,F.([1,-1]))
    @test tr(id(Y)) == TensorCategories.left_trace(id(Y))
    @test F(tr(id(Y))) == F(-3)
end

@testset "Enclosure equality versus numerical overlap" begin
    K = ArbField(64)
    C = six_j_category(K,ones(Int,1,1,1))
    set_one!(C,1)
    U = one(C)
    f,g,h = [morphism(U,U,[matrix(K,1,1,[K(s)])])
             for s in ["0 +/- 1","1.5 +/- 1","3 +/- 1"]]
    # Overlap is nontransitive numerical compatibility, not equality.
    @test overlaps(f,g) && overlaps(g,h) && !overlaps(f,h)
    @test f != g && g != h && f != h && f == f

    L = AcbField(32)
    D = six_j_category(L,ones(Int,1,1,1))
    set_one!(D,1)
    set_spherical!(D,[L(1)])
    set_braiding!(D,[identity_matrix(L,1) for _ in 1:1,_ in 1:1,_ in 1:1])
    # Nondegeneracy is tested at the precision of the coefficient field.
    @test is_modular(D)
    U = one(D)
    a = L("0.4142135624 +/- 2.72e-11")*L(0,1)
    h = morphism(U,U,[matrix(L,1,1,[a])])
    @test overlaps(L(h),a)
    A = U⊕U
    bad = morphism(A,A,[matrix(L,2,2,[a,L(1),L(0),a])])
    @test_throws ArgumentError L(bad)
end

# Maschke's theorem gives semisimplicity when the characteristic does not
# divide |G|; see EGNO (2015), Remark 4.2.14. Fusion additionally requires
# splitting; see Mäurer--Thiel, arXiv:2406.13438v2, Section 2.1.
@testset "Weak fusion and splitting over nonclosed fields" begin
    R = representation_category(GF(2), cyclic_group(3))
    # x^2+x+1 is irreducible over F2, so the two-dimensional simple has
    # endomorphism field F4: this category is weak fusion but not fusion.
    @test is_semisimple(R) && is_weak_fusion(R)
    @test is_braided(R)
    @test !is_split_semisimple(R) && !is_fusion(R)
    # Rigidity gives Hom(1,X⊗X*) ≅ End(X). For the two-dimensional simple X,
    # End(X)=F4 has F2-dimension two, so the generic split fallback must not
    # choose an arbitrary first basis vector. The representation backend's
    # explicit duality remains available.
    X = only(filter(X -> int_dim(X) == 2,simples(R)))
    @test coev(X) isa TensorCategories.GroupRepresentationMorphism
    @test ev(X) isa TensorCategories.GroupRepresentationMorphism
    @test_throws ArgumentError invoke(coev,Tuple{Object},X)
    @test_throws ArgumentError invoke(ev,Tuple{Object},X)
    # The unweighted sum of squared norms is likewise the split formula; a
    # weak fusion category needs the arbitrary-field categorical dimension.
    @test_throws ArgumentError dim(R)
    # F- and R-symbol arrays use scalar multiplicity spaces. They cannot encode
    # the F4 endomorphisms of X as if Schur's lemma gave End(X)=F2.
    @test_throws ArgumentError six_j_category(R)
    @test_throws ArgumentError TensorCategories.six_j_symbols(R)
    @test_throws ArgumentError TensorCategories.skeletal_braiding(R)
    @test is_fusion(representation_category(GF(5), cyclic_group(2)))
    M = representation_category(GF(3), cyclic_group(3))
    @test !is_weak_fusion(M) && !is_fusion(M)

    # EGNO Example 8.5.4 identifies Z(Vec_G) with Yetter--Drinfeld modules.
    # For G=C2 in characteristic two it contains Rep_F2(C2), which is not
    # semisimple by Maschke. Thus a split fusion input can have a center that
    # is neither multifusion nor modular when its global dimension vanishes.
    N = zeros(Int,2,2,2)
    N[1,1,1] = N[1,2,2] = N[2,1,2] = N[2,2,1] = 1
    C = six_j_category(GF(2),N)
    set_one!(C,1)
    Z = center(C)
    RZ = centralizer(C,simples(C))
    @test is_fusion(C) && iszero(dim(C))
    @test !is_semisimple(Z) && !is_weak_fusion(Z)
    @test !is_multifusion(Z) && !is_modular(Z)
    @test !is_semisimple(RZ) && !is_weak_fusion(RZ)
    @test !is_multifusion(RZ) && !is_fusion(RZ)
end

# Exact ordinary-characteristic representations exercise a different backend
# from the finite-field MeatAxe. See EGNO (2015), Section 2.3, Lemma 1.5.2,
# and Example 8.1.2.
@testset "Ordinary group representations over QQ" begin
    G = symmetric_group(3)
    C = representation_category(QQ,G)
    U = one(C)
    # The trivial representation exists over any field; constructing it must
    # not require a generating set for the infinite matrix group GL(1,QQ).
    @test int_dim(Hom(U,U)) == 1
    @test is_braided(C) && is_braided(representation_category(GF(3),G))
    @test Representation(C,g -> identity_matrix(QQ,1)) == U

    # Biproducts with zero have exactly two injections and projections in
    # either order, and their diagonal terms resolve the identity.
    for (A,B) in ((U,zero(C)),(zero(C),U))
        D,i,p = direct_sum(A,B)
        @test length(i) == length(p) == 2
        @test i[1]∘p[1] + i[2]∘p[2] == id(D)
    end

    # The standard two-dimensional representation is the sum-zero subspace
    # of the permutation representation on three letters.
    P = Representation(C,gens(G),[
        matrix(QQ,3,3,[Int(j == i^g) for i in 1:3,j in 1:3])
        for g in gens(G)])
    aug = morphism(P,U,matrix(QQ,3,1,[1,1,1]))
    X,inc = kernel(aug)
    @test int_dim(X) == 2 && is_simple(X) && int_dim(End(X)) == 1
    @test is_zero(aug ∘ inc) && int_dim(Hom(U,X)) == 0
    @test is_split_semisimple(C) && is_fusion(C)

    # Character theory gives std tensor std = 1 + sign + std, so Schur's
    # lemma gives a three-dimensional endomorphism algebra.
    @test int_dim(End(X⊗X)) == 3
    ok,f = is_isomorphic(U⊕X,X⊕U)
    @test ok && inv(f) ∘ f == id(U⊕X)

    # QQ[C3] has a nonsplit simple QQ(zeta_3), whose endomorphism field has
    # dimension two over QQ rather than one.
    H = cyclic_group(3)
    D = representation_category(QQ,H)
    Y = Representation(D,gens(H),[matrix(QQ,2,2,[0,1,-1,-1])])
    @test int_dim(End(Y)) == 2 && is_simple(Y)
    @test !is_split_semisimple(D) && !is_fusion(D)
    @test sort(int_dim.(simples(D;backend=:hecke))) == [1,2]
    @test sort([(int_dim(S),m) for (S,m) in
                composition_factors(regular_representation(D);backend=:hecke)]) ==
          [(1,1),(2,1)]
    L = splitting_field(D)
    @test degree(L) == 2
    # Enumeration needs a rational/Schur-index-aware backend; importing GAP's
    # absolutely irreducible list over a different field would be incorrect.
    @test_throws ArgumentError simples(C)
    @test_throws ArgumentError is_simple(X;backend=:gap)

    # The no-field constructor uses the abelian closure, where GAP's ordinary
    # irreducibles give all simples. This restores the established
    # characteristic-zero path without relying on an undocumented rational
    # irreducible-module interface.
    A = representation_category(G)
    @test rep(G) == A
    @test rep(QQ,G) == C
    @test int_dim.(simples(A)) == [1,1,2]
    reg = regular_representation(A)
    @test int_dim(reg) == 6
    @test sort([(int_dim(S),m) for (S,m) in composition_factors(reg)]) ==
          [(1,1),(1,1),(2,2)]
end

function audit_jordan_representation(C, n)
    J = identity_matrix(base_ring(C), n)
    for i in 1:n-1
        J[i, i+1] = 1
    end
    Representation(C, gens(base_group(C)), [J])
end

# A modular representation's matrices act on its underlying vector space.
# Categorical dimension is reduced in the base field and cannot determine the
# matrix size; compare EGNO (2015), Definition 4.7.11.
@testset "Modular representation matrix dimensions" begin
    F = GF(5)
    R = representation_category(F,cyclic_group(5))
    J3 = audit_jordan_representation(R,3)
    X = J3 ⊗ J3
    @test matrix(morphism(X,X,identity_matrix(F,9))) ==
          identity_matrix(F,9)
    @test_throws ErrorException morphism(X,X,identity_matrix(F,4);check=true)
    @test sum(int_dim(S)*m for (S,m) in decompose(X)) == 9
end

# Jordan--Hölder factors, socle types, and indecomposable summands are
# different notions; see EGNO (2015), §1.5, and GAP Reference Manual 69.5/69.7.
@testset "Representation simplicity and composition factors" begin
    R = representation_category(GF(5),cyclic_group(5))
    J1, J2 = audit_jordan_representation(R,1), audit_jordan_representation(R,2)
    @test is_simple(J1)
    # J2 is a nonsplit length-two extension with only one factor type.
    @test !is_simple(J2) && !is_simple(J1 ⊕ J1)
    @test !is_simple(J2;backend=:hecke)
    cf = composition_factors(J2)
    @test length(cf) == 1 && cf[1][2] == 2 && int_dim(cf[1][1]) == 1
    cfh = composition_factors(J2;backend=:hecke)
    @test length(cfh) == 1 && cfh[1][2] == 2 && int_dim(cfh[1][1]) == 1
    @test_throws ArgumentError decompose(J2;backend=:hecke)
    soc = simple_subobjects(J2)
    @test length(soc) == 1 && int_dim(only(soc)) == 1

    T = representation_category(GF(5),cyclic_group(1))
    @test only(simple_subobjects(one(T) ⊕ one(T))) isa GroupRepresentation

    R2 = representation_category(GF(2),cyclic_group(3))
    W = Representation(R2,gens(base_group(R2)),[matrix(GF(2),[0 1;1 1])])
    # x²+x+1 is irreducible: a simple need not be absolutely simple.
    @test is_simple(W) && int_dim(End(W)) == 2

    F3 = GF(3)
    G = matrix_group([matrix(F3,[1 1;0 1]),matrix(F3,[1 0;0 -1])])
    RG = representation_category(F3,G)
    Y = Representation(RG,gens(G),matrix.(gens(G)))
    @test length(composition_factors(Y)) == 2
    S = only(simple_subobjects(Y))
    # The sign socle embeds, while the trivial head is only a quotient.
    @test int_dim(Hom(S,Y)) == 1 && int_dim(Hom(Y,S)) == 0
end

@testset "Splitting fields and representation type" begin
    H = cyclic_group(3)
    C = representation_category(GF(2),H)
    L = splitting_field(C)
    @test degree(L) == 2
    @test !is_split_semisimple(C)
    @test is_split_semisimple(representation_category(L,H))
    @test is_finite_representation_type(C)
    @test sort(int_dim.(indecomposables(C))) == [1,2]

    M = representation_category(GF(5),cyclic_group(5))
    @test is_finite_representation_type(M)
    @test splitting_field(M) == GF(5)
    @test_throws ArgumentError indecomposables(M)

    I = representation_category(GF(2),symmetric_group(4))
    @test !is_finite_representation_type(I)
    @test_throws ArgumentError indecomposables(I)
end

# Finite-dimensional vector spaces are simple exactly in integer dimension
# one. Categorical dimension is reduced modulo the characteristic.
@testset "Vector-space simplicity in positive characteristic" begin
    V = VectorSpaces(GF(5))
    @test is_simple(VectorSpaceObject(V,1))
    @test !is_simple(VectorSpaceObject(V,6))
    @test dim(VectorSpaceObject(V,6)) == GF(5)(1)
end

# Deduplication must preserve categorical parents and coefficient fields; it
# mutates only its temporary list, never the contained morphisms.
@testset "Identity-preserving deduplication" begin
    C = graded_vector_spaces(GF(5),cyclic_group(2))
    f = id(C[one(base_group(C))])
    B = TensorCategories.unique_without_hash([f,f])
    @test length(B) == 1
    @test parent(only(B)) === C
    @test base_ring(matrix(only(B))) === base_ring(only(B))
end

# Schur's lemma has no converse in a nonsemisimple category. In the abelian
# arrow category, k --id→ k is a brick but contains the proper subobject 0→k.
@testset "Simplicity in the arrow category" begin
    V = VectorSpaces(GF(5))
    A = ArrowCategory(V)
    U = one(V)
    E = ArrowObject(A,id(U))
    S = ArrowObject(A,zero_morphism(zero(V),U))
    T = ArrowObject(A,zero_morphism(U,zero(V)))
    @test int_dim(End(E)) == 1
    @test !is_simple(E)
    @test is_simple(S) && is_simple(T)
    @test !is_simple(zero(A))
    @test object_type(A) == ArrowObject &&
          TensorCategories.morphism_type(A) == ArrowMorphism
end

# Semisimplification quotients Hom spaces by negligible morphisms; see
# Etingof--Ostrik, arXiv:1801.04409v4, Definition 2.1 and Proposition 2.4.
@testset "Semisimplification scalar coordinates" begin
    F = GF(5)
    R = representation_category(F, cyclic_group(5))
    J2 = audit_jordan_representation(R, 2)
    Q = Semisimplification(R)
    X = semisimplify(J2, Q)
    generator = matrix(J2(gens(base_group(R))[1]))
    f = semisimplify(morphism(J2, J2, generator), Q)

    # End(J2)=F5[N]/(N^2), and tr(Nh)=0 for every h. Thus [N]=0 and
    # the group generator [1+N] is the quotient identity.
    @test f == id(X)
    @test express_in_basis(f, [id(X)]) == [F(1)]
    @test F(f) == F(id(X)) == F(1)
    n = f - id(X)
    @test is_zero(n) && F(n) == F(0)
    @test F(F(3) * id(X) + n) == F(3)

    P = semisimplify(audit_jordan_representation(R, 5), Q)
    # Every endomorphism of J5 has trace divisible by five, so J5 becomes zero.
    @test F(id(P)) == F(0)

    U, inc, proj = direct_sum(X, X)
    # Projection onto one summand is not a scalar endomorphism of X⊕X.
    @test_throws ArgumentError F(inc[1] ∘ proj[1])
    # A nonzero map between different objects has no scalar value.
    @test_throws ArgumentError F(inc[1])
end

# A representation morphism is an intertwiner; see Etingof et al.,
# Introduction to Representation Theory (2011), Definition 1.13.
@testset "Representation morphism validation" begin
    F = GF(5)
    R = representation_category(F, cyclic_group(5))
    J2 = audit_jordan_representation(R, 2)
    bad = matrix(F, [1 0; 0 0])
    @test_throws ArgumentError morphism(J2, J2, bad; check = true)
    # Expensive equivariance validation is opt-in for performance.
    @test matrix(morphism(J2, J2, bad)) == bad

    S = representation_category(F, cyclic_group(2))
    T = Representation(S, gens(base_group(S)), [identity_matrix(F, 2)])
    @test_throws ArgumentError morphism(J2, T, identity_matrix(F, 2))
    @test_throws ArgumentError morphism(J2, J2, identity_matrix(GF(7), 2))
    @test_throws ErrorException morphism(J2, J2, identity_matrix(F, 1))
end

# Monomorphisms and epimorphisms need not split outside semisimple categories;
# see EGNO (2015), Sections 1.3--1.4.
@testset "Nonsplit monomorphisms and epimorphisms" begin
    R = representation_category(GF(5), cyclic_group(5))
    J1, J2 = audit_jordan_representation(R, 1), audit_jordan_representation(R, 2)
    p = only(basis(Hom(J2, J1)))
    i = only(basis(Hom(J1, J2)))
    @test is_epimorphism(p) && !is_monomorphism(p)
    @test is_monomorphism(i) && !is_epimorphism(i)
    # Splittings would decompose the indecomposable Jordan block J2.
    @test_throws ErrorException right_inverse(p)
    @test_throws ErrorException left_inverse(i)
    @test is_monomorphism(zero_morphism(zero(R), J1))
    @test is_epimorphism(zero_morphism(J1, zero(R)))
end
