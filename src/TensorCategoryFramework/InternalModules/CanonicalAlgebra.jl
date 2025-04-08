#=----------------------------------------------------------
    Construct the canonical algebra Hom(1,1) in C ⊠ Crev 
----------------------------------------------------------=#

function canonical_algebra(C::SixJCategory)
    C_rev = reversed_monoidal_category(C)

    C_C_rev = C ⊠ C_rev 

    A, = direct_sum([_deligne_product(C_C_rev,C[i],dual(C_rev[i])) for i ∈ 1:C.simples])
    
    u = Hom(one(C_C_rev),A)[1] 

    m = zero_morphism(zero(C_C_rev),A)

    bases = [[basis(Hom(C[i]⊗C[j],C[k])) for k ∈ 1:C.simples] for i ∈ 1:C.simples, j ∈ 1:C.simples]
    dual_bases = [[basis(Hom(dual(C[j])⊗dual(C[i]), dual(C[k]))) for k ∈ 1:C.simples] for i ∈ 1:C.simples, j ∈ 1:C.simples]

    dual_bases = [[cannonical_algebra_dual_basis(Bi, dual_Bi, C[i],C[j],C[k]) for (Bi,dual_Bi,k) ∈ zip(bases[i,j],dual_bases[i,j],1:C.simples)] for i ∈ 1:C.simples, j ∈ 1:C.simples]

    for i ∈ 1:C.simples, j ∈ 1:C.simples 
        B = bases[i,j]
        dual_B = dual_bases[i,j]
        dom = _deligne_product(C_C_rev, C[i]⊗C[j], dual(C[i]⊗C[j]))

        ϕ = vertical_direct_sum([length(B[k]) == 0 ? 
            zero_morphism(dom, _deligne_product(C_C_rev, C[k], dual(C[k]))) :
            sum([_deligne_product(C_C_rev, fi,gi) for (fi,gi) ∈ zip(B[k],dual_B[k])]) for k ∈ 1:C.simples])
        m = horizontal_direct_sum(m, ϕ)
    end

    AlgebraObject(C_C_rev, A, m, u)
end 

function _deligne_product(C::SixJCategory, f::SixJMorphism, g::SixJMorphism)
    f_times_g = zero_morphism(C) 

    n = parent(f).simples 
    m = parent(g).simples
    mats = Array{MatElem,1}(undef, n*m)
    for i ∈ 1:n, j ∈ 1:m 
        mat = kronecker_product(matrix(f[i]), matrix(g[j]))
        mats[(i-1)*m + j] = mat
    end

    dom_dims = [s for (s,_) ∈ size.(mats)]
    cod_dims = [s for (_,s) ∈ size.(mats)]

    X = direct_sum([v^k for (v,k) ∈ zip(simples(C), dom_dims)])[1]
    Y = direct_sum([v^k for (v,k) ∈ zip(simples(C), cod_dims)])[1]
    return morphism(X, Y, mats)
end

function _deligne_product(C::SixJCategory, X::SixJObject, Y::SixJObject)
    domain(_deligne_product(C,id(X), id(Y)))
end

function _canonical_algebra_pairing(f::Morphism, g::Morphism, X,Y,Z)
    K = base_ring(f) 
    
    K(compose(
        coev(X),
        (id(X) ⊗ coev(Y)) ⊗ id(dual(X)),
        inv_associator(X, Y, dual(Y)) ⊗ id(dual(X)),
        associator(X ⊗ Y, dual(Y), dual(X)),
        f ⊗ g, 
        spherical(Z) ⊗ id(dual(Z)),
        ev(dual(Z))
    ))
end

function cannonical_algebra_dual_basis(bases::Vector{<:Morphism}, dual_bases::Vector{<:Morphism}, X,Y,Z)
    if length(bases) == 0 
        return dual_bases 
    end

    if length(bases) > 1 
        error("Not yet implemented")
    end

    f, = bases 
    g, = dual_bases
    [inv(_canonical_algebra_pairing(f,g,X,Y,Z)) * g]
end
