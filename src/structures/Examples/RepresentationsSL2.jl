#=----------------------------------------------------------
    Representation category of Uq(𝕤𝕝₂) 
----------------------------------------------------------=#

struct UqSl2Representations <: category
    base_ring::Field
    q::RingElem
end

struct UqSl2rep <: Object
    parent::UqSl2Representations
    components::Dict{ZZRingElem, ZZRingElem}
end

struct UqSl2repMorphism <: Morphism
    domain::UqSl2rep
    codomain::UqSl2rep
    m::Dict{ZZRingElem, <:MatElem}
end

function Morphism(X::UqSl2rep, Y::UqSl2rep, m::Dict{T, <:MatElem}) where T  
    UqSl2repMorphism(X,Y, order(Dict(ZZ(v) => n for (v,n) ∈ m)))
end

function Morphism(X::UqSl2rep, Y::UqSl2rep, ms::Pair{T, <:MatElem}) where T
    UqSl2repMorphism(X,Y, Dict(ms...))
end

#=----------------------------------------------------------
    getter 
----------------------------------------------------------=#

matrices(f::UqSl2repMorphism) = collect(values(f.m))
matrix(f::UqSl2repMorphism) = diagonal_matrix(matrices(f))

#=----------------------------------------------------------
    Functionality 
----------------------------------------------------------=#

getindex(X::UqSl2rep, k) = X.components[k] 

function zero_morphism(X::UqSl2rep, Y::UqSl2rep)
    return Morphism(X,Y, Dict(k => zero_matrix(MatElem, base_ring(X), X[k], Y[k]) for k ∈ keys(X.components) ∩ keys(Y.components)))
end

function direct_sum(X::UqSl2rep, Y::UqSl2rep)
    S = UqSl2rep(parent(X), merge(+, X.components, Y.components))
    ix_mats = matrices(zero_morphism(X,S))
    iy_mats = matrices(zero_morphism(Y,S))
    px_mats = matrices(zero_morphism(S,X))
    py_mats = matrices(zero_morphism(S,Y))

    for i ∈ keys(S)
        x = i ∈ keys(X.components) ? X[i] : 0
        y = i ∈ keys(Y.components) ? Y[i] : 0
        for j ∈ 1:x 
            ix_mats[i][j,j] = 1
            px_mats[i][j,j] = 1
        end
        for j ∈ 1:y 
            iy_mats[i][j,j+x] = 1
            py_mats[i][j+x,j] = 1
        end
    end

    ix = Morphism(X,S, Dict(k => m for (k,m) ∈ zip(keys(X), ix_mats)))
    px = Morphism(S,X, Dict(k => m for (k,m) ∈ zip(keys(X), px_mats)))
    iy = Morphism(Y,S, Dict(k => m for (k,m) ∈ zip(keys(Y), iy_mats)))
    py = Morphism(S,Y, Dict(k => m for (k,m) ∈ zip(keys(Y), px_mats)))

    return S,[ix,iy],[px,py]
end

function direct_sum(f::UqSl2repMorphism, g::UqSl2repMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    F = base_ring(dom)
    m = zero_morphism(dom,cod).m
    for i ∈ 1:parent(dom).simples
        mf,nf = size(f.m[i])
        mg,ng = size(g.m[i])
        z1 = zero(MatrixSpace(F,mf,ng))
        z2 = zero(MatrixSpace(F,mg,nf))
        m[i] = [f.m[i] z1; z2 g.m[i]]
    end

    return Morphism(dom,cod, m)
end

