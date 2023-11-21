#=----------------------------------------------------------
    Representation category of Uq(𝕤𝕝₂) 
----------------------------------------------------------=#

struct UqSl2Representations <: Category
    base_ring::Field
    q::RingElem
end

struct UqSl2rep <: Object
    parent::UqSl2Representations
    components::AbstractVector
end

struct UqSl2repMorphism <: Morphism
    domain::UqSl2rep
    codomain::UqSl2rep
    m::AbstractVector
end

function Morphism(X::UqSl2rep, Y::UqSl2rep, m::AbstractArray)
    UqSl2repMorphism(X,Y,m)
end

function Morphism(X::UqSl2rep, Y::UqSl2rep, m::Dict{Int, <:MatElem}) 
    UqSl2repMorphism(X,Y, sparsevec(m))
end

function Morphism(X::UqSl2rep, Y::UqSl2rep, ms::Pair{Int, <:MatElem}...) 
    UqSl2repMorphism(X,Y, Dict(ms...))
end

#=----------------------------------------------------------
    getter 
----------------------------------------------------------=#

matrices(f::UqSl2repMorphism) = f.m
matrix(f::UqSl2repMorphism) = diagonal_matrix(findnz(f.m)[2])

function getindex(C::UqSl2Representations, k::Int) 
    UqSl2rep(C,sparsevec([k+1],[1])) 
end

==(f::UqSl2repMorphism, g::UqSl2repMorphism) = 
    domain(f) == domain(g) &&
    codomain(f) == codomain(g) &&
    findnz(f.m) == findnz(g.m)
#=----------------------------------------------------------
    Functionality 
----------------------------------------------------------=#

getindex(X::UqSl2rep, k) = k ≤ length(X.components) ? X.components[k+1] : 0
getindex(f::UqSl2repMorphism, k) = k ≤ length(f.m) ? f.m[k+1] : zero_matrix(base_ring(f), domain(f)[k], codomain(f)[k]) 

function id(X::UqSl2rep)
    ind, vals = findnz(X.components)
    mats = sparsevec(ind, [diagonal_matrix(base_ring(X)(1), v) for v ∈ vals])
    Morphism(X,X,mats)
end

function zero_morphism(X::UqSl2rep, Y::UqSl2rep)
    if X == Y == zero(parent(X)) 
        return Morphism(X,Y, sparsevec(MatElem[]))
    end
    return Morphism(X,Y, Dict(k => zero_matrix(MatElem, base_ring(X), X[k-1], Y[k-1]) for k ∈ findnz(X.components)[1] ∪ findnz(Y.components)[1]))
end

zero(C::UqSl2Representations) = UqSl2rep(C,sparsevec([]))
one(C::UqSl2Representations) = UqSl2rep(C, sparsevec([1]))

function direct_sum(X::UqSl2rep, Y::UqSl2rep)
    X_length = length(X.components)
    Y_length = length(Y.components)

    if X_length > Y_length 
        S_components = deepcopy(X.components) 
        S_components[1:Y_length] .+= Y.components
    else
        S_components = deepcopy(Y.components)
        S_components[1:X_length] .+= X.components
    end

    S = UqSl2rep(parent(X), S_components)
    ix_mats = matrices(zero_morphism(X,S))
    iy_mats = matrices(zero_morphism(Y,S))
    px_mats = matrices(zero_morphism(S,X))
    py_mats = matrices(zero_morphism(S,Y))

    for i ∈ findnz(S.components)[1]
        x = i ≤ length(X.components) ? X[i] : 0
        y = i ≤ length(Y.components) ? Y[i] : 0
        for j ∈ 1:x 
            ix_mats[i][j,j] = 1
            px_mats[i][j,j] = 1
        end
        for j ∈ 1:y 
            iy_mats[i][j,j+x] = 1
            py_mats[i][j+x,j] = 1
        end
    end

    ix = Morphism(X,S, ix_mats)
    px = Morphism(S,X, px_mats)
    iy = Morphism(Y,S, iy_mats)
    py = Morphism(S,Y, py_mats)

    return S,[ix,iy],[px,py]
end

function direct_sum(f::UqSl2repMorphism, g::UqSl2repMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    F = base_ring(dom)

    mats = Dict{Int,MatElem}(k => diagonal_matrix(f[k], g[k]) for k ∈ findnz(f.m)[1] ∩ findnz(g.m)[1])

    not_in_g = [k => f[k-1] for k ∈ findnz(f.m)[1] if k ∉ findnz(g.m)[1]]
    not_in_f = [k => g[k-1] for k ∈ findnz(g.m)[1] if k ∉ findnz(f.m)[1]]
    length(not_in_f) ≠ 0 ? push!(mats, not_in_f...) : nothing
    length(not_in_g) ≠ 0 ? push!(mats, not_in_g...) : nothing

    return Morphism(dom,cod, mats)
end

function *(λ, f::UqSl2repMorphism)
    ind, vals = findnz(f.m)
    Morphism(domain(f), codomain(f), sparsevec(ind, λ .* vals))
end

function +(f::UqSl2repMorphism, g::UqSl2repMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g)
    mats = Dict(i => f.m[i] + g.m[i] for i ∈ findnz(f.m)[1])
    return Morphism(domain(f), codomain(f), mats)
end

function compose(f::UqSl2repMorphism, g::UqSl2repMorphism)
    mats = Dict(k => f[k-1]*g[k-1] for k ∈ findnz(f.m)[1] ∪ findnz(g.m)[1])
    Morphism(domain(f), codomain(g), mats)
end

function is_simple(X::UqSl2rep)
    return sum(X.components) == 1
end

function decompose(X::UqSl2rep)
    C = parent(X)
    return [(UqSl2rep(C, sparsevec(Dict(i => 1))), v) for (i,v) ∈ zip(findnz(X.components)...)]
end

function inv(f::UqSl2repMorphism)
    mats = Dict(i => inv(m) for (i,m) ∈ zip(findnz(f.m)...))
    return Morphism(codomain(f), domain(f), sparsevec(mats))
end

#=----------------------------------------------------------
    tensor_product 
----------------------------------------------------------=#

function tensor_product(X::UqSl2rep, Y::UqSl2rep)
    N = maximum(keys(X.components)) + maximum(keys(Y.components)) 
    T_components = sparsevec(zeros(Int,N))
    for (i,v) ∈ zip(findnz(X.components)...), (j,w) ∈ zip(findnz(Y.components)...)
        for k ∈ clebsch_gordan_rule(i-1,j-1)
            T_components[k+1] += v*w
        end
    end
    return UqSl2rep(parent(X), T_components)
end

function clebsch_gordan_rule(m::Int, n::Int)
    collect(abs(m-n):2:m+n) 
end


function tensor_product(f::UqSl2repMorphism, g::UqSl2repMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)

    C = parent(dom)

    h = zero_morphism(zero(C), zero(C))

    for (i,m) ∈ zip(findnz(f.m)...), (j,n) ∈ zip(findnz(g.m)...)
        A = kronecker_product(m,n)
        d1,d2 = size(A)
        #if d1*d2 == 0 continue end
        for k ∈ clebsch_gordan_rule(i-1,j-1)
            m = zero_morphism(C[k]^(d1),C[k]^(d2)).m
            m[k+1] = A

            h = h ⊕ Morphism(C[k]^(d1),C[k]^(d2), m)
        end
    end
    return h
end

#=----------------------------------------------------------
    Associator
----------------------------------------------------------=#

function simples_associator(C::UqSl2Representations, i::Int, j::Int, k::Int)
    #println("$lx, $ly, $lz, $lw")

    K = base_ring(C)
    mats = []
    q = quantum(C.q + inv(C.q), 2*(i+j+k))

    for w ∈ 1:i+j+k
        li = intersect(clebsch_gordan_rule(i,j),clebsch_gordan_rule(w,k))
        lj = intersect(clebsch_gordan_rule(i,w),clebsch_gordan_rule(j,k))
        gr = length(li)

        if gr == 0 continue end

        li =  

        push!(mats, w => matrix(K,gr,gr,[tl_six_j_symbol(q,j,i,w,k,n,m) for m in li, n in lj]))
    end
    dom = C[i] ⊗ C[j] ⊗ C[k]
    return Morphism(dom,dom,sparsevec(Dict(mats...)))
end

"""
    associator(X::UqSl2rep, Y::UqSl2rep, Z::UqSl2rep)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
 function associator(X::UqSl2rep, Y::UqSl2rep, Z::UqSl2rep)
    @assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"

    C = parent(X)

    dom = X⊗Y⊗Z

    if zero(C) == dom
        return zero_morphism(zero(C),zero(C))
    end

    F = base_ring(C)
    n = X.components[end] + Y.components[end] + Z.components[end]


    #---------------------------------
    # associators on simple objects
    #---------------------------------
    if is_simple(X) && is_simple(Y) && is_simple(Z)
        i = findnz(X.components)[1][1]
        j = findnz(Y.components)[1][1]
        k = findnz(Z.components)[1][1]
        return simples_associator(C,i-1,j-1,k-1)
    end

    #---------------------------------
    # associators for arbitrary objects
    #---------------------------------

    X_summands = vcat([[s for l ∈ 1:k] for (s,k) ∈ decompose(X)]...)
    Y_summands = vcat([[s for l ∈ 1:k] for (s,k) ∈ decompose(Y)]...)
    Z_summands = vcat([[s for l ∈ 1:k] for (s,k) ∈ decompose(Z)]...)

    #=-------------------------------------------------
        Distribution 
    -------------------------------------------------=#

    # Before
    distr_before = distribute_left(X_summands, Y) ⊗ id(Z)
    distr_before = (direct_sum([distribute_right(Xᵢ,Y_summands) for Xᵢ ∈ X_summands]...)⊗id(Z)) ∘ distr_before
    distr_before = distribute_left([Xᵢ⊗Yⱼ for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:], Z) ∘ distr_before
    distr_before = direct_sum([distribute_right(Xᵢ⊗Yⱼ,Z_summands) for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:]...) ∘ distr_before
    
    # After
    distr_after = id(X)⊗distribute_left(Y_summands, Z)
    distr_after = (id(X)⊗direct_sum([distribute_right(Yⱼ,Z_summands) for Yⱼ ∈ Y_summands]...)) ∘ distr_after
    distr_after = distribute_left(X_summands, Y⊗Z) ∘ distr_after
    YZ_arr = [Yⱼ⊗Zₖ for  Zₖ ∈ Z_summands, Yⱼ ∈ Y_summands][:]
    distr_after = direct_sum([distribute_right(Xᵢ, YZ_arr) for Xᵢ ∈ X_summands]) ∘ distr_after

    #-----------------------------------
    # Associator morphism
    #-----------------------------------
    m = zero_morphism(zero(C),zero(C))
    for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands
        m = m ⊕ associator(x,y,z)
    end

    return inv(distr_after) ∘ m ∘ distr_before
end

#=----------------------------------------------------------
    Hom Spaces  
----------------------------------------------------------=#

function Hom(X::UqSl2rep, Y::UqSl2rep)
    @assert parent(X) == parent(Y) "Mismatching parents"

    Xi, Yi = collect(X.components), collect(Y.components)
    F = base_ring(X)

    d = sum([x*y for (x,y) ∈ zip(Xi,Yi)])

    if d == 0 return HomSpace(X,Y,UqSl2repMorphism[], VectorSpaces(F)) end

    basis = [zero_morphism(X,Y).m for i ∈ 1:d]
    next = 1
    for k ∈ 1:minimum([length(Xi), length(Yi)])

        for i ∈ 1:Xi[k], j ∈ 1:Yi[k]
            basis[next][k][i,j] = 1
            next = next + 1
        end
    end
    basis_mors = [Morphism(X,Y,m) for m ∈ basis]
    return HomSpace(X,Y,basis_mors, VectorSpaces(F))
end
#=----------------------------------------------------------
    pretty printing 
----------------------------------------------------------=#

function show(io::IO, C::UqSl2Representations)
    print(io, """Representation category of Uq(𝔰𝔩₂) at q = $(C.q)""")
end

function show(io::IO, X::UqSl2rep)
    indices, values = findnz(X.components)
    str = ""
    if length(indices) == 0 
        print(io, "0")
        return
    end
    i,v = popfirst!(indices), popfirst!(values)

    if v == 1
        str *=  "V$(i-1)"
    else
        str *= "$v⋅V$(i-1)"
    end

    for (i,v) ∈ zip(indices, values)
        str *= " ⊕ "
        if v == 1
            str *=  "V$(i-1)"
        else
            str *= "$v⋅V$(i-1)"
        end
    end
    
    print(io, str)
end

function show(io::IO, f::UqSl2repMorphism)
    print(io, """Morphism between representations of Uq(𝔰𝔩₂) with 
    Domain: $(domain(f))
    Codomain: $(codomain(f))""")
end
