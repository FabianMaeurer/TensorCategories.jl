
struct Cocycle{N}
    group::GAPGroup
    m::Union{Nothing,Dict{NTuple{N,G},T}} where {G<:GroupElem,T<:FieldElem}
end

"""
    Cocylce(G::GAPGroup, m::Dict{NTuple{N,G}, T})

Return a ```N```-cocylce of ```G```. By now the condition is not checked.
"""
function Cocycle(G::GAPGroup, m::Dict{NTuple{N,S},T}) where {S<:GroupElem,T<:FieldElem,N}
    return Cocycle{N}(G,m)
end

function Cocycle(G::GAPGroup, N::Int, f::Function)
    Cocycle(G,Dict(x => f(x...) for x ∈ Base.product([G for i ∈ 1:N]...)))
end

trivial_3_cocycle(G) = Cocycle{3}(G,nothing)

(c::Cocycle{N})(x...) where N = c.m == nothing ? 1 : c.m[x]

function cyclic_group_3cocycle(G::GAPGroup, F::Field, ξ::FieldElem)
    g = G[1]
    n = order(G)
    D = Dict((g^i,g^j,g^k) => ξ^(div(i*(j+k - rem(j+k,n)),n)) for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n)
    return Cocycle{3}(G,D)
end

function show(io::IO, c::Cocycle{N}) where N
    print(io, "$N-Cocycle of $(c.group)")
end
