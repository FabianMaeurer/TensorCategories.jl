
struct Cocycle{N}
    group::Group
    F::Field
    m::Union{Nothing,Dict{NTuple{N,G},T}} where {G<:GroupElem,T<:FieldElem}
end

"""
    Cocycle(G::Group, m::Dict{NTuple{N,G}, T})

Return a ```N```-cocylce of ```G```. By now the condition is not checked.
"""
function Cocycle(G::Group, m::Dict{NTuple{N,S},T}) where {S<:GroupElem,T<:FieldElem,N}
    return Cocycle(G,parent(collect(values(m))[1]),m)
end

function Cocycle(G::Group, N::Int, f::Function)
    Cocycle(G,Dict(x => f(x...) for x ∈ Base.product([G for i ∈ 1:N]...)))
end

trivial_3_cocycle(G,F) = Cocycle{3}(G,F,nothing)
trivial_cocycle(G,F,k) = Cocycle{k}(G,F,nothing)

(c::Cocycle{N})(x...) where N = c.m === nothing ? c.F(1) : c.m[x]

function cyclic_group_3cocycle(G::Group, F::Field, ξ::FieldElem)
    g = G[1]
    n = order(G)
    D = Dict((g^i,g^j,g^k) => ξ^(div(i*(j+k - rem(j+k,n)),n)) for i ∈ 0:n-1, j ∈ 0:n-1, k ∈ 0:n-1)
    return Cocycle(G,D)
end

function show(io::IO, c::Cocycle{N}) where N
    print(io, "$N-Cocycle of $(c.group)")
end
