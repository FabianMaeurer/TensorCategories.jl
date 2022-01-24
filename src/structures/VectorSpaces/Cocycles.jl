
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
