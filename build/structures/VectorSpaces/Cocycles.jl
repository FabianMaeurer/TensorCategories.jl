
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

trivial_3_cocycle(G) = Cocycle{3}(G,nothing)

(c::Cocycle{N})(x...) where N = c.m == nothing ? 1 : c.m[x]
