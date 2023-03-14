struct SoergelDiagram
    domain::Vector{RingElem}
    codomain::Vector{RingElem}
    vertices::Int
    decorations::Vector{SoergelDecoration}
    edges::Matrix{Int}
    contact::Matrix{Int}
end


abstract type SoergelDecoration end
struct UnivalentVertex <: SoergelDecoration end
struct TrivalentVertex <: SoergelDecoration end
struct MultivalentVertex <: SoergelDecoration end
struct BoxVertex <: SoergelDecoration 
    polynomial::PolyElem
end
    

function frobenius_unit(D::SoergelDiagram, subgraph::Vector{})
end