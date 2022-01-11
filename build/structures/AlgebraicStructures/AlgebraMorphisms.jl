struct AlgebraMorphism{T}
    domain::Algebra
    codomain::Algebra
    m::Dict{R,S} where {R<:AlgebraElem,S<:AlgebraElem}
end
