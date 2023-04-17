
#=----------------------------------------------------------
    Grothendieck Ring 
----------------------------------------------------------=#

"""
    grothendieck_ring(C::Category)

Return the grothendieck ring of the multiring category ``C``.
"""
function grothendieck_ring(C::Category, simples = simples(C))
    @assert is_multiring(C) "C is required to be tensor"

    m = multiplication_table(C,simples)

    A = ℕRing(ZZ.(m), ZZ.(coefficients(one(C),simples)))
    try 
        names = simples_names(C)
        A.basis_names = names
    catch
    end

    return A, Decategorification(C,A,simples)
end

# TODO: Needs Refinement
mutable struct Decategorification 
    domain::Category
    codomain::ℕRing
    simples::Vector{CategoryObject}
end

function (D::Decategorification)(X::CategoryObject)
    @assert D.domain = parent(X)
    coeffs = coefficients(X,simples)
    D.codomain(coeffs)
end

function show(io::IO, D::Decategorification)
    print(io, "Decategorification of $(D.domain)")
end


#=----------------------------------------------------------
    Generic Grothendieck Group 
----------------------------------------------------------=#

mutable struct GrothendieckGroup 
    category::Category
    is_ring::Bool
    base_ring::Ring
    objects::Vector{CategoryObject}
end

mutable struct GrothendieckGroupElem
    parent::GrothendieckGroup
    class::Vector{CategoryObject}
end

parent(x::GrothendieckGroupElem) = x.parent
equivalence_class(x::GrothendieckGroupElem) = x.class
representative(x::GrothendieckGroupElem) = equivalence_class(x)[1]

function set_class!(x::GrothendieckGroupElem, v::Vector{CategoryObject}) 
    x.class = v
end

function fuse_classes!(x::GrothendieckGroupElem, y::GrothendieckGroupElem)
    if !isempty(equivalence_class(x) ∩ equivalence_class(y)) 
        union_class = equivalence_class(x) ∪ equivalence_class(y)
        set_class!(x, union_class)
        set_class!(y, union_class)
        return x
    elseif is_isomorphic(representative(x), representative(y))[1]
        union_class = equivalence_class(x) ∪ equivalence_class(y)
        set_class!(x, union_class)
        set_class!(y, union_class)
        return x
    end
    error("Not the same element in the Grothendieck group")
end

function ==(x::GrothendieckGroupElem, y::GrothendieckGroupElem)
    @assert parent(x) == parent(y)
    if !isempty(equivalence_class(x) ∩ equivalence_class(y)) 
        union_class = equivalence_class(x) ∪ equivalence_class(y)
        set_class!(x, union_class)
        set_class!(y, union_class)
        return true
    elseif is_isomorphic(representative(x), representative(y))[1]
        union_class = equivalence_class(x) ∪ equivalence_class(y)
        set_class!(x, union_class)
        set_class!(y, union_class)
        return true
    end
    return false
end