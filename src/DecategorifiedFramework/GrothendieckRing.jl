
#=----------------------------------------------------------
    Grothendieck Ring 
----------------------------------------------------------=#

"""
    grothendieck_ring(C::Category)

Return the grothendieck ring of the multiring category ``C``.
"""
function split_grothendieck_ring(C::Category, simples = nothing)
    @assert is_multiring(C) "C is required to be multiring"

    m = if simples !== nothing 
         multiplication_table(C,simples)
    else
        multiplication_table(C)
    end

    one_coeffs = if simples !== nothing 
        coefficients(one(C),simples)
    else
        coefficients(one(C))
    end

    A = ℕRing(ZZ.(m), ZZ.(one_coeffs))
    try 
        names = simples_names(C)
        A.basis_names = names
    catch
    end

    

    if is_rigid(C)
        invol = [findfirst(j -> sum(m[i,j,Bool.(one_coeffs)]) > 0, 1:length(one_coeffs)) for i ∈ 1:length(one_coeffs)]
        set_attribute!(A, :involution, invol)
    end

    return A
end

# TODO: Needs Refinement
mutable struct Decategorification 
    domain::Category
    codomain::ℕRing
    simples::Vector{Object}
end

function (D::Decategorification)(X::Object)
    @assert D.domain == parent(X)
    coeffs = coefficients(X,D.simples)
    D.codomain(ZZ.(coeffs))
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
    objects::Vector{Object}
end

mutable struct GrothendieckGroupElem
    parent::GrothendieckGroup
    class::Vector{Object}
end

parent(x::GrothendieckGroupElem) = x.parent
equivalence_class(x::GrothendieckGroupElem) = x.class
representative(x::GrothendieckGroupElem) = equivalence_class(x)[1]

function set_class!(x::GrothendieckGroupElem, v::Vector{Object}) 
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