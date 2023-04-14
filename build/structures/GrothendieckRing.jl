
#=----------------------------------------------------------
    Grothendieck Ring 
----------------------------------------------------------=#

"""
    grothendieck_ring(C::Category)

Return the grothendieck ring of the multiring category ```C```.
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