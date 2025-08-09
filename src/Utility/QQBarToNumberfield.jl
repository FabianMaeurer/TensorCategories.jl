#=----------------------------------------------------------
    Transfer algebraic numbers into a number_field
    along a complex embedding       
----------------------------------------------------------=#

function preimage(e::AbsSimpleNumFieldEmbedding, x, deg = 80; tol = 10^(-10))

    # Get the number field in question
    K = number_field(e)

    # Get the Complex Field 
    CC = parent(e(K(1)))

    min = minpoly(x, deg) 

    # roots in domain 
    rs = roots(K, min)

    # Roots in Complex Field 
    complex_rs = e.(rs) 

    # find closest root 
    i = argmin(abs.(e.(rs) .- CC(x)))

    if abs(e(rs[i]) - CC(x)) > tol
        error("Preimage not in $K")
    end

    return rs[i]
end

function minpoly(x::AcbFieldElem, deg = 128)
    minpoly(guess(QQBarField(), x, deg))
end