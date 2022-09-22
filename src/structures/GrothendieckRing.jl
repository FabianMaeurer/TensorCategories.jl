
"""
    grothendieck_ring(C::Category)

Return the grothendieck ring of the multiring category ```C```.
"""
function grothendieck_ring(C::Category, simples = simples(C))
    @assert ismultiring(C) "C is required to be tensor"

    m = multiplication_table(C,simples)

    #Z = Integers{Int64}()
    Z = QQ

    A = AlgAss(Z, Z.(m), Z.(coefficients(one(C),simples)))
    function to_gd(X)
        coeffs = coefficients(X,simples)
        #z = AlgAssElem{Int64, AlgAss{Int64}}(A)
        z = A(coeffs)
        #z.coeffs = coeffs
        return z
    end
    return A,to_gd
end
