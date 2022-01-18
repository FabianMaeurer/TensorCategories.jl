function groethendieck_ring(C::MultiTensorCategory, simples = simples(C))
    @assert issemisimple(C) "C is required to be semi-simple"

    m = multiplication_table(C,simples)

    Z = Integers{Int64}()

    A = AlgAss(Z, m, coefficients(one(C),simples))
    function to_gd(X)
        coeffs = coefficients(X,simples)
        z = AlgAssElem{Int64, AlgAss{Int64}}(A)
        z.coeffs = coeffs
        return z
    end
    return A,to_gd
end