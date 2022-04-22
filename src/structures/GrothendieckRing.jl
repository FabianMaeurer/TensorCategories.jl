function grothendieck_ring(C::Category, simples = simples(C))
    @assert isfusion(C) "C is required to be tensor"

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
