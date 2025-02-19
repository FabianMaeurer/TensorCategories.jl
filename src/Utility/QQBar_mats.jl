function *(m::MatElem{QQBarFieldElem}, n::MatElem{QQBarFieldElem})
    a,b = size(m)
    c,d = size(n)
    @assert b == c "Mismatching dimensions"
    if a*b*c*d == 0
        return zero_matrix(QQBarField(),a,d)
    end

    deg = lcm([degree.(collect(m))[:]; degree.(collect(n))[:]])
    h = lcm([height_bits.(collect(m))[:]; height_bits.(collect(n))[:]])
    prec = 16*deg*h
    CC = ComplexField(prec)
    m2 = change_base_ring(CC,m)
    n2 = change_base_ring(CC,n)
    ret = m2*n2
    matrix(QQBarField(), [guess(QQBarField(), a, deg) for a ∈ ret])
end