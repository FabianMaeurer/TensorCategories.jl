function Fibonacci(a::Int = 1)
    C = RingCategory(QQBar, ["𝟙", "τ"])

    _,x = QQ["x"]
    a = roots(x^2+x-1, QQBar)[a]

    M = zeros(Int, 2,2,2)

    M[1,1,:] = [1,0]
    M[1,2,:] = M[2,1,:] = [0,1]
    M[2,2,:] = [1,1]

    set_tensor_product!(C, M)

    set_associator!(C,2,2,2,2,matrix(QQBar, [a 1; -a -a]))
    set_name!(C, "Fibonacci fusion category")
    set_one!(C, [1,0])
    return C
end