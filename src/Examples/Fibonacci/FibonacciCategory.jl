function Fibonacci(K::Ring, a::Int = 1)

    C = six_j_category(K, ["𝟙", "τ"])

    _,x = K["x"]
    a = roots(x^2+x-1)[a]

    M = zeros(Int, 2,2,2)

    M[1,1,:] = [1,0]
    M[1,2,:] = M[2,1,:] = [0,1]
    M[2,2,:] = [1,1]

    set_tensor_product!(C, M)

    #set_associator!(C,2,2,2,1, matrix(K, 1,1, [a]))
    set_associator!(C,2,2,2,2, matrix(K, [a a; 1 -a]))
    set_name!(C, "Fibonacci fusion category")
    set_one!(C, [1,0])
    return C
end

function Fibonacci(a::Int = 1) 
    _,x = QQ[:x]
    K,ϕ = number_field(x^2 + x - 1, "ϕ")
    Fibonacci(K, a)
end