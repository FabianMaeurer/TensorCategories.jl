function E6subfactor()
    C24,ξ = CyclotomicField(24,"ξ₂₄")
    _,x = C24["x"]
    r3 = -ξ^6 + ξ^2
    r2 = -ξ^5 + ξ^3 + ξ
    K,v = NumberField(x^2 - r3 - 1, "√(1+√3)")
    d = 1 + r3
    i = ξ^6
    E6 = RingCategory(K,["𝟙","y","x"])

    mult = Array{Int,3}(undef,3,3,3)
    mult[1,1,:] = [1,0,0]
    mult[1,2,:] = mult[2,1,:] = [0,1,0]
    mult[1,3,:] = mult[3,1,:] = [0,0,1]
    mult[2,2,:] = [1,0,0]
    mult[2,3,:] = mult[3,2,:] = [0,0,1]
    mult[3,3,:] = [1,1,2]
    
    set_tensor_product!(E6,mult)

    set_associator!(E6,3,2,3,[matrix(K,1,1,[1]), matrix(K,1,1,[-1]), matrix(K,[1 0; 0 -1])])
    set_associator!(E6,2,3,2,[zero_matrix(K,0,0), zero_matrix(K,0,0), matrix(K,1,1,[-1])])
    set_associator!(E6,3,3,2,[matrix(K,1,1,[1]), matrix(K,1,1,[1]), matrix(K,[0 i; -i 0])])
    set_associator!(E6,2,3,3,[matrix(K,1,1,[1]), matrix(K,1,1,[1]), matrix(K,[0 1; 1 0])])
    F1 = inv(r2)*ξ^7*matrix(K,[1 i; 1 -i])
    Fy = inv(r2)*ξ^7*matrix(K,[i 1; -i 1])

    k = inv(r2*v)
    # Fx = matrix(K,[inv(d) inv(d)  k k k -k;
    #       inv(d) -inv(d) k k -k k;
    #       k*ξ^(-10) k*ξ^(-10) inv(r2*d)*ξ^(-5) 1//2*ξ^4 inv(r2*d)*ξ^(-10) 1//2*ξ^(-11);
    #       k*ξ^(-4) k*ξ^(-4) 1//2*ξ^10 inv(r2*d)*ξ 1//2*ξ^10 inv(r2*d)*ξ^(-11);
    #       k*ξ^(-4) k*ξ^(8) inv(r2*d)*ξ 1//2*ξ^10 inv(r2*d)*ξ^(-11) 1//2*ξ^10;
    #       k*ξ^(-10) k*ξ^2 1//2*ξ^4 inv(r2*d)*ξ^(-5) 1//2*ξ^(-8) inv(r2*d)*ξ^(-5)])
    a = 1//4*(1-r3)

    Fx = matrix(K, [-2*a -2*a a*ξ^2 a*ξ^8 a*ξ^8 a*ξ^2;
                    -2*a 2*a a*ξ^2 a*ξ^8 -a*ξ^8 -a*ξ^2;
                    1 1 -1//2*(ξ^2-1) 1//2*ξ^10 1//2*(ξ^(-4)+i) 1//2*ξ^4;
                    1 1 1//2*ξ^4 1//2*(ξ^(-4)+i) 1//2*ξ^10 -1//2*(ξ^4);
                    1 -1 -1//2*(ξ^2-1) 1//2*(ξ^(-4)+i) -1//2*(ξ^(-4)+i) -1//2*ξ^4;
                    -1 1 -1//2*ξ^4 -1//2*(ξ^(-4)+i) 1//2*ξ^10 -1//2*(ξ^2-1)])

    set_associator!(E6,3,3,3,[F1,Fy,Fx])
    #set_associator!(E6, [transpose(m) for m in E6.ass])
    set_name!(E6, "E6 subfactor fusion category")
    set_one!(E6, [1,0,0])
    return E6
end