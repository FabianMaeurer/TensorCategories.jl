#=----------------------------------------------------------
    https://web.math.ucsb.edu/~zhenghwa/data/research/pub/exotic-08.pdf 
----------------------------------------------------------=#

function E6subfactor()
    #K,ξ = CyclotomicField(24,"ξ₂₄") 
    K = QQBar
    r3 = sqrt(K(3))
    r2 = sqrt(K(2))
    ξ = root_of_unity(K,24)
    i = sqrt(K(-1))
    E6 = RingCategory(K,["𝟙","x","y"])

    mult = Array{Int,3}(undef,3,3,3)
    mult[1,1,:] = [1,0,0]
    mult[1,2,:] = mult[2,1,:] = [0,1,0]
    mult[1,3,:] = mult[3,1,:] = [0,0,1]
    mult[2,2,:] = [1,2,1]
    mult[2,3,:] = mult[3,2,:] = [0,1,0]
    mult[3,3,:] = [1,0,0]
    
    set_tensor_product!(E6,mult)

    set_associator!(E6,2,3,2,[matrix(K,1,1,[1]), matrix(K,[1 0; 0 -1]), matrix(K,1,1,[-1])])
    set_associator!(E6,3,2,3,[zero_matrix(K,0,0), matrix(K,1,1,[-1]), zero_matrix(K,0,0)])
    set_associator!(E6,2,2,3,[matrix(K,1,1,[1]), matrix(K,[0 i; -i 0]), matrix(K,1,1,[1])])
    set_associator!(E6,3,2,2,[matrix(K,1,1,[1]), matrix(K,[0 1; 1 0]), matrix(K,1,1,[1])])
    F1 = inv(r2)*ξ^7*matrix(K,[1 1; i -i])
    Fy = inv(r2)*ξ^7*matrix(K,[i -i; 1 1])

    d = (-1 + r3)//2
    k = inv(r2*sqrt(d))
    Fx = matrix(K,[d d -d//2*ξ^2 -d//2*ξ^8 -d//2*ξ^8 -d//2*ξ^2;
           d -d -d//2*ξ^2 -d//2*ξ^8 d//2*ξ^8 d//2*ξ^2;
           1 1 -inv(K(2))*(ξ^2 - 1) inv(K(2))*ξ^10 inv(K(2))*(ξ^(-3) + i) inv(K(2))*ξ^4;
           1 1 inv(K(2))*ξ^4 inv(K(2))*(ξ^(-3) + i) inv(K(2))*ξ^10 -inv(K(2))*(ξ^2 - 1);
           1 -1 -inv(K(2))*(ξ^2 - 1) inv(K(2))*ξ^10 -inv(K(2))*(ξ^(-3) + i) -inv(K(2))*ξ^4;
           -1 1 inv(K(2))*ξ^4 -inv(K(2))*(ξ^(-3) + i) inv(K(2))*ξ^10 -inv(K(2))*(ξ^2 - 1)])
    a = inv(K(4))*(1-r3)

    Fx = matrix(K,[0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 inv(K(2))*i 0 inv(K(2)) 0;
                0 0 0 0 0 0;])

    # Fx = matrix(K, [-2*a -2*a a*ξ^2 a*ξ^8 a*ξ^8 a*ξ^2;
    #                 -2*a 2*a a*ξ^2 a*ξ^8 -a*ξ^8 -a*ξ^2;
    #                 1 1 -inv(K(2))*(ξ^2-1) inv(K(2))*ξ^10 inv(K(2))*(ξ^(-4)+i) inv(K(2))*ξ^4;
    #                 1 1 inv(K(2))*ξ^4 inv(K(2))*(ξ^(-4)+i) inv(K(2))*ξ^10 -inv(K(2))*(ξ^2-1);
    #                 1 -1 -inv(K(2))*(ξ^2-1) inv(K(2))*(ξ^10) -inv(K(2))*(ξ^(-4)+i) -inv(K(2))*ξ^4;
    #                 -1 1 -inv(K(2))*ξ^4 -inv(K(2))*(ξ^(-4)+i) inv(K(2))*ξ^10 -inv(K(2))*(ξ^2-1)])

    set_associator!(E6,2,2,2,[F1,Fx,Fy])
    set_associator!(E6, [transpose(m) for m in E6.ass])
    set_name!(E6, "E6 subfactor fusion category")
    set_one!(E6, [1,0,0])
    return E6
end