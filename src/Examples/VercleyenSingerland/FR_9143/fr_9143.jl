
function cat_fr_9143() 
    K = QQBar
    r2 = sqrt(K(2))
    z6 = root_of_unity(K,6) #6th root of unity

    function modifier(x::Expr)
        replace!(x, :r2, r2)
        replace!(x, :z6, z6)
        replace!(x, :K, K)
        return x
    end

    dic = include(modifier, joinpath(@__DIR__, "fr_9143_associator.jl"))#This imports a dictionary of associators.
    C = SixJCategory(K,["g0", "g3", "t2", "t1", "t0", "g4" , "g2", "g5", "g1"])
    #g0,g3,t2,t1,t0,g4,g2,g5,g1=simples(C)
    M = zeros(Int64,9,9,9)
    #Objects gi generate the sub Z/6Z
    # We have gi.tj=t(i+j) modulo ;
    # tj.gi=t(j-i) modulo 3
    # ti.tj=(i-j).(g0+g3)
    A=(9,7,2,6,8,1)
    B=(5,4,3)
    for i in 1:6
        for j in 1:6
            M[A[i],A[j],A[rem(i+j-1,6)+1]] = 1 
        end
    end
    for i in 1:6
        for j in 1:3
            M[A[i],B[j],B[rem(i+j-1,3)+1]]=1
            M[B[j],A[i],B[rem(j-i+5,3)+1]]=1
        end
    end
    for i in 1:3
        for j in 1:3
            M[B[i],B[j],A[rem(i-j+2,6)+1]]=1
            M[B[i],B[j],A[rem(i-j+5,6)+1]]=1
        end
    end
    set_tensor_product!(C,M)
    #Now we import the associators given. 
    kk=collect(keys(dic))
    for i in 1:9
        for j in 1:9
            for k in 1:9
                for l in 1:9
                    L=findall(x -> x[1]==i && x[2]==j && x[3]==k && x[4]==l, kk)
                    if length(L)==1 #In the case the associator is a 1x1 matrix
                        #println((i,j,k,l))
                        #print(L)
                        C.ass[i,j,k,l]=matrix(K,1,1,[dic[kk[L[1]]]])
                    elseif length(L)>0  #Here we have a 2x2 matrix and need to translate the dictionary into our matrices; It is not clear if this ordering works for all examples of associators
                        LL=kk[L]
                        #println(L)
                        #println("$i, $j, $k, $l")
                        LL=sort(LL, by = x -> x[6])
                        #print(L)
                        LL=sort(LL, by = x -> x[5])
                        #print(L)
                        M=matrix(K,2,2,[dic[l] for l in LL])
                        C.ass[i,j,k,l]=M
                    end
                end
            end
        end
    end    
    set_one!(C,[1,0,0,0,0,0,0,0,0])
    TensorCategories.set_spherical!(C, [K(1) for s ∈ simples(C)])
    return C
end