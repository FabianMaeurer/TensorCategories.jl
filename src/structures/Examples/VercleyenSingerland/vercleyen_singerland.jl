#The Fusion ring and Associators were computed by Vercleyen and Singerland in https://arxiv.org/abs/2205.15637

#One needs to include the wanted associator of the other files ???
#In total there are 97 different associators; So far only 10 can be used; For the other one also has to check pentagox; Maybe the import does not work for all
function cat_fr_8122(n::Int) #n gives the number of associator
    K,ζ = CyclotomicField(24)
    
    dic=include(joinpath(@__DIR__, "asso_$n.jl"))#This imports a dictionary of associators.

    #(1,7,8) are class of one in D3/Z3; (2,3,4) are other class (5) is element t_1, (6) is other one
    C = RingCategory(K,["e", "a", "b", "aba", "t", "s", "ba", "ab"])


    M = zeros(Int64,8,8,8)

    #Multiplication table of dihedral group
    #(Maybe make this shorter; I dont know how)

    M[1,1,:]=[1,0,0,0,0,0,0,0]
    M[1,2,:]=[0,1,0,0,0,0,0,0]
    M[1,3,:]=[0,0,1,0,0,0,0,0]
    M[1,4,:]=[0,0,0,1,0,0,0,0]
    M[1,7,:]=[0,0,0,0,0,0,1,0]
    M[1,8,:]=[0,0,0,0,0,0,0,1]
    M[2,1,:]=[0,1,0,0,0,0,0,0]
    M[2,2,:]=[1,0,0,0,0,0,0,0]
    M[2,3,:]=[0,0,0,0,0,0,0,1]
    M[2,4,:]=[0,0,0,0,0,0,1,0]
    M[2,7,:]=[0,0,0,1,0,0,0,0]
    M[2,8,:]=[0,0,1,0,0,0,0,0]
    M[3,1,:]=[0,0,1,0,0,0,0,0]
    M[3,2,:]=[0,0,0,0,0,0,1,0]
    M[3,3,:]=[1,0,0,0,0,0,0,0]
    M[3,4,:]=[0,0,0,0,0,0,0,1]
    M[3,7,:]=[0,1,0,0,0,0,0,0]
    M[3,8,:]=[0,0,0,1,0,0,0,0]
    M[4,1,:]=[0,0,0,1,0,0,0,0]
    M[4,2,:]=[0,0,0,0,0,0,0,1]
    M[4,3,:]=[0,0,0,0,0,0,1,0]
    M[4,4,:]=[1,0,0,0,0,0,0,0]
    M[4,7,:]=[0,0,1,0,0,0,0,0]
    M[4,8,:]=[0,1,0,0,0,0,0,0]
    M[7,1,:]=[0,0,0,0,0,0,1,0]
    M[7,2,:]=[0,0,1,0,0,0,0,0]
    M[7,3,:]=[0,0,0,1,0,0,0,0]
    M[7,4,:]=[0,1,0,0,0,0,0,0]
    M[7,7,:]=[0,0,0,0,0,0,0,1]
    M[7,8,:]=[1,0,0,0,0,0,0,0]
    M[8,1,:]=[0,0,0,0,0,0,0,1]
    M[8,2,:]=[0,0,0,1,0,0,0,0]
    M[8,3,:]=[0,1,0,0,0,0,0,0]
    M[8,4,:]=[0,0,1,0,0,0,0,0]
    M[8,7,:]=[1,0,0,0,0,0,0,0]
    M[8,8,:]=[0,0,0,0,0,0,1,0]

    #The classes permute the special particles; their product is sum of elemens of a class
    for i in (1,7,8)
        M[i,5,5]=1
        M[5,i,5]=1
        M[i,6,6]=1
        M[6,i,6]=1

        M[5,5,i]=1
        M[6,6,i]=1
    end

    for i in (2,3,4)
        M[i,5,6]=1
        M[5,i,6]=1
        M[i,6,5]=1
        M[6,i,5]=1

        M[5,6,i]=1
        M[6,5,i]=1
    end

    set_tensor_product!(C,M)

    #Now we import the associators given. 
    kk=collect(keys(dic))
    for i in 1:8
        for j in 1:8
            for k in 1:8
                for l in 1:8
                    L=findall(x -> x[1]==i && x[2]==j && x[3]==k && x[4]==l, kk)
                    if length(L)==1 #In the case the associator is a 1x1 matrix
                        #println((i,j,k,l))
                        #print(L)
                        C.ass[i,j,k,l]=matrix(K,1,1,[dic[kk[L[1]]]])
                    elseif length(L) > 0 #Here we have a 3x3 matrix and need to translate the dictionary into our matrices; It is not clear if this ordering works for all examples of associators
                        LL=kk[L]
                        #println(L)
                        LL=sort(LL, by = x -> x[6])
                        #print(L)
                        LL=sort(LL, by = x -> x[5])
                        #print(L)
                        M=matrix(K,3,3,[dic[l] for l in LL])
                        C.ass[i,j,k,l]=M
                    end
                end
            end
        end
    end

    set_one!(C,[1,0,0,0,0,0,0,0])

    TensorCategories.set_spherical!(C, [K(1) for s ∈ simples(C)])

    return C
end


function cat_fr_9143() #n gives the number of associator
    K = QQBar
    r2=sqrt(K(2))
    z6=(K(-1))^(K(1)//3) #6th root of unity

    dic=include(joinpath(@__DIR__, "asso_9143_1.jl"))#This imports a dictionary of associators.
    C = RingCategory(K,["g0", "g3", "t2", "t1", "t0", "g4" , "g2", "g5", "g1"])
    #g0,g3,t2,t1,t0,g4,g2,g5,g1=simples(C)
    M = zeros(Int64,9,9,9)
    #Objects gi generate the subgroup Z/6Z
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