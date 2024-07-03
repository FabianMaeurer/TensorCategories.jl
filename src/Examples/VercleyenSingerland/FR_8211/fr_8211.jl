#The Fusion ring and Associators were computed by Vercleyen and Singerland in https://arxiv.org/abs/2205.15637

#One needs to include the wanted associator of the other files ???
#In total there are 97 different associators; So far only 10 can be used; For the other one also has to check pentagox; Maybe the import does not work for all
@doc raw""" 

    cat_fr_8122(n::Int64)

Categorification of fusion ring FR8211. `n` chooses one of 96 possibily 
equivalent sets of associators.
"""
function cat_fr_8122(n::Int) #n gives the number of associator
    K = QQBar
    
    #(1,7,8) are class of one in D3/Z3; (2,3,4) are other class (5) is element t_1, (6) is other one
    C = six_j_category(K,["e", "a", "b", "aba", "t", "s", "ba", "ab"])


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

    ζ = root_of_unity(QQBar, 24)

    function modifier(x::Expr)
        replace!(x, :ζ, ζ)
        return x
    end

    dic = include(modifier, joinpath(@__DIR__, "asso_$n.jl"))


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

    set_canonical_spherical!(C)

    return C
end

