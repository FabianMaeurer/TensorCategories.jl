""" 

    HaagerupH3([p1 = 1, p2 = 2])

Build the Haagerup ℋ₃ subfactor category. The category is build as SixJCategory. The associators are taken from the paper 

https://arxiv.org/pdf/1906.01322

where p1,p2 = ±1 are parameters for the different possible sets of associators.
"""
function HaagerupH3(K::Field = QQ; p1 = 1, p2 = 1)
    _,x = QQ["x"]
    K,_ = NumberField(x^16 - 3380*x^12 + 118368*x^10 + 814294*x^8 - 68093376*x^6 + 572623596*x^4 + 12977778528*x^2 + 1803785841)

    # _,x = K["x"]
    # r13 = roots(x^2-13)[2]

    r13 = sqrt(K(13))

    H = SixJCategory(K,["𝟙", "α", "α∗", "ρ", "αρ", "α∗ρ"])

    mult = Array{Int,3}(undef,6,6,6)

    mult[1,1,:] = [1,0,0,0,0,0]
    mult[2,1,:] = mult[1,2,:] = [0,1,0,0,0,0]
    mult[3,1,:] = mult[1,3,:] = [0,0,1,0,0,0]
    mult[4,1,:] = mult[1,4,:] = [0,0,0,1,0,0]
    mult[5,1,:] = mult[1,5,:] = [0,0,0,0,1,0]
    mult[6,1,:] = mult[1,6,:] = [0,0,0,0,0,1]

    mult[2,2,:] = [0,0,1,0,0,0]
    mult[2,3,:] = mult[3,2,:] = [1,0,0,0,0,0]
    mult[2,4,:] = mult[6,2,:] = [0,0,0,0,1,0]
    mult[2,5,:] = mult[4,2,:] = [0,0,0,0,0,1]
    mult[2,6,:] = mult[5,2,:] = [0,0,0,1,0,0]

    mult[3,3,:] = [0,1,0,0,0,0]
    mult[3,4,:] = mult[5,3,:] = [0,0,0,0,0,1]
    mult[3,5,:] = mult[6,3,:] = [0,0,0,1,0,0]
    mult[3,6,:] = mult[4,3,:] = [0,0,0,0,1,0]

    mult[4,4,:] = mult[5,5,:] = mult[6,6,:] = [1,0,0,1,1,1]
    mult[4,5,:] = mult[5,6,:] = mult[6,4,:] = [0,0,1,1,1,1]
    mult[4,6,:] = mult[6,5,:] = mult[5,4,:] = [0,1,0,1,1,1]

    set_tensor_product!(H,mult)

    set_one!(H,[1,0,0,0,0,0])

    #=----------------------------------------------------------
      load the associators from the file.
      Format:
        associator(H[i], H[j], H[k])[l] = [l,i,j,k,:,:]
    ----------------------------------------------------------=#    
    function modifier(x::Expr)
        replace!(x, :r13, r13)
        replace!(x, :p1, p1)
        replace!(x, :p2, p2)
        return x
    end

    dic = include(modifier, joinpath(@__DIR__, "Haagerup_associator.jl"))


    kk=collect(keys(dic))
    for i in 1:6
        for j in 1:6
            for k in 1:6
                for l in 1:6
                    L=findall(x -> x[2]==i && x[3]==j && x[4]==k && x[1]==l, kk)
                    if length(L)==1 #In the case the associator is a 1x1 matrix
                        #println((i,j,k,l))
                        #print(L)
                        H.ass[i,j,k,l]=matrix(K,1,1,[dic[kk[L[1]]]])
                    elseif length(L) > 0 #Here we have a 3x3 matrix and need to translate the dictionary into our matrices; It is not clear if this ordering works for all examples of associators
                        LL=kk[L]
                        #println(L)
                        LL=sort(LL, by = x -> x[5])
                        #print(L)
                        LL=sort(LL, by = x -> x[6])
                        #print(L)
                        n = Int(sqrt(length(L)))
                        M=matrix(K,n,n,[dic[l] for l in LL])
                        H.ass[i,j,k,l]=M
                    end
                end
            end
        end
    end

    set_name!(H, "Fusion category from Haagerup ℋ₃ subfactor")
    return H
end