@doc raw""" 

    haagerup_H1()

Return the fusion category ℋ₁. 
"""
function haagerup_H1()
    load_fusion_category(joinpath(artifact"Haagerup_H1", "haagerup_H1"))
end

@doc raw""" 

    haagerup_H2()

Return the fusion category ℋ₂. 
"""
function haagerup_H2()
    C = anyonwiki(6,1,2,8,2,0,1)
    set_name!(C, "Fusion category ℋ₂")
    set_simples_names!(C, ["𝟙", "α", "α∗", "ρ", "αρ", "α∗ρ"])
    C
end

@doc raw""" 

    haagerup_H3()

Return the fusion category ℋ₃. 
"""
function haagerup_H3()
    C = anyonwiki(6,1,2,8,1,0,1)
    set_name!(C, "Fusion category ℋ₃")
    set_simples_names!(C, ["𝟙", "α", "α∗", "ρ", "αρ", "α∗ρ"])
    C
end

function haagerup_H3(K::QQBarField)
    C = haagerup_H3()
    e = C.embedding
    CC = AcbField()
    extension_of_scalars(C,K, embedding = x -> argmin(y -> abs(e(x) - CC(y)), roots(QQBarField(), minpoly(x))))
end

function haagerup_H3_center()
    path = joinpath(artifact"center_haagerup")
    C = load_fusion_category(joinpath(path, "center_haagerup"))
    CC = AcbField() 
    setfield!(C, :embedding, complex_embedding(base_ring(C), CC(1.29 + 0.25*im)))
    C
end

function numeric_unitary_center_H3()
    path = joinpath(artifact"NumericUnitaryCenterH3", "NumericUnitaryCenterH3")
    F_path = joinpath(path, "F_symbols.csv")
    R_path = joinpath(path, "R_symbols.csv")

    Z = load_numeric_fusion_category(F_path,R_path)
    set_simples_names!(Z, ["𝟙", "π₁", "π₂","σ₀", "σ₁", "σ₂", "μ₁", "μ₂", "μ₃", "μ₄", "μ₅", "μ₆"])
    Z
end


function unitary_haagerup_H3()
    _,x = QQ[:x]
    K,a = number_field(1 + x^2 - x^4 + x^6 + x^8)

    ass = include(joinpath(@__DIR__, "UnitaryH3.jl"))
    ass = Dict( k => K(v) for (k,v) ∈ ass)

    ass_mats = dict_to_associator(6, K, ass)

    mult = multiplication_table_from_F_symbols(ass_mats)

    H = six_j_category(K, mult,  ["𝟙", "α", "α∗", "ρ", "αρ", "α∗ρ"])

    set_associator!(H, transpose.(ass_mats))
    set_pivotal!(H, K.([1,1,1,1,1,1]))
    set_one!(H, [1,0,0,0,0,0])
    set_name!(H, "Fusion category ℋ₃")
    H.embedding = complex_embedding(K, AcbField()(-0.908677010511985310037597243 + 0.41749980906223339057410 * im))
    H
end
""" 

    unitary_haagerup_H3_wolf([p1 = 1, p2 = 2])

Build the Haagerup ℋ₃ subfactor category. The category is build as SixJCategory. The associators are taken from the paper 

https://arxiv.org/pdf/1906.01322

where p1,p2 = ±1 are parameters for the different possible sets of associators.
"""
function unitary_haagerup_H3_wolf(K::Ring = QQ; p1 = 1, p2 = 1)
    _,x = QQ["x"]
    if K == QQ
        K,_ = number_field(x^16 - 4*x^14 + 13*x^12 + 4*x^10 + 53*x^8 + 4*x^6 + 13*x^4 - 4*x^2 + 1)
    end
    #K,_ = number_field(x^16 - 3380*x^12 + 118368*x^10 + 814294*x^8 - 68093376*x^6 + 572623596*x^4 + 12977778528*x^2 + 1803785841)

    # _,x = K["x"]
    # r13 = roots(x^2-13)[2]

    r13 = sqrt(K(13))

    try
        sqrt(6*(1 + r13))
        sqrt((-3 + r13)//2)
    catch
        r13 = -r13
    end
    
    H = six_j_category(K,["𝟙", "α", "α∗", "ρ", "αρ", "α∗ρ"])

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

    dic = include(modifier, joinpath(@__DIR__, "Haagerup_H3_associator.jl"))


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

function haagerup_H3(s::Symbol; p1 = 1, p2 = 1)
    if s == :splitting_field
        _,x = QQ[:x]
        K,a = number_field()

        return haagerup_H3(K, p1 = p1, p2 = p2)
    end

    error("unknown keyword $s")
end