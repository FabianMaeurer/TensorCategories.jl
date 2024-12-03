#=----------------------------------------------------------
    Build the G-Crossed product 𝒞 ⋊ G of a Fusion category 
    with a G-action on 𝒞. 
----------------------------------------------------------=#

function gcrossed_product(C::SixJCategory, T::GTensorAction)
    S = simples(C)
    G = group(T)
    K = base_ring(C) 

    irreducibles = ["($s,$g)" for  g ∈ elements(permutation_group(G)), s ∈ simples_names(C)][:]

    elements_of_G = elements(G)
    CxG = six_j_category(K, irreducibles)
    
    m,n = length(S), length(elements_of_G)

    mult = zeros(Int,n*m,n*m,n*m)

    for i1 ∈ 1:m, j1 ∈ 1:n, i2 ∈ 1:m, j2 ∈ 1:n
        g,h = elements_of_G[[j1,j2]]
        X = S[i1] ⊗ (T(g)(S[i2]))
        Y = g * h
        
        Y_ind = findfirst(==(Y),elements_of_G)

        for k ∈ 1:m
            mult[(i1-1)*n + j1, (i2-1)*n + j2, (k-1)*n + Y_ind] = X.components[k] 
        end
    end

    ass = Array{MatElem,4}(undef, n*m, n*m, n*m, n*m)

    for i1 ∈ 1:m, j1 ∈ 1:n, i2 ∈ 1:m, j2 ∈ 1:n, i3 ∈ 1:m, j3 ∈ 1:n
        g1,g2,g3 = elements_of_G[[j1,j2,j3]]

        X = S[i1] ⊗ (T(g1)(S[i2])) ⊗ (T(g1*g2)(S[i3]))
        Y = g1 * g2 * g3

        for k ∈ 1:m, l ∈ 1:n
            ass[(i1-1)*n + j1, (i2-1)*n + j2, (i3-1)*n + j3, (k-1)*n + l] = C.ass[i1,i2,i3,k]
        end
    end

    set_tensor_product!(CxG, mult)
    set_associator!(CxG, ass)

    one_coeffs = zeros(Int,n*m)
    one_C = one(C).components
    
    for i ∈ 1:m
        one_coeffs[(i-1)*n + 1] = one_C[i]
    end

    set_one!(CxG, one_coeffs)

    try 
        spheric = [C.spherical[i] for j ∈ 1:n, i ∈ 1:m][:]
        set_spherical(CxG, spheric)
    catch
    end

    set_name!(CxG, "Crossed product of $(C.name) and $G")

    return CxG
end


function ⋊(C::SixJCategory, G)
    gcrossed_product(C,G)
end

function gcrossed_product(C::SixJCategory, G::GAPGroup)
    # Define a canonical G-action on C. Might be trivial

    action = gtensor_action(C,G)
   
    return gcrossed_product(C, action)
end

