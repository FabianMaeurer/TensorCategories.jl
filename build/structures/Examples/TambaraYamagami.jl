#=------------------------------------------------
    Construct arbitrary Tambara-Yamagami fusion
    categories.
------------------------------------------------=#

struct BilinearForm 
    group::GAPGroup
    base_ring::Field
    root_of_unity::FieldElem
    map::Function
end

function (B::BilinearForm)(x::GroupElem, y::GroupElem) 
    B.map(x,y,B.root_of_unity)
end

#=-------------------------------------------------
    ToDO: Centers of graded fusion categories (2009). Gelaki, Naidu, Nikhshych
            https://msp.org/ant/2009/3-8/ant-v3-n8-p05-s.pdf
-------------------------------------------------=#
function TambaraYamagami(A::GAPGroup, χ = nothing)
    n = Int(order(A))
    @assert is_abelian(A)
    
    m = Int(exponent(A))

    #K,ξ = CyclotomicField(8*m, "ξ($(8*m))")
    K = QQBar

    a = sqrt(K(n))

    if χ === nothing
        χ = nondegenerate_bilinear_form(A, root_of_unity(K,m))
    end

    els = elements(A)

    mult = Array{Int,3}(undef,n+1,n+1,n+1)
    for i ∈ 1:n, j ∈ 1:n
        k = findfirst(e -> e == els[i]*els[j], els)
        mult[i,j,:] = [l == k ? 1 : 0 for l ∈ 1:n+1]
        mult[j,i,:] = mult[i,j,:]
    end

    for i ∈ 1:n
        mult[i,n+1,:] = [[0 for i ∈ 1:n]; 1]
        mult[n+1,i,:] = [[0 for i ∈ 1:n]; 1]
    end

    mult[n+1,n+1,:] = [[1 for i ∈ 1:n]; 0]

    TY = RingCategory(K, mult, [["a$i" for i ∈ 1:n]; "m"])

    zero_mat = matrix(K,0,0,[])
    for i ∈ 1:n
        set_associator!(TY, n+1, i, n+1, [[χ(els[i],els[j])*matrix(id(TY[j])) for j ∈ 1:n]; zero_mat])
        for j ∈ 1:n
            set_associator!(TY, i, n+1, j, matrices(χ(els[i],els[j])*id(TY[n+1])))
        end
    end
    set_associator!(TY, n+1, n+1, n+1, [[zero_mat for _ ∈ 1:n]; inv(a)*matrix(K,[inv(χ(els[i],els[j])) for i ∈ 1:n, j ∈ 1:n])])
    set_one!(TY, [1; [0 for _ ∈ 1:n]])
    set_spherical!(TY, [K(1) for _ ∈ 1:n+1])

    set_name!(TY, "Tambara-Yamagami fusion category over $A")
    return TY
end


#=------------------------------------------------
    Calculation taken from 
    https://mathoverflow.net/questions/374021/is-there-a-non-degenerate-quadratic-form-on-every-finite-abelian-group
------------------------------------------------=#

function nondegenerate_bilinear_form(G::GAPGroup, ξ::FieldElem = CyclotomicField(Int(exponent(G)))[2])
    @assert is_abelian(G)

    if order(G) == 1
        return BilinearForm(G,parent(ξ),ξ,(x,y,_) -> one(parent(ξ)))
    end

    function nondeg(x::GroupElem, y::GroupElem, ξ::FieldElem)
        G = parent(x)
        m = GAP.gap_to_julia(GAP.Globals.AbelianInvariants(G.X))
        χ(m::Int) = isodd(m) ? 2 : 1

        x_exp = GAP.gap_to_julia(GAP.Globals.IndependentGeneratorExponents(G.X, x.X))
        y_exp = GAP.gap_to_julia(GAP.Globals.IndependentGeneratorExponents(G.X, y.X))

        return ξ^(ZZ(sum([mod(χ(mₖ)*aₖ*bₖ,mₖ) for (mₖ,aₖ,bₖ) ∈ zip(m,x_exp,y_exp)])))
    end

    return BilinearForm(G,parent(ξ),ξ,nondeg)
end

#-------------------------------------------------------------------------------
#   Examples
#-------------------------------------------------------------------------------

function Ising()
    #F,ξ = CyclotomicField(16, "ξ₁₆")
    F = QQBar
    a = sqrt(F(2))
    C = RingCategory(F,["𝟙", "χ", "X"])
    M = zeros(Int,3,3,3)

    M[1,1,:] = [1,0,0]
    M[1,2,:] = [0,1,0]
    M[1,3,:] = [0,0,1]
    M[2,1,:] = [0,1,0]
    M[2,2,:] = [1,0,0]
    M[2,3,:] = [0,0,1]
    M[3,1,:] = [0,0,1]
    M[3,2,:] = [0,0,1]
    M[3,3,:] = [1,1,0]

    set_tensor_product!(C,M)

    set_associator!(C,2,3,2, matrices(-id(C[3])))
    set_associator!(C,3,2,3, matrices((id(C[1]))⊕(-id(C[2]))))
    z = zero(MatrixSpace(F,0,0))
    set_associator!(C,3,3,3, [z, z, inv(a)*matrix(F,[1 1; 1 -1])])

    set_one!(C,[1,0,0])

    set_spherical!(C, [F(1) for s ∈ simples(C)])

    G = abelian_group(PcGroup, [2])
    χ = nondegenerate_bilinear_form(G,F(-1))

    # C = TambaraYamagami(G)

    # set_simples_name!(C,["𝟙","χ","X"])

    set_name!(C, "Ising fusion category")

    # set one of the four possible braidings 
    # http://arxiv.org/abs/2010.00847v1 (Ex. 4.13)
    ξ = root_of_unity(F,16)

    α = root_of_unity(F,8)

    braid = Array{MatElem,3}(undef, 3,3,3)
    a,b = elements(G)
    braid[1,1,:] = χ(a,a).*matrices(id(C[1]))
    braid[1,2,:] = braid[2,1,:] = χ(a,b).*matrices(id(C[2]))
    braid[2,2,:] = χ(b,b).*matrices(id(C[1]))

    braid[1,3,:] = braid[3,1,:] = matrices(id(C[3]))
    braid[2,3,:] = braid[3,2,:] = ξ^4 .* matrices(id(C[3]))

    braid[3,3,:] = α .* matrices((id(C[1]) ⊕ (inv(ξ^4) * id(C[2]))))

    set_braiding!(C,braid)
    return C
end

