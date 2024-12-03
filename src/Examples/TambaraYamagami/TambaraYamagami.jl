
#=-------------------------------------------------
    ToDO: Centers of graded fusion categories (2009). Gelaki, Naidu, Nikhshych
            https://msp.org/ant/2009/3-8/ant-v3-n8-p05-s.pdf
-------------------------------------------------=#

""" 

    tambara_yamagami(A::GAPGroup)

Construct ``TY(A,τ,χ)`` over ℚ̅ where ``τ = √|A|`` and ``χ`` is a generic non-degenerate bilinear form.  
"""
function tambara_yamagami(A::GAPGroup) 
    m = Int(exponent(A))
    _, x = QQ[:x]
    K = splitting_field([x^m + 1, x^2 - order(A)])
    tambara_yamagami(K, A)
end

function tambara_yamagami(A::Int64...)
    tambara_yamagami(abelian_group(PcGroup, collect(A)))
end

function tambara_yamagami(K::Ring, A::Int64...)
    tambara_yamagami(K, abelian_group(PcGroup, collect(A)))
end


""" 

    tambara_yamagami(K::ring, A::GAPGroup)

Construct ``TY(A,τ,χ)`` over ``K`` where ``τ = √|A|`` and ``χ`` is a generic non-degenerate bilinear form.  
"""
function tambara_yamagami(K::Ring, A::GAPGroup) 
    # n = Int(order(A))     
    m = Int(exponent(A))
    sqrt_n = sqrt(K(Int(order(A))))
    χ = nondegenerate_bilinear_form(A, K)
    tambara_yamagami(K,A,sqrt_n,χ)
end

""" 

    tambara_yamagami(K::Ring, A::GAPGroup, τ::RingElem)

Construct ``TY(A,τ,χ)`` over ``K`` where ``χ`` is a generic non-degenerate bilinear form.  
"""
function tambara_yamagami(K::Ring, A::GAPGroup, sqrt_n::RingElem) 
    χ = nondegenerate_bilinear_form(A, K)
    tambara_yamagami(K,A,sqrt_n,χ)
end

""" 

    tambara_yamagami(K::Ring, A::GAPGroup, τ::RingElem)

Construct ``TY(A,τ,χ)`` over ``K`` where ``τ = √|A|``.  
"""
function tambara_yamagami(K::Ring, A::GAPGroup, χ::BilinearForm)
    # n = Int(order(A))     
    # m = Int(exponent(A))
    sqrt_n = sqrt(K(Int(order(A))))
    tambara_yamagami(K,A,sqrt_n,χ)
end

""" 

    tambara_yamagami(K::Ring, A::GAPGroup, τ::RingElem, χ::BilinearForm)

Construct the Category ``TY(A,τ,χ)``. 
"""
function tambara_yamagami(K::Ring, A::GAPGroup, τ::RingElem, χ::BilinearForm)
    n = Int(order(A))
    @assert is_abelian(A)
    
    m = Int(exponent(A))

    #K,ξ = cyclotomic_field(8*m, "ξ($(8*m))")
    
    try 
        a = sqrt(K(n))
    catch
        error("Base field needs to contain a square root of ord(A)")
    end
    a = τ 
    if χ === nothing
        χ = nondegenerate_bilinear_form(A, K)
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

    TY = six_j_category(K, mult, [["a$i" for i ∈ 1:n]; "m"])

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

    #Try to set a braiding. 
    #Ref: https://arxiv.org/pdf/2010.00847v1.pdf (Thm 4.9)
    # try 
    #     braid = Array{MatElem,3}(undef, n+1, n+1, n+1)
    #     q = x -> sqrt(χ(x,inv(x)))
        
    #     a = sqrt(inv(τ) * sum([q(x) for x ∈ A]))
    #     for i ∈ 1:n
    #         braid[i,n+1,:] = braid[n+1,i,:] = matrices(q(els[i]) * id(TY[n+1]))

    #         braid[n+1,n+1,1:n] = a .* [matrix(inv(q(els[k])) * id(TY[k])) for k ∈ 1:n]
    #         braid[n+1,n+1,n+1] = zero_matrix(K,0,0)
    #         for j ∈ i:n
    #             braid[i,j,:] = braid[j,i,:] = matrices(χ(els[i],els[j]) * id(TY[i]⊗TY[j]))
    #         end
    #     end
    #     set_braiding!(TY, braid)
    # catch e
    #     print(e)
    # end

    set_name!(TY, "Tambara-Yamagami fusion category over $A")
    return TY
end

@alias TY tambara_yamagami

#=------------------------------------------------
    Calculation taken from 
    https://mathoverflow.net/questions/374021/is-there-a-non-degenerate-quadratic-form-on-every-finite-abelian-group
------------------------------------------------=#

function nondegenerate_bilinear_form(G::GAPGroup, K::Field)
    @assert is_abelian(G)

    if order(G) == 1
        return BilinearForm(G,K,root_of_unity(K,2),Dict(x => one(K) for x ∈ G))
    end

    m = Int.(GAP.gap_to_julia(GAP.Globals.AbelianInvariants(G.X)))
    χ(m::Int) = isodd(m) ? 1 : 2
    
    x_exp = [GAP.gap_to_julia(GAP.Globals.IndependentGeneratorExponents(G.X, x.X)) for x ∈ G]
    
    M = Int(exponent(G))
    
    roots = [root_of_unity(K, χ(n)*n) for n ∈ m] 

    images = Dict(x => ((prod([r^aₖ^2 for (r,mₖ,aₖ) ∈ zip(roots,m,c)]))) for (x,c) ∈ zip(G, x_exp))
    

    return BilinearForm(G,K,root_of_unity(K,M),images)
end

function trivial_bilinear_form(G::GAPGroup, K::Field)
    @assert is_abelian(G)

    BilinearForm(G,K,root_of_unity(K,2), Dict(x => one(K) for x ∈ G))
end

#-------------------------------------------------------------------------------
#   Examples
#-------------------------------------------------------------------------------

""" 

    Ising()

Construct the Ising category over ``ℚ(√2)``.
"""
function Ising()
    _,x = QQ[:x]
    K,r2 = number_field(x^2-2, "√2")
    Ising(K, r2, 1)
end

""" 

    Ising(F::Ring)

Construct the Ising category over ``F``.
"""
function Ising(F::Ring)
    Ising(F, sqrt(F(2)))
end

""" 

    Ising(F::Ring, τ::RingElem)

Construct the Ising category with specific ``τ = √2``.
"""
function Ising(F::Ring, τ::RingElem)
    Ising(F,τ,1)
end

""" 

    Ising(F::Ring, q::Int)

Construct the braided Ising category over ``F`` where q = ±1 defined the braiding defined by ±i. 
"""
function Ising(F::Ring, q::Int)
    Ising(F,sqrt(F(2)),q)
end

""" 

    Ising(F::Ring, τ::RingElem, q::Int)

Construct the Ising fusion category where ``τ = √2`` a root and `q ∈ {1,-1}` specifies the braiding if it exists.
"""
function Ising(F::Ring, sqrt_2::RingElem, q::Int)
    #F,ξ = cyclotomic_field(16, "ξ₁₆")

    a = sqrt_2 
    C = six_j_category(F,["𝟙", "χ", "X"])
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
    z = zero(matrix_space(F,0,0))
    set_associator!(C,3,3,3, [z, z, inv(a)*matrix(F,[1 1; 1 -1])])

    set_one!(C,[1,0,0])

    set_spherical!(C, [F(1) for s ∈ simples(C)])

    
    # C = tambara_yamagami(G)

    # set_simples_name!(C,["𝟙","χ","X"])

    set_name!(C, "Ising fusion category")

    # set one of the four possible braidings 
    # http://arxiv.org/abs/2010.00847v1 (Ex. 4.13)
    
    try 

        G = abelian_group(PcGroup, [2])
        χ = nondegenerate_bilinear_form(G,F)

        ξ = q * root_of_unity(F,4)

        α = sqrt(inv(a)*(1 + ξ)) 

        braid = Array{MatElem,3}(undef, 3,3,3)
        a,b = elements(G)
        braid[1,1,:] = χ(a,a).*matrices(id(C[1]))
        braid[1,2,:] = braid[2,1,:] = χ(a,b).*matrices(id(C[2]))
        braid[2,2,:] = χ(b,b).*matrices(id(C[1]))

        braid[1,3,:] = braid[3,1,:] = matrices(id(C[3]))
        braid[2,3,:] = braid[3,2,:] = ξ .* matrices(id(C[3]))

        braid[3,3,:] = α .* matrices((id(C[1]) ⊕ (inv(ξ) * id(C[2]))))

        set_braiding!(C,braid)
    catch e
    end
    return C
end


function root_of_unity(K::Field, n::Int)
    Kx, x = K[:x]
    rs = roots(x^n - 1)
    divs = filter(e -> e != n, divisors(n))

    for r ∈ rs
        if K(1) ∈ r.^divs
            continue
        else
            return r
        end
    end
    error("There is no $n-th root of unity in the field")
end