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

function TambaraYamagami(A::GAPGroup, χ = nothing)
    n = Int(order(A))
    @assert is_abelian(A)
    
    m = Int(exponent(A))

    K,ξ = CyclotomicField(4*m, "ξ($(m^4))")
    _,x = K["x"]

    if is_square(n) 
        a = K(ZZ(√n))
    else
        a = roots(x^2-n)[1]
    end

    if χ === nothing
        χ = nondegenerate_bilinear_form(A, ξ^4)
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
    return TY
end


#=------------------------------------------------
    Calculation taken from 
    https://mathoverflow.net/questions/374021/is-there-a-non-degenerate-quadratic-form-on-every-finite-abelian-group
------------------------------------------------=#

function nondegenerate_bilinear_form(G::GAPGroup, ξ::FieldElem = CyclotomicField(exponent(G))[2])
    @assert is_abelian(G)

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

    
        
        
