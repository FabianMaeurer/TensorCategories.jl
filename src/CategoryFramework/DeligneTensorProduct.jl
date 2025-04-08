

function tensor_product(C::SixJCategory, D::SixJCategory, names1::Vector{String} = String[], names2::Vector{String} = String[])
    @assert is_multifusion(C)
    @assert is_multifusion(D)

    F = base_ring(C)
    try 
        F = parent(base_ring(C)(1)*base_ring(D)(1))
    catch
        throw(ErrorException("Cannot coerce to common base field"))
    end

    try 
        names1 = simples_names(C)
    catch
    end
    
    try 
        names2 = simples_names(D)
    catch
    end

    S = simples(C)
    T = simples(D)

    m,n = length.([S,T])

    mult = zeros(Int,n*m,n*m,n*m)

    for i1 ∈ 1:m, j1 ∈ 1:n, i2 ∈ 1:m, j2 ∈ 1:n
        X = S[i1]⊗S[i2]
        Y = T[j1]⊗T[j2]
        for k ∈ 1:m, l ∈ 1:n
            mult[(i1-1)*n + j1, (i2-1)*n + j2, (k-1)*n + l] = X.components[k] * Y.components[l] 
        end
    end

    ass = Array{MatElem,4}(undef, n*m, n*m, n*m, n*m)

    # for i1 ∈ 1:m, j1 ∈ 1:n, i2 ∈ 1:m, j2 ∈ 1:n, i3 ∈ 1:m, j3 ∈ 1:n
    #     X = S[i1]⊗S[i2]⊗S[i3]
    #     Y = T[j1]⊗T[j2]⊗T[j3]
    #     for k ∈ 1:m, l ∈ 1:n
    #         ass[(i1-1)*n + j1, (i2-1)*n + j2, (i3-1)*n + j3, (k-1)*n + l] = kronecker_product(change_base_ring(F, skel_C.ass[i1,i2,i3,k]), change_base_ring(F,skel_D.ass[j1,j2,j3,l]))
    #     end
    # end
 

    CD = six_j_category(F, mult, ["$s ⊠ $t" for t ∈ T, s ∈ S][:])

    N = n*m


    function _six_j_symbol(i,j,k,l) 
        i1,j1 = divrem(i-1,n) .+ [1,1]
        i2,j2 = divrem(j-1,n) .+ [1,1]
        i3,j3 = divrem(k-1,n) .+ [1,1]
        i4,j4 = divrem(l-1,n) .+ [1,1]
         kronecker_product(change_base_ring(F, six_j_symbol(C,i1,i2,i3,i4)),change_base_ring(F,six_j_symbol(D, j1,j2,j3,j4)))
    end

    set_tensor_product!(CD, mult)
    set_associator!(CD, Array{MatElem,4}(undef,N,N,N,N))

    set_attribute!(CD, :six_j_symbol, _six_j_symbol )

   # set_associator!(CD, ass)

    one_coeffs = zeros(Int,n*m)
    one_C = one(C).components
    one_D = one(D).components
    for i ∈ 1:m, j ∈ 1:n
        one_coeffs[(i-1)*n + j] = one_C[i]*one_D[j]
    end
    set_one!(CD, one_coeffs)

    try 
        spheric = [C.pivotal[i]*D.pivotal[j] for j ∈ 1:n, i ∈ 1:m][:]
        set_pivotal!(CD, spheric)
    catch e
        print(e)
    end

    return CD
end


⊠(C::Category, D::Category) = tensor_product(C,D)






# struct DeligneProduct <: Category
#     C::Category 
#     D::Category 
#     base_ring::Field
# end

# abstract type DeligneProdObject <: Object end

# struct SimpleDeligneProdObject <: DeligneProdObject
#     X::Object 
#     Y::Object 
#     parent::DeligneProduct 
# end

# struct CompositeDeligneProdObject <: DeligneProdObject
#     summands::Vector{SimpleDeligneProdObject}
#     parent::DeligneProduct
# end

# abstract type DeligneProdMorphism <: Morphism end

# struct SimpleDeligneProdMorphism <: DeligneProdMorphism 
#     domain::SimpleDeligneProdObject 
#     codomain::SimpleDeligneProdObject 
#     m::Tuple{Morphism,Morphism}
# end

# struct CompositeDeligneProdMorphism <: DeligneProdMorphism
#     domain::CompositeDeligneProdMorphism
#     codomain::CompositeDeligneProdMorphism
#     m::Matrix{SimpleDeligneProdMorphism}
# end

# #=------------------------------------------------
#     Constructors
# ------------------------------------------------=#
# function DeligneProduct(C::Category, D::Category)
#     @assert base_ring(C) == base_ring(D) "Missmatching base fields"
#     @assert is_semisimple(C) "Categories have to be semisimple"
#     @assert is_semisimple(D) "Categories have to be semisimple"
#     DeligneProduct(C,D,base_ring(C))
# end    

# ⊠(C::Category, D::Category) = DeligneProduct(C,D)

# function DeligneProdMorphism(f::Morphism, g::Morphism)
#     dom = domain(f) ⊠ domain(g)
#     cod = codomain(f) ⊠ codomain(g)
#     SimpleDeligneProdMorphism(dom, cod, (f,g))
# end

# ⊠(f::Morphism, g::Morphism) = DeligneProdMorphism(f,g)

# function DeligneProdObject(X::Object, Y::Object)
#     SimpleDeligneProdObject(X,Y,parent(X) ⊠ parent(Y))
# end

# ⊠(X::Object, Y::Object) = DeligneProdObject(X,Y)

# #=------------------------------------------------
#     simple objects
# ------------------------------------------------=#

# function simples(C::DeligneProduct) 
#     C_simples = simples(C.C)
#     D_simples = simples(C.D)
#     [DeligneProdObject(x,y) for x ∈ C_simples, y ∈ D_simples][:]
# end

# #=------------------------------------------------
#     Implied Structures
# ------------------------------------------------=#

# one(C::DeligneProduct) = one(C.C) ⊠ one(C.D)

# zero(C::DeligneProduct) = zero(C.C) ⊠ zero(C.D)



# #=------------------------------------------------
#     Checks
# ------------------------------------------------=#

# is_semisimple(C::DeligneProduct) = is_semisimple(C.C) && is_semisimple(C.D)
# is_abelian(C::DeligneProduct) = is_abelian(C.C) && is_abelian(C.D)
# is_linear(C::DeligneProduct) = is_linear(C.C) && is_linear(C.D)
# is_monoidal(C::DeligneProduct) = is_monoidal(C.C) && is_monoidal(C.D)
# isfinite(C::DeligneProduct) = isfinte(C.C) && isfinite(C.D)
# is_tensor(C::DeligneProduct) = is_tensor(C.C) && is_tensor(C.D)
# is_multitensor(C::DeligneProduct) = is_multitensor(C.C) && is_multitensor(C.D)
# is_fusion(C::DeligneProduct) = is_fusion(C.C) && is_fusion(C.D)
# is_multifusion(C::DeligneProduct) = is_multifusion(C.C) && is_multifusion(C.D)

# #=------------------------------------------------
#     Homspaces
# ------------------------------------------------=#

# function Hom(X::SimpleDeligneProdObject, Y::SimpleDeligneProdObject)
#     hom_X1_Y1 = basis(Hom(X.X,Y.X))
#     hom_X2_Y2 = basis(Hom(X.Y,Y.Y))
    
#     hom_basis = [f⊠g for f ∈ hom_X1_Y1, g ∈ hom_X2_Y2][:]

#     return HomSpace(X,Y,hom_basis, VectorSpaces(base_ring(X)))
# end