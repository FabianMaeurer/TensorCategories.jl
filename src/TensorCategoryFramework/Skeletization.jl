function SixJCategory(C::Category, names::Vector{String} = ["X$i" for i ∈ 1:length(simples(C))])
    @assert is_multifusion(C)

    if typeof(C) == SixJCategory 
        return C
    end
    
    S = simples(C)
    n = length(S)
    F = base_ring(C)

    mult = Array{Int,3}(undef,n,n,n)

    for i ∈ 1:n, j ∈ 1:n
        mult[i,j,:] = [int_dim(Hom(S[k],S[i]⊗S[j])) for k ∈ 1:n]
    end

    # Define SixJCategory
    skel_C = SixJCategory(F,names)

    set_tensor_product!(skel_C,mult)

    set_associator!(skel_C, six_j_symbols(C, S))

    if multiplicity(C) > 1
        @warn "6j-Symbols might be wrong since multiplicity is greater than one"
    end
    # Try to set spherical
    try 
        set_spherical!(skel_C,[F(spherical(s)) for s ∈ S])
    catch
    end

    if is_braided(C)
        set_braiding!(skel_C, skeletal_braiding(C,S))
    end
    set_one!(skel_C, [int_dim(Hom(one(C), s)) for s ∈ S])

    set_name!(skel_C, "Skeletization of $C")

    return skel_C
end

function six_j_symbols(C::Category, S = simples(C))
    @assert is_multifusion(C)

    N = length(S)
    C_morphism_type = morphism_type(C)
    F = base_ring(C) 

    ass = Array{MatElem}(undef,N,N,N,N)
    homs = Dict()

    one_indices = findall(s -> int_dim(Hom(s,one(C))) > 0 , S)

    for (i,j,k,l) ∈ Base.product(1:N, 1:N, 1:N, 1:N)

        X,Y,Z,W = S[[i,j,k,l]]

        # Build a basis for Hom((X⊗Y)⊗Z,W)
        B_XY_Z_W = C_morphism_type[]
        for n ∈ 1:N
            V = S[n]

            H_XY_V = if (XYV = (X⊗Y,V)) ∈ keys(homs) 
                        homs[XYV] 
                    elseif !isempty([i,j] ∩ one_indices) && n ∈ [i,j] && (n ∉ one_indices || i == j)
                        push!(homs, XYV => [id(XYV[1])])[XYV]
                    else
                        push!(homs, XYV => basis(Hom(XYV[1], V)))[XYV]
                    end

            H_VZ_W = if (VZW = (V⊗Z,W)) ∈ keys(homs) 
                        homs[VZW] 
                    elseif !isempty([n,k] ∩ one_indices) && l ∈ [n,k] && (l ∉ one_indices || n == k)
                        push!(homs, VZW => [id(VZW[1])])[VZW]
                    else
                        push!(homs, VZW => basis(Hom(VZW[1], W)))[VZW]
                    end

            B = [f ∘ (g ⊗ id(Z)) for g ∈ H_XY_V, f ∈ H_VZ_W][:]
            if length(B) == 0 continue end
            B_XY_Z_W = [B_XY_Z_W; B]
        end

        # Build a basis for Hom(X⊗(Y⊗Z),W)
        B_X_YZ_W = C_morphism_type[]
        for n ∈ 1:N
            V = S[n]

            H_YZ_V = if (YZV = (Y⊗Z,V)) ∈ keys(homs)
                        homs[YZV]
                    elseif !isempty([j,k] ∩ one_indices) && n ∈ [j,k] && (n ∉ one_indices || j == k)
                        push!(homs, YZV => [id(YZV[1])])[YZV]
                    else
                        push!(homs, YZV => basis(Hom(YZV[1], V)))[YZV]
                    end

            H_XV_W = if (XVW = (X⊗V,W)) ∈ keys(homs)
                        homs[XVW]
                    elseif !isempty([i,n] ∩ one_indices) && l ∈ [i,n] && (l ∉ one_indices || i == n)
                        push!(homs, XVW => [id(XVW[1])])[XVW]
                    else
                        push!(homs, XVW => basis(Hom(XVW[1], W)))[XVW]
                    end
                    
            B = [f ∘ (id(X) ⊗ g) for g ∈ H_YZ_V, f ∈ H_XV_W][:]
            if length(B) == 0 continue end
            B_X_YZ_W = [B_X_YZ_W; B]
        end
        
        # Express the asociator in the corresponding basis
        a = associator(X,Y,Z)
        
        associator_XYZ_W = hcat([express_in_basis(f ∘ a, B_XY_Z_W) for f ∈ B_X_YZ_W]...)

        ass[i,j,k,l] = matrix(F, length(B_X_YZ_W), length(B_XY_Z_W), associator_XYZ_W)
            
    end

    return ass           
end

function skeletal_braiding(C::Category, S = simples(C))
    @assert is_braided(C)
    
    N = length(S)
    C_morphism_type = morphism_type(C)
    F = base_ring(C) 
    braid = Array{MatElem}(undef,N,N,N)

    for (i,j,l) ∈ Base.product(1:N,1:N,1:N)
        X,Y,W = S[[i,j,l]]
        # Basis for Hom(X⊗Y,W)
        B_XY_W = basis(Hom(X⊗Y,W))

        # Basis for Hom(Y⊗X,W)
        B_YX_W = basis(Hom(Y⊗X,W))

        braid_XY_W = hcat([express_in_basis(f ∘ braiding(X,Y), B_XY_W) for f ∈ B_YX_W]...)
        braid[i,j,l] = matrix(F, length(B_XY_W), length(B_YX_W), braid_XY_W)
    end
    return braid
end
# function six_j_symbols(C::Category, S = simples(C))
#     @assert is_multifusion(C)

#     n = length(S)
#     ass = Array{MatElem,4}(undef,n,n,n,n)

#     objects = []
#     isos = []

#     for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
#         X = (S[i] ⊗ S[j]) ⊗ S[k]
#         Y = S[i] ⊗ (S[j] ⊗ S[k])

#         _, ij_decomposition_morphism, ij_incl, ij_proj = direct_sum_decomposition(S[i] ⊗ S[j], S)
#         _, jk_decomposition_morphism, jk_incl, jk_proj =  direct_sum_decomposition(S[j] ⊗ S[k], S)

#         distribute_before = (inv(ij_decomposition_morphism) ⊗ id(S[k])) ∘ inv(distribute_left([domain(i) for i ∈ ij_incl], S[k]))

#         distribute_after = distribute_right(S[i], [domain(i) for i ∈ jk_incl]) ∘ (id(S[i]) ⊗ jk_decomposition_morphism)

#         # @show codomain(distribute_after) == domain(distribute_before)

#         summands = vcat([[x for _ ∈ 1:k] for (x,k) ∈ decompose(X, S)]...)

#         Z, before, incl, _ = direct_sum_decomposition(X,S)
#         Z2, after, _, proj = direct_sum_decomposition(Y,S)
        
#         before = inv(before) ∘ distribute_before ∘ before
#         after = after ∘ distribute_after ∘ inv(after)
#         # if X ∈ objects
#         #     before = isos[findfirst(e -> e == X, objects)]
#         # else
#         #     _,before = is_isomorphic(Z,X)
#         #     isos = [isos; before]
#         #     objects = [objects; X]
#         # end

#         # if Y ∈ objects
#         #     after = isos[findfirst(e -> e == Y, objects)]
#         # else
#         #     _,after = is_isomorphic(Z2,Y)
#         #     isos = [isos; after]
#         #     objects = [objects; Y]
#         # end

    
#         ass_mor =  after ∘ associator(S[i],S[j],S[k]) ∘ before 

#         for l ∈ 1:n
#             before_incl_l = filter(f -> domain(f) == S[l], incl)
#             after_proj_l = filter(f -> codomain(f) == S[l], proj)
            
#             if length(before_incl_l) == 0
#                 ass[i,j,k,l] = zero_matrix(base_ring(C),0,0)
#                 continue
#             end

#             m = length(before_incl_l)
#             F = base_ring(C)
#             ass[i,j,k,l] = zero_matrix(F,m,m)
#             for q ∈ 1:m, p ∈ 1:m
#                 ass[i,j,k,l][p,q] = F(after_proj_l[q] ∘ ass_mor ∘ before_incl_l[p])
#             end
    
#         end
#     end
#     return ass
# end

# struct SkeletizationFunctor <: AbstractFunctor
#     domain::Category
#     codomain::Category
# end

# function SkeletizationFunctor(C::Category)
#     SkeletizationFunctor(C,SixJCategory(C))
# end

# function (F::SkeletizationFunctor)(X::Object)
#     return SixJObject(codomain(F), [int_dim(Hom(X,s)) for s ∈ simples(parent(X))])
# end

# function (F::SkeletizationFunctor)(f::Morphism)
#     X = domain(f)
#     Y = codomain(f)

#     S = simples(domain(F))
#     n = length(S)

#     _, before, before_incl, before_proj = direct_sum_decomposition(X, S)
#     _, after,  after_incl,  after_proj  = direct_sum_decomposition(Y, S)

#     m = Vector{MatElem}(undef,n)

#     for l ∈ 1:n
#         before_incl_l = filter(f -> domain(f) == S[l], before_incl)
#         before_proj_l = filter(f -> codomain(f) == S[l], before_proj)
#         after_incl_l = filter(f -> domain(f) == S[l], after_incl)
#         after_proj_l = filter(f -> codomain(f) == S[l], after_proj)
        
#         if length(before_incl_l) == 0
#             m[l] = zero_matrix(base_ring(f),0,0)
#             continue
#         end
#         m[l] = matrix(direct_sum(after_proj_l) ∘ f ∘ direct_sum(before_incl_l))
#     end

#     return Morphism(F(X),F(Y),m)
# end

# function show(io::IO, F::SkeletizationFunctor)
#     print(io,"Skeletal projection of $(domain(F))")
# end


#=----------------------------------------------------------
    Skeletization Functor 
----------------------------------------------------------=#

