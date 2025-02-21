function six_j_category(C::Category, names::Vector{String} = ["X$i" for i ∈ 1:length(simples(C))])
    F = six_j_category(simples(C), names)
    set_name!(F, "Skeletization of $C")
    return F
end

function six_j_category(S::Vector{<:Object}, names::Vector{String} = ["X$i" for i ∈ 1:length(S)])
    C = parent(S[1])
    @assert is_ring(C)

    if typeof(C) == SixJCategory 
        return C
    end
    
    #S = simples(C)
    n = length(S)
    F = base_ring(C)

    mult = Array{Int,3}(undef,n,n,n)

    # prods = [X ⊗ Y for X ∈ S, Y ∈ S]
    # homs = [basis(Hom(prods[i,j],Z)) for i ∈ 1:length(S), j ∈ 1:length(S), Z ∈ S]
    # for i ∈ 1:n, j ∈ 1:n
    #     mult[i,j,:] = [length(homs[i,j,k]) for k ∈ 1:n]
    # end

    # Define SixJCategory
    skel_C = six_j_category(F,names)

    # Extract 6j-Symbols
    ass = six_j_symbols(C, S)

    # Recover multiplication table 
    one_index = findfirst(s -> int_dim(Hom(one(C),s)) > 0, S)
    set_one!(skel_C, [i == one_index for i ∈ 1:n])
    mult = [size(ass[i,j,one_index,k], 1) for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n]

    set_tensor_product!(skel_C,mult)

    set_associator!(skel_C, ass)

    if multiplicity(C) > 1
        @warn "6j-Symbols might be wrong since multiplicity is greater than one"
    end
    # Try to set spherical
    try 
        set_spherical!(skel_C,[F(1) for s ∈ S])
        sp = [dim(S[i]) * inv(dim(skel_C[i])) for i ∈ 1:length(S)]
        set_spherical!(Skel_C, sp)
    catch
    end

    if is_braided(C)
        set_braiding!(skel_C, skeletal_braiding(C,S))
    end
    

    return skel_C
end

function six_j_symbols(C::Category, S = simples(C))
    @assert is_semisimple(C)

    N = length(S)
    C_morphism_type = morphism_type(C)
    F = base_ring(C) 

    ass = Array{MatElem}(undef,N,N,N,N)

    one_indices = findall(s -> int_dim(Hom(s,one(C))) > 0 , S)
    one_components = simple_subobjects(one(C))

    # Set unitors to identity
    S[one_indices] = one_components

    prods = [X ⊗ Y for X ∈ S, Y ∈ S]
    homs = [basis(Hom(prods[i,j],Z)) for i ∈ 1:N, j ∈ 1:N, Z ∈ S]
    associators = [associator(X,Y,Z) for X ∈ S, Y ∈ S, Z ∈ S]



    for (i,j,k,l) ∈ Base.product(1:N, 1:N, 1:N, 1:N)

        #set trivial associators
        if !isempty([i,j,k] ∩ one_indices) 
            n = sum([length(homs[i,j,v]) * length(homs[v,k,l]) for v ∈ 1:N])
            ass[i,j,k,l] = identity_matrix(F,n)
            continue
        end

        # Build a basis for Hom((X⊗Y)⊗Z,W)
        B_XY_Z_W = C_morphism_type[]
        for n ∈ 1:N
            V = S[n]

            H_XY_V = homs[i,j,n]

            H_VZ_W = homs[n,k,l]

            B = [f ∘ (g ⊗ id(S[k])) for f ∈ H_VZ_W, g ∈ H_XY_V][:]
            if length(B) == 0 continue end
            B_XY_Z_W = [B_XY_Z_W; B]
        end

        # Build a basis for Hom(X⊗(Y⊗Z),W)
        B_X_YZ_W = C_morphism_type[]
        for n ∈ 1:N
            V = S[n]

            H_YZ_V = homs[j,k,n]

            H_XV_W = homs[i,n,l]
                    
            B = [f ∘ (id(S[i]) ⊗ g) for f ∈ H_XV_W, g ∈ H_YZ_V][:]
            if length(B) == 0 continue end
            B_X_YZ_W = [B_X_YZ_W; B]
        end
        
        # Express the asociator in the corresponding basis
        a = associators[i,j,k]

        associator_XYZ_W = hcat([express_in_basis(f ∘ a, B_XY_Z_W) for f ∈ B_X_YZ_W]...)
       
        ass[i,j,k,l] = matrix(F, length(B_XY_Z_W), length(B_X_YZ_W),  associator_XYZ_W)
            
    end

    

    return ass           
end

function skeletal_spherical(C::Category, Homs)
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


#=----------------------------------------------------------
    Gauge Transform  
----------------------------------------------------------=#

function gauge_transform(C::SixJCategory, bases::Dict)

    S = simples(C)
    F = base_ring(C)
    N = length(S)

    old_bases = [basis(Hom(x⊗y,z)) for x ∈ S, y ∈ S, z ∈ S]
    new_bases = [(i,j,k) ∈ keys(bases) ? bases[(i,j,k)] : old_bases[i,j,k] for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N]

    new_ass = Array{MatElem,4}(undef,N,N,N,N)

    N = length(S)
    C_morphism_type = morphism_type(C)

    for (i,j,k,l) ∈ Base.product(1:N, 1:N, 1:N, 1:N)

        # Build a basis for Hom((X⊗Y)⊗Z,W)
        B_XY_Z_W = C_morphism_type[]
        for n ∈ 1:N
            V = S[n]

            H_XY_V = new_bases[i,j,n]

            H_VZ_W = new_bases[n,k,l]

            B = [f ∘ (g ⊗ id(S[k])) for f ∈ H_VZ_W, g ∈ H_XY_V][:]
            if length(B) == 0 continue end
            B_XY_Z_W = [B_XY_Z_W; B]
        end


        # Build a basis for Hom(X⊗(Y⊗Z),W)
        B_X_YZ_W = C_morphism_type[]
        for n ∈ 1:N
            V = S[n]

            H_YZ_V = new_bases[j,k,n]

            H_XV_W = new_bases[i,n,l]
                    
            B = [f ∘ (id(S[i]) ⊗ g) for f ∈ H_XV_W, g ∈ H_YZ_V][:]
            if length(B) == 0 continue end
            B_X_YZ_W = [B_X_YZ_W; B]
        end


        a = associator(S[[i,j,k]]...)
        m = hcat([express_in_basis(f ∘ a, B_XY_Z_W) for f ∈ B_X_YZ_W]...)
        new_ass[i,j,k,l] =  matrix(F, length(B_X_YZ_W), length(B_XY_Z_W), m)
    end

    new_C = deepcopy(C)
    new_C.ass = new_ass 
    new_C.base_ring = F
    return new_C 
end


function with_gauge_freedom(C::SixJCategory)
    K = base_ring(C)
    
    one_index = C.one 
    N = length(simples(C)) 

    # Take all combinations not fixed by unitors 
    indices = [(i,j) for i ∈ 1:N, j ∈ 1:N if isempty([i,j] ∩ one_index)]

    # Bases for Hom spaces Hom(ij,k)
    homs = [(i,j,k) => basis(Hom(C[i]⊗C[j], C[k])) for (i,j) ∈ indices, k ∈ 1:N if C.tensor_product[i,j,k] > 0]

    R,x = polynomial_ring(K, sum(length.(homs)))
    L = fraction_field(R)
    y = L.(x)
    homs = [k => [popfirst!(y) * (f ⊗ L) for f ∈ H] for (k,H) ∈ homs]

    return gauge_transform(C ⊗ L, Dict(homs))
end

function _subst_gauge!(C::SixJCategory, x::Vector{<:FracFieldElem})
    N = length(simples(C))
    L = base_ring(C)

    for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N, l ∈ 1:N
        C.ass[i,j,k,l] = matrix(L, size(C.ass[i,j,k,l])..., [numerator(a)(x...)//denominator(a)(x...) for a in C.ass[i,j,k,l]])
    end
end

function _subst_gauge!(C::SixJCategory, ind::Tuple{Int,Int,Int}, x::Vector{<:FracFieldElem})
    N = length(simples(C))
    L = base_ring(C)

    for l ∈ 1:N
        i,j,k = ind
        C.ass[i,j,k,l] = matrix(L, size(C.ass[i,j,k,l])..., [numerator(a)(x...)//denominator(a)(x...) for a in C.ass[i,j,k,l]])
    end
end

function _subst_gauge!(C::SixJCategory, ind::Vector{Int}, x::Vector{<:FracFieldElem})
    N = length(simples(C))
    L = base_ring(C)

    for i ∈ ind, j ∈ ind, k ∈ ind, l ∈ 1:N
        C.ass[i,j,k,l] = matrix(L, size(C.ass[i,j,k,l])..., [numerator(a)(x...)//denominator(a)(x...) for a in C.ass[i,j,k,l]])
    end
end