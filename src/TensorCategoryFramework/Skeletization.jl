function six_j_category(C::Category, names::Vector{String} = simples_names(C); multiplication_table = nothing)
    F = six_j_category(simples(C), names, multiplication_table = multiplication_table)
    set_name!(F, "Skeletization of $C")
    return F
end

function six_j_category(S::Vector{<:Object}, names::Vector{String} = simples_names(parent(S[1])); multiplication_table = nothing)
    C = parent(S[1])
    @assert is_ring(C)

    if typeof(C) == SixJCategory 
        return C
    end
    
    #S = simples(C)
    n = length(S)
    F = base_ring(C)


    # prods = [X ⊗ Y for X ∈ S, Y ∈ S]
    # homs = [basis(Hom(prods[i,j],Z)) for i ∈ 1:length(S), j ∈ 1:length(S), Z ∈ S]
    # for i ∈ 1:n, j ∈ 1:n
    #     mult[i,j,:] = [length(homs[i,j,k]) for k ∈ 1:n]
    # end

    # Define SixJCategory
    skel_C = six_j_category(F,names)

    # Extract 6j-Symbols
    ass = six_j_symbols(C, S, multiplication_table)

    # Recover multiplication table 
    one_index = findfirst(s -> int_dim(Hom(one(C),s)) > 0, S)
    set_one!(skel_C, [i == one_index for i ∈ 1:n])

    mult = multiplication_table !== nothing ? multiplication_table : [size(ass[i,j,one_index,k], 1) for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n]

    set_tensor_product!(skel_C,mult)

    set_associator!(skel_C, ass)

    # if multiplicity(C) > 1
    #     @warn "6j-Symbols might be wrong since multiplicity is greater than one"
    # end
    # Try to set spherical
    try 
        set_pivotal!(skel_C,[F(1) for s ∈ S])
        sp = [dim(S[i]) * inv(dim(skel_C[i])) for i ∈ 1:length(S)]
        set_pivotal!(Skel_C, sp)
    catch
    end

    if is_braided(C)
        set_braiding!(skel_C, skeletal_braiding(C,S))
    end
    

    return skel_C
end

function six_j_symbols(C::Category, S = simples(C), mult = nothing)
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

    if mult !== nothing 
        global homs = [mult[i,j,k] > 0 ? basis(Hom(prods[i,j],S[k])) : C_morphism_type[] for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N]
    else
        global homs = [basis(Hom(prods[i,j],S[k])) for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N]
    end

    associators = [associator(X,Y,Z) for X ∈ S, Y ∈ S, Z ∈ S]



    @threads for (i,j,k,l) ∈ collect(Base.product(1:N, 1:N, 1:N, 1:N))

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

function six_j_symbols_of_construction(C::Category, S = simples(C), mult = nothing; log = nothing)
    @assert is_semisimple(C)

    if log !== nothing && !isdir(joinpath(@__DIR__, "$log"))
        mkdir(joinpath(@__DIR__, "$log"))
    end 

    N = length(S)
    C_morphism_type = morphism_type(category(C))
    F = base_ring(C) 

   # mult = multiplication_table(C)

    ass = Array{MatElem}(undef,N,N,N,N)

    one_indices = findall(s -> int_dim(Hom(s,one(C))) > 0 , S)
    one_components = simple_subobjects(one(C))

    # Set unitors to identity
    S[one_indices] = one_components

    # prods = [X ⊗ Y for X ∈ S, Y ∈ S]

    # if mult !== nothing 
    #     global homs = [mult[i,j,k] > 0 ? morphism.(basis(Hom(prods[i,j],S[k]))) : C_morphism_type[] for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N]
    # else
    #     global homs = [morphism.(basis(Hom(prods[i,j],S[k]))) for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N]
    # end  
    
    homs = multiplicity_spaces(C)
    homs = Dict(k => morphism.(basis(v)) for (k,v) in homs)
    missed = [(i,j,k) => C_morphism_type[] for i in 1:N, j in 1:N, k in 1:N if !haskey(homs, (i,j,k))]
    if length(missed) > 0
        push!(homs, missed...)
    end
    #associators = [associator(object.((X,Y,Z))...) for X ∈ S, Y ∈ S, Z ∈ S]



    @threads for i ∈ 1:N 
        for j ∈ 1:N, k ∈ 1:N
            a = associator(object.(S[[i,j,k]])...)

            for l ∈ 1:N
                if log !== nothing && isfile(joinpath(@__DIR__, "$log/$(i)_$(j)_$(k)_$(l)"))
                    _m =  load(joinpath(@__DIR__, "$log/$(i)_$(j)_$(k)_$(l)"))
                    ass[i,j,k,l] = matrix(F, size(_m,1),size(_m,2), collect(_m)) 
                    continue 
                end

                #set trivial associators
                if !isempty([i,j,k] ∩ one_indices) 
                    n = sum([length(homs[i,j,v]) * length(homs[v,k,l]) for v ∈ 1:N])
                    ass[i,j,k,l] = identity_matrix(F,n)

                    if log !== nothing 
                    save( joinpath(@__DIR__, "$log/$(i)_$(j)_$(k)_$(l)"),ass[i,j,k,l])
                end
                    continue
                end

                # Build a basis for Hom((X⊗Y)⊗Z,W)
                B_XY_Z_W = C_morphism_type[]
                for n ∈ 1:N
                    V = S[n]

                    H_XY_V = homs[i,j,n]
                    H_VZ_W = homs[n,k,l]

                    B = [f ∘ (g ⊗ id(object(S[k]))) for f ∈ H_VZ_W, g ∈ H_XY_V][:]
                    if length(B) == 0 continue end
                    B_XY_Z_W = [B_XY_Z_W; B]
                end

                # Build a basis for Hom(X⊗(Y⊗Z),W)
                B_X_YZ_W = C_morphism_type[]
                for n ∈ 1:N
                    V = S[n]

                    H_YZ_V = homs[j,k,n]

                    H_XV_W = homs[i,n,l]
                            
                    B = [f ∘ (id(object(S[i])) ⊗ g) for f ∈ H_XV_W, g ∈ H_YZ_V][:]
                    if length(B) == 0 continue end
                    B_X_YZ_W = [B_X_YZ_W; B]
                end
                

                # Express the asociator in the corresponding basis
                #a = associators[i,j,k]
                associator_XYZ_W = hcat([express_in_basis(f ∘ a, B_XY_Z_W) for f ∈ B_X_YZ_W]...)
            
                ass[i,j,k,l] = matrix(F, length(B_XY_Z_W), length(B_X_YZ_W),  associator_XYZ_W) 

                if log !== nothing && !isempty(associator_XYZ_W)
                    @show (i,j,k,l)
                    save(joinpath(@__DIR__, "$log/$(i)_$(j)_$(k)_$(l)"),associator_XYZ_W)
                end
            end
        end
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

    homs = multiplicity_spaces(C) 
    homs = Dict(k => (basis(v)) for (k,v) in homs)
    missed = [(i,j,k) => C_morphism_type[] for i in 1:N, j in 1:N, k in 1:N if !haskey(homs, (i,j,k))]
    if length(missed) > 0
        push!(homs, missed...)
    end

    Threads.@threads for i ∈ 1:N
        for j ∈ 1:N, l ∈ 1:N
            X,Y,W = S[[i,j,l]]
            # Basis for Hom(X⊗Y,W)
            B_XY_W = (homs[i,j,l])

            # Basis for Hom(Y⊗X,W)
            B_YX_W = (homs[j,i,l])

            braid_XY_W = hcat([express_in_basis(f ∘ braiding(X,Y), B_XY_W) for f ∈ B_YX_W]...)
            braid[i,j,l] = matrix(F, length(B_XY_W), length(B_YX_W), braid_XY_W)
        end
    end
   return braid
end

function skeletal_braiding_of_construction(C::Category, S = simples(C), mult = nothing)
    @assert is_braided(C)
    
    N = length(S)
    C_morphism_type = morphism_type(C)
    F = base_ring(C) 
    braid = Array{MatElem}(undef,N,N,N)

    homs = multiplicity_spaces(C) 
    homs = Dict(k => (basis(v)) for (k,v) in homs)
    missed = [(i,j,k) => C_morphism_type[] for i in 1:N, j in 1:N, k in 1:N if !haskey(homs, (i,j,k))]
    if length(missed) > 0
        push!(homs, missed...)
    end

    for (i,j,l) ∈ Base.product(1:N,1:N,1:N)
        X,Y,W = S[[i,j,l]]
        # Basis for Hom(X⊗Y,W)
        B_XY_W = morphism.(basis(homs[i,j,l]))

        # Basis for Hom(Y⊗X,W)
        B_YX_W = morphism.(basis(homs[j,i,l]))

        braid_XY_W = hcat([express_in_basis(f ∘ morphism(braiding(X,Y)), B_XY_W) for f ∈ B_YX_W]...)
        braid[i,j,l] = matrix(F, length(B_XY_W), length(B_YX_W), braid_XY_W)
    end
    return braid
end

#=----------------------------------------------------------
    Gauge Transform  
----------------------------------------------------------=#

# function unitary_gauge(C::SixJCategory)
#     @assert multiplicity(C) == 1
#     N = rank(C) 
#     K = base_ring(C)
    
#     if  ! is_normal(K) 
#         K,emb = normal_closure(K)
#         C = extension_of_scalars(C,K,embedding = emb)
#     end

#     conjg = complex_conjugation(K)

    

#     D = six_j_category(R, multiplication_table(C))
#     D.one = C.one 
#     D.ass = [change_base_ring(R,m) for m in C.ass]
    
#     trafo = Dict((i,j,k) => popfirst!(g) * Hom(D[]))
# end

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

    non_trivial_indices = findall(!=(0), multiplication_table(C))

    R,g = polynomial_ring(K, sum(multiplication_table(C)))
    L = fraction_field(R)
    D = six_j_category(L, multiplication_table(C))
    D.one = C.one 
    D.ass = [change_base_ring(L,m) for m in C.ass]

    # Bases for Hom spaces Hom(ij,k)
    @show homs = [(i,j,k) => [popfirst!(g) * f for f ∈ basis(Hom(D[i]⊗D[j], D[k]))] for (i,j,k) ∈ Tuple.(non_trivial_indices)]

    return gauge_transform(D, Dict(homs))
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


#=----------------------------------------------------------
    Export Skeletal categories as dicts 
----------------------------------------------------------=#

function save_F_symbols(C::SixJCategory, file::String)  

    F = F_symbols(C)
    S,T = typeof(F).parameters
    pol = polynomial(QQ,collect(coefficients(base_ring(C).pol)))

    open(file, "w") do io
        write(io, "# Field with defining polynomial $pol \n# relative to the basis 1,...,x^$(degree(pol)-1)\n\n ")

        write(io,"Dict{$S,$T}(")

        write(io, join(["\t$k => $(coefficients(v))" for (k,v) in F], ",\n"))

        write(io, "\n)")
    end
    nothing
end

function save_R_symbols(C::SixJCategory, file::String)  

    R = R_symbols(C)
    S,T = typeof(R).parameters
    pol = polynomial(QQ,collect(coefficients(base_ring(C).pol)))

    open(file, "w") do io
        write(io, "# Field with defining polynomial $pol \n# relative to the basis 1,...,x^$(degree(pol)-1)\n\n ")
        write(io,"Dict{$S,$T}(")

        write(io, join(["\t$k => $(coefficients(v))" for (k,v) in R], ",\n"))

        write(io, "\n)")
    end
    nothing
end

function save_P_symbols(C::SixJCategory, file::String)  

    P = P_symbols(C)
    S,T = typeof(P).parameters
    pol = polynomial(QQ,collect(coefficients(base_ring(C).pol)))

    open(file, "w") do io
        write(io, "# Field with defining polynomial $pol \n# relative to the basis 1,...,x^$(degree(pol)-1)\n\n ")
        write(io,"Dict{$S,$T}(")

        write(io, join(["\t$k => $(coefficients(v))" for (k,v) in P], ",\n"))

        write(io, "\n)")
    end
    nothing
end