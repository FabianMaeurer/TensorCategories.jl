#=----------------------------------------------------------
    Load fusion categories from the Anyon Wiki by 
    Gert Vercleyen. 
----------------------------------------------------------=#

associator_path = joinpath(@__DIR__, "AnyonWikiData//")
pivotal_path = joinpath(@__DIR__, "AnyonWikiData/PivotalStructures/")
anyon_path = joinpath(artifact"AnyonWiki", "AnyonWiki")




@doc raw""" 

    anyonwiki(r,m,n,i,a,b,p)

Load the fusion category from the list of multiplicity free fusion categories of rank ≤ 7 with index (r,m,n,i,a,b,p).
"""
function anyonwiki(rank::Int, 
                    multiplicity::Int, 
                    non_self_dual::Int,
                    fusion_ring::Int, 
                    associator::Int, 
                    braiding::Int, 
                    pivotal::Int)

    K,rt = load_anyonwiki_number_field(rank,multiplicity,non_self_dual,fusion_ring,associator,braiding, pivotal)

    cat_code = "$(rank)_$(multiplicity)_$(non_self_dual)_$(fusion_ring)_$(associator)_$(braiding)_$(pivotal)"
    cat_string = "cat_$cat_code.jl"

    C = six_j_category(K, ["𝟙"; String["X$i" for i in 2:rank]])

    ass = include(joinpath(anyon_path, "algebraic_F_symbols/$cat_string"))

    ass = Dict(k => K == QQ ? K(v...) : K(v) for (k,v) in ass)

    ass = dict_to_associator(rank, K, ass)

    set_tensor_product!(C, multiplication_table_from_F_symbols(ass))
    set_associator!(C, transpose.(ass))
    set_one!(C, [i == 1 for i in 1:rank])

    if braiding != 0 
        braid = include(joinpath(anyon_path, "algebraic_R_symbols/$cat_string"))
        braid = Dict(k =>  K == QQ ? K(v...) : K(v) for (k,v) in braid)
        braid = dict_to_braiding(rank, K, braid)
        set_braiding!(C, braid)
    end

    piv = include(joinpath(anyon_path, "algebraic_P_symbols/$cat_string"))

    piv = Dict(k =>  K == QQ ? K(v...) : K(v) for (k,v) in piv)

    piv = [K(piv[[p]]) for p in 1:rank]
    set_pivotal!(C, piv)

    setfield!(C, :embedding, rt)
    set_name!(C, "Fusion Category $cat_code")
    return C
end

@doc raw""" 

    anyonwiki_center(i,j,k,l,m,n,o)

Return the center of the fusion category with index (i,j,k,l,m,n,o) from the database.
"""
function anyonwiki_center(i,j,k,l,m,n,o)
    path = anyonwiki_center_artifact_path(i,j,k,l,m,n,o)

    C = load_fusion_category(path)
    set_name!(C, replace(C.name, "Skeletization" => "Skeletonization"))
    C
end

function anyonwiki_center_artifact_path(i,j,k,l,m,n,o)
    try 
        if i ≤ 4 
            path = @artifact_str "AnyonWikiCenter1to4"
            path = joinpath(path, "center_$(i)_$(j)_$(k)_$(l)_$(m)_$(n)_$(o)")
            open(path)
            return path

        elseif i == 5 
            path = @artifact_str "AnyonWikiCenter5"
            path = joinpath(path, "center_$(i)_$(j)_$(k)_$(l)_$(m)_$(n)_$(o)")
            open(path)
            return path
        end
    catch 
        error("There is no center saved for a fusion category with index $((i,j,k,l,m,n,o))")
    end
end

function dict_to_associator(ass::Dict)
    N = length(filter(e -> all(e[[1,2]] .== 1), keys(ass)))
    dict_to_associator(N, parent(first(ass)[2]), ass)
end

function dict_to_associator(N::Int, K::Field, ass::Dict)
    # Transform associator dict to Matrices 

    ass_matrices = Array{MatElem,4}(undef,N,N,N,N)

    groups = group_dict_keys_by(e -> e[1:4], ass)

    for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N, d ∈ 1:N 
        if !haskey(groups, [a,b,c,d])
            ass_matrices[a,b,c,d] = zero_matrix(K,0,0)
            continue
        end
        D = groups[[a,b,c,d]]
        abc_d = collect(keys(D))
        l = Int(sqrt(length(abc_d)))
        if length(first(keys(ass))) == 6 
            abc_d = sort(abc_d, by = v -> v[6])
            abc_d = sort(abc_d, by = v -> v[5])
        else
            abc_d = sort(abc_d, by = v -> v[[8,5,10,9,7,6]])
        end
        M = matrix(K,l,l, [D[v] for v ∈ abc_d])
        ass_matrices[a,b,c,d] = transpose(M)
    end

    ass_matrices 
end

function anyonwiki_keys(n::Int = 7)
    d = include(joinpath(@__DIR__, "base_field_generators.jl"))

    return [k for (k,v) in sort(d) if k[1] ≤ n]
end

function group_dict_keys_by(f::Function, D::Dict)
    groups = Dict()
    for (k,v) ∈ D 
        f_k = f(k)
        if f_k ∈ keys(groups)
            push!(groups[f_k], k => v)
        else 
            push!(groups, f_k => Dict(k => v))
        end
    end
    return groups 
end

function dict_to_braiding(ass::Dict)
    N = length(filter(e -> all(e[1] == 1), keys(ass)))
    dict_to_braiding(N, parent(first(ass)[2]), ass)
end

function dict_to_braiding(N::Int, K::Field, braid::Dict)
    # Transform associator dict to Matrices 
    braiding_array = Array{MatElem,3}(undef,N,N,N)

    for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N
        ab_c = filter(e -> e[[1,2,3]] == [a,b,c], collect(keys(braid)))
        l = Int(sqrt(length(ab_c)))
        sort!(ab_c)
        M = matrix(K,l,l, [braid[v] for v ∈ ab_c])
        braiding_array[a,b,c] = transpose(M)
    end
    braiding_array
end

function anyonwiki_from_numerics(rank::Int, 
                    multiplicity::Int, 
                    non_self_dual::Int,
                    fusion_ring::Int, 
                    associator::Int, 
                    braiding::Int, 
                    pivotal::Int, 
                    K::Union{Ring,Nothing} = nothing ; minimal = false)
    
    field_dict = include(joinpath(@__DIR__, "base_field_generators.jl"))

    # if ! minimal && K === nothing 
    #     if n ∈ ZZ.(readlines(joinpath(anyon_path, "cycloPos.txt")))
    #         return anyonwiki_cyclotomic(n)
    #     end
    # end

    K,e = load_anyonwiki_number_field(rank, multiplicity, non_self_dual,fusion_ring, associator, braiding, pivotal)


    data = []
    CC = AcbField(2048)
    Q = QQBarField()
    # deg = typeof(K) <: Union{<:NumField,QQField} ? degree(K) : 24
    # if typeof(K) <: Union{<:NumField,QQField}
    #     es = complex_embeddings(K) 
    #     r = load_anyonwiki_associator_root(n)
    #     r2 = load_anyonwiki_pivotal_root(n) 
    #     for f ∈ es 
    #         try 
    #             preimage(f,r)
    #             preimage(f,r2)
    #             global e = f
    #             break 
    #         catch 
    #         end
    #     end
    # end

    cat_string = "cat_$(rank)_$(multiplicity)_$(non_self_dual)_$(fusion_ring)_$(associator)_$(braiding)_$(pivotal).jl"
    
    # Get associator dict with strings
    _ass_dict = include(joinpath(@__DIR__, "AnyonWikiData/new_6j_Symbols/$cat_string"))
    # ass_dict = Dict{Vector{Int}, elem_type(K)}()
    # # Convert data to numberfield
    # for (k, str) ∈ _ass_dict 
    #     x = string_to_acb(CC, str)
    #     x = preimage(e, x, degree(K))
    #     ass_dict[k] = K(x)
    # end

    ass_dict = apply_preimage_to_anyon_file(e, joinpath(@__DIR__, "AnyonWikiData/new_6j_Symbols/$cat_string"))

    # open(joinpath(anyon_path, "new_6j_symbols/cat_$n")) do f 

    #     while ! eof(f)

    #         s = readline(f) 
    #         s = filter(!=('\"'), s)
    #         s = split(s, " ")

    #         indices = Int[eval(Meta.parse(a)) for a ∈ s[1:end-3]]

    #         found = false 

    #         while ! found
    #             x_real = CC(s[end-2] * "+/- 1e-510")
    #             x_imag = CC(s[end][1:end-2] * "+/- 1e-510") * CC(im) 
    #             try 
    #             global x = guess(Q, x_real + x_imag, deg)
    #             found = true
    #             catch  err
    #                 if precision(CC) > 5000 error(err) end
    #                 CC = AcbField(2*precision(CC)) 
    #             end
    #         end
    #         if typeof(K) <: Union{<:NumField,QQField}
    #             append!(data, [Any[indices; preimage(e,x)]])
    #         else 
    #             append!(data, [Any[indices; K(x)]])       
    #         end
    #     end
    # end
    N = rank
    
    # Construct multiplication table from associators
    mult = zeros(Int,N,N,N)
    
    for (_,a,b,c) ∈ filter(e -> e[1] == 1, keys(ass_dict)) 
        
        mult[a,b,c] = 1 
    end

    C = six_j_category(K, mult)
    set_name!(C, "Fusion category ($(cat_string[5:end-3])) from AnyonWiki")

    # Build matrices for 6j-Symbols
    for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N, d ∈ 1:N
        abc_d = filter(e -> e[[1,2,3,4]] == [a,b,c,d], collect(keys(ass_dict)))
        l = Int(sqrt(length(abc_d)))
        abc_d = sort(abc_d, by = v -> v[6])
        abc_d = sort(abc_d, by = v -> v[5])
        
        M = matrix(K,l,l, [ass_dict[v] for v ∈ abc_d])
        set_associator!(C,a,b,c,d,M)  
    end
    set_one!(C, 1)
    

    set_pivotal!(C, load_anyonwiki_pivotal(rank::Int, 
        multiplicity::Int, 
        non_self_dual::Int,
        fusion_ring::Int, 
        associator::Int, 
        braiding::Int, 
        pivotal::Int, K, e))
    
    if braiding != 0
        set_braiding!(C, anyonwiki_braiding(rank::Int, 
            multiplicity::Int, 
            non_self_dual::Int,
            fusion_ring::Int, 
            associator::Int, 
            braiding::Int, 
            pivotal::Int, K, e))
    end
    C
end

function string_to_acb(CC::AcbField, str::String)
    re,co = split(str, "+")
    if co == "0"
        x = CC(re * "+/- 2e-510")
    else
        x = CC(re * "+/- 2e-510") + CC(co[1:end-2] * "+/- 2e-510")*CC(im)
    end
end

function load_anyonwiki_number_field(rank::Int, 
    multiplicity::Int, 
    non_self_dual::Int,
    fusion_ring::Int, 
    associator::Int, 
    braiding::Int, 
    pivotal::Int)

    field_dict = include(joinpath(@__DIR__, "base_field_generators.jl"))

    data = field_dict[[rank, multiplicity, non_self_dual, fusion_ring, associator, braiding, pivotal]]

    if typeof(data) == Int 
        if data == 0 
            return QQ, complex_embedding(rationals_as_number_field()[1], 1)
        end

        K,z = cyclotomic_field(data, "z$(data)")

        emb = complex_embeddings(K)[1]

        return K, emb
    else 
        p,str = data
        K,a = number_field(polynomial(QQ,p))
        
        CC = AcbField(2048)
        root = string_to_acb(CC,str)
        emb = complex_embedding(K,root)
        return K, emb
    end
    # line = 1
    # open(joinpath(anyon_path, "base_field_generators.jl")) do f 

    #     while ! eof(f) 

    #         s = readline(f) 
            
    #         if line == n 
    #             global _,x = QQ[:x]
    #             min_pol = eval(Meta.parse(s))

    #             global K1 = min_pol == 1 ? QQ : number_field(min_pol)[1]
    #             break
    #         end
            
    #         line += 1
    #     end
    # end
    
    # line = 1
    # open(joinpath(pivotal_path, "Polynomials")) do f 

    #     while ! eof(f) 

    #         s = readline(f) 
            
    #         if line == n 
    #             global _,x = QQ[:x]
    #             min_pol = eval(Meta.parse(s))

    #             global K2 =  min_pol == 1 ? QQ : number_field(min_pol)[1]
    #             break
    #         end
            
    #         line += 1
    #     end
    # end
    
    # K1 == K2 && return K1
    # K1 == QQ && return K2
    # K2 == QQ && return K1 

    # simplify(splitting_field([defining_polynomial(K1), defining_polynomial(K2)]))[1]
end

# function load_anyonwiki_associator_root(n::Int)
#     line = 1
#     CC = AcbField(2048)
#     Q = QQBarField()
#     deg = 24
#     open(joinpath(associator_path, "Roots.dat")) do f 

#         while ! eof(f) 
#             s = readline(f) 
#             if line == n
#                 s = filter(!=('\"'), s)
#                 s = split(s, " ")
#                 x_real = CC(s[1] * " +/- 1e-512")
#                 x_imag = CC(s[end][1:end-2] * "+/- 1e-512") * CC(im) 

#                 return guess(Q,x_real + x_imag,deg)
#             end
#             line += 1
#         end
#     end
# end

# function load_anyonwiki_pivotal_root(n::Int)
#     line = 1
#     CC = AcbField(2048)
#     Q = QQBarField()
#     deg = 24
#     open(joinpath(pivotal_path, "Roots.dat")) do f 

#         while ! eof(f) 
#             s = readline(f) 
#             if line == n
#                 s = filter(!=('\"'), s)
#                 s = split(s, " ")
#                 x_real = CC(s[1] * " +/- 1e-512")
#                 x_imag = CC(s[end][1:end-2] * "+/- 1e-512") * CC(im) 

#                 return guess(Q,x_real + x_imag,deg)
#             end
#             line += 1
#         end
#     end
# end

function apply_preimage_to_anyon_file(e::AbsSimpleNumFieldEmbedding, file::String)
    str = read(file, String)
    reg = r"\"[^\"]*\""
    matches = unique(collect([m.match[2:end-1] for m in eachmatch(reg,str)]))
    sub = Dict(String(k) =>  preimage(e, string_to_acb(AcbField(2048), String(k)), degree(number_field(e))) for k in matches)
    dict = include(file)

    new_dict = Dict(k => sub[v] for (k,v) ∈ dict)
end

function load_anyonwiki_pivotal(rank::Int, 
    multiplicity::Int, 
    non_self_dual::Int,
    fusion_ring::Int, 
    associator::Int, 
    braiding::Int, 
    pivotal::Int,
    K, e = complex_embeddings(K)[1])

    CC = AcbField(2048)
    Q = QQBarField()

    cat_string = "cat_$(rank)_$(multiplicity)_$(non_self_dual)_$(fusion_ring)_$(associator)_$(braiding)_$(pivotal).jl"

    piv = include(joinpath(@__DIR__, "AnyonWikiData/new_P_Symbols/$cat_string"))

    [K(preimage(e,string_to_acb(CC, piv[x]), degree(K))) for x in sort(collect(keys(piv)))]
    # deg = typeof(K) <: Union{<:NumField,QQField} ? degree(K) : 24

    # piv = elem_type(K)[]
    # open(joinpath(pivotal_path, "pivots_cat_$n")) do f 

    #     while ! eof(f) 

    #         s = readline(f) 
    #         s = filter(!=('\"'), s)
    #         s = split(s, " ")

    #         x_real = CC(s[2] * " +/- 1e-512")
    #         x_imag = CC(s[4][1:end-2] * "+/- 1e-512") * CC(im) 

    #         x = guess(Q, x_real + x_imag, deg)

    #         if typeof(K) <: Union{<:NumField,QQField}
    #             push!(piv, preimage(e, x))
    #         else 
    #             push!(piv, K(x))    
    #         end
    #     end
    # end
    # return piv 
end


function load_anyonwiki_attributes(n::Int) 
    symbols = [ :braided, :unitary, :spherical, :ribbon, :modular]
    open(joinpath(@__DIR__, "AnyonWikiData/cat_properties.txt")) do f 

        line = 1
        while ! eof(f) 
            s = readline(f) 

            if line == n+1
                s = filter(!=('\"'), s)
                s = split(s, ", ")

                return Dict(k => parse(Bool, t) for (t, k) in zip(s,symbols))
            end
            line += 1
        end
    end
end


function anyonwiki_cyclotomic(n::Int)
    r = anyonwiki_cyclotomic_root(n)
    
    if r > 2 
        eval(Meta.parse("global K, z$r = cyclotomic_field($r, \"z$r\")"))
    else
        global K = QQ 
    end
    
    anyonwiki_cyclotomic(n, K)
end

# function anyonwiki_cyclotomic(n::Int, K::Field)
    
#     r = anyonwiki_cyclotomic_root(n)
#     root = root_of_unity(K,r)
#     eval(Meta.parse("z$r = $(root)"))


#     associator_array = anyonwiki_cyclotomic_associator(n,K)
   
#     mult = multiplication_table_from_F_symbols(associator_array)


#     C = six_j_category(K, mult)
#     set_name!(C, "Fusion category number $n from AnyonWiki")
#     set_associator!(C, associator_array)
#     set_one!(C, [i == 1 for i ∈ 1:C.simples])

#     pivotals = include(joinpath(anyon_path, "cyclic_P_Symbols/cat_$n.jl"))

#     set_pivotal!(C, K.([pivotals[[i]] for i ∈ 1:length(pivotals)]))

#     attrs = load_anyonwiki_attributes(n)

#     if attrs[:braided] 
       
#         braiding_array = anyonwiki_cyclotomic_braiding(n, K)
#         set_braiding!(C, braiding_array)
#     end

#     return C
    
# end

function anyonwiki_finite(K::FqField,i,j,k,l,m,n,o)
    C = anyonwiki(i,j,k,l,m,n,o)

    extension_of_scalars(C,K)
    # r = anyonwiki_cyclotomic_root(n)

    # f_symbol_string = readlines(joinpath(anyon_path, "cyclic_6j_Symbols/cat_$n.jl"))[4]

    # ind = findall(r"\//\d+", f_symbol_string)
    # inverse = lcm([parse(Int, f_symbol_string[i[3:end]]) for i in ind])

    # K = finite_prime_field_with_root_of_unity(r, inverse + 1)

    # anyonwiki_cyclotomic(n,K)
end

function multiplication_table_from_F_symbols(ass::Array{<:MatElem,4})
    # Build multiplication_table
    N, = size(ass)

    mult = zeros(Int,N,N,N)
    
    for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N 
        mult[i,j,k] = size(ass[1,i,j,k])[1]
    end
    return mult
end

# function anyonwiki_braiding(rank::Int, 
#     multiplicity::Int, 
#     non_self_dual::Int,
#     fusion_ring::Int, 
#     associator::Int, 
#     braiding::Int, 
#     pivotal::Int,
#     K::Field, emb)

#     CC = AcbField(2048)

#     cat_string = "cat_$(rank)_$(multiplicity)_$(non_self_dual)_$(fusion_ring)_$(associator)_$(braiding)_$(pivotal).jl"
    
#     _braiding_dict = include(joinpath(anyon_path, "new_R_Symbols/$cat_string"))
#     braiding_dict = Dict{Vector{Int}, elem_type(K)}()
#     # for (k, str) ∈ _braiding_dict 
#     #     x = string_to_acb(CC, @show str)
#     #     x = preimage(emb, x, degree(K))
#     #     braiding_dict[k] = K(x)
#     # end

#     braiding_dict = apply_preimage_to_anyon_file(emb, joinpath(anyon_path, "new_R_Symbols/$cat_string"))
#     N = rank

#     braiding_array = Array{MatElem}(undef, N,N,N)

#     for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N
#         ab_c = filter(e -> e[[1,2,3]] == [a,b,c], collect(keys(braiding_dict)))
#         l = Int(sqrt(length(ab_c)))

#         M = matrix(K,l,l, [braiding_dict[v] for v ∈ ab_c])
#         braiding_array[a,b,c] = M
#     end

#     return braiding_array
# end

# function anyonwiki_cyclotomic_root(n::Int)
#     cyclopos = [parse(Int, s) for s ∈ readlines(joinpath(anyon_path, "cycloPos.txt"))]

#     attrs = load_anyonwiki_attributes(n)

#     if n ∉ cyclopos 
#         error("Category number $n is not cyclotomic")
#     end

#     f_symbol_string = readlines(joinpath(anyon_path, "cyclic_6j_Symbols/cat_$n.jl"))[4]
    
#     p_symbol_string = readlines(joinpath(anyon_path, "cyclic_P_Symbols/cat_$n.jl"))[4]

#     r_symbol_string = if attrs[:braided]
#         readlines(joinpath(anyon_path, "cyclic_R_Symbols/cat_$n.jl"))[4]
#     else 
#         ""
#     end

#     # read root
#     i = findfirst(r"z\d*", f_symbol_string)
#     j = findfirst(r"z\d*", p_symbol_string)
#     k = if attrs[:braided]
#         findfirst(r"z\d*", r_symbol_string)
#     else 
#         nothing
#     end 

#     if i === nothing && j === nothing && k == nothing
#         return 1 
#     else
#         roots = [l === nothing ? 1 : parse(Int, s[l][2:end]) for (s,l) ∈ zip([f_symbol_string,p_symbol_string, r_symbol_string], [i,j,k])]

#         return r = lcm(roots...)
#     end
# end

# function anyonwiki_cyclotomic_associator(n::Int, K::Field)
#     r = anyonwiki_cyclotomic_root(n)
#         root = root_of_unity(K,r)
#     eval(Meta.parse("z$r = $(root)"))

#     # read 6j-symbols 
#     associator_dict = include(joinpath(anyon_path, "cyclic_6j_Symbols/cat_$n.jl"))

#     N = maximum([a[1] for a ∈ keys(associator_dict)])

#     associator_array = Array{MatElem}(undef, N,N,N,N)
    
#     # Build matrices for 6j-Symbols
#     for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N, d ∈ 1:N
#         abc_d = filter(e -> e[[1,2,3,4]] == [a,b,c,d], collect(keys(associator_dict)))
#         l = Int(sqrt(length(abc_d)))
#         abc_d = sort(abc_d, by = v -> v[6])
#         abc_d = sort(abc_d, by = v -> v[5])
        
#         M = matrix(K,l,l, [associator_dict[v] for v ∈ abc_d])
#         associator_array[a,b,c,d] = M
#     end
#     return associator_array 
# end

function finite_prime_field_with_root_of_unity(n::Int, lower_bound = 2)
    p = next_prime(maximum([n,lower_bound-1])) 
    while gcd(n, p - 1) < n
        p = next_prime(p)
    end
    return GF(p)
end

function anyonwiki_center(n::Int; cyclotomic = false)
    if cyclotomic 
        load(joinpath(@__DIR__, "AnyonWikiData/CyclotomicCenters/CyclotomicCenter_$n.mrdi"))
    else 
        load(joinpath(@__DIR__, "AnyonWikiData/Centers/Center_$n.mrdi"))
    end
end


# function load_anyonwiki_fusion_rules(n::Int)

#     data = []

#     open(joinpath(associator_path, "cat_$n")) do f 

#         while ! eof(f)

#             s = readline(f) 
#             s = filter(!=('\"'), s)
#             s = split(s, " ")

#             push!(data, Int[eval(Meta.parse(a)) for a ∈ s[1:4]])
#         end
#     end

#     N = maximum([a[1] for a ∈ data])
    
#     # Construct multiplication table from associators
#     mult = zeros(Int,N,N,N)
    
#     for (_,a,b,c) ∈ filter(e -> e[1] == 1, data) 
        
#         mult[a,b,c] = 1 
#     end

#     return mult 
# end

function fusion_ring_name(m::Array{Int,3})
    r = size(m,1)

    Iᵣ = [i == j ? 1 : 0 for i ∈ 1:r, j ∈ 1:r]

    # Group simples by One, self dual, and not self dual
    one = findfirst(i -> all([m[i,:,:] == m[:,i,:] == Iᵣ]), 1:r)
    self_dual = findall(i -> i != one && m[i,i,one] == 1, 1:r)
    non_self_dual = findall(i -> i != one && m[i,i,one] == 0, 1:r)
    
    n = length(non_self_dual)

    # Compute FPdims to sort groups
    fpdims = maximum.(filter(isreal, eigenvalues(QQBarField(), matrix(QQ, r,r, m[i,:,:]))) for i ∈ 1:r)


    #Compute dual pairs
    dual_pairs = Tuple.(unique(Set.(Tuple.(findall(==(1), m[:,:,1] .- Iᵣ)))))

    # First canonical ordering
    self_dual = sort(self_dual, by =  e -> fpdims[e])
    pairs = sort(dual_pairs, by = e -> fpdims[e[1]])
    non_self_dual = vcat(dual_pairs...)

    # get all permutations fixing the ordering rules 
    self_dual_perms = 
        if !isempty(self_dual)
            Sₘ = symmetric_group(length(self_dual))
            elements(stabilizer(Sₘ, fpdims[self_dual], permuted)[1])
        else 
            [symmetric_group(1)[0]]
        end
    
    pairs_perms = 
        if !isempty(dual_pairs) 
            Sₖ = symmetric_group(length(dual_pairs))
            elements(stabilizer(Sₖ, fpdims[getindex.(dual_pairs,1)], permuted)[1])
        else
            [symmetric_group(1)[0]]
        end

    # Create the signature for every permutation 
    signatures = []

    base = maximum(m) + 1
    
    fixed_order = [one; self_dual; vcat(collect.(dual_pairs)...)]

    for p1 ∈ self_dual_perms, p2 ∈ pairs_perms 
        
        binary_choice = Base.product([fpdims[i] == fpdims[j] ? [true,false] : [false] for (i,j) ∈ dual_pairs]...)

        # Add the permutations inside the pairs where possible
        perms = [[one; permuted(self_dual, p1); vcat([rev ? collect(reverse(v)) : collect(v) for (v,rev) ∈ zip(permuted(dual_pairs, p2), bool)]...)] for bool ∈ binary_choice]

        for p3 ∈ perms
            val = ""
            for i ∈ p3, j ∈ p3, k ∈ p3 
                val = val * "$(m[i,j,k])"
            end

            push!(signatures, parse(ZZRingElem, val, base))
        end
    end
    return r,n, maximum(signatures)
end

#=----------------------------------------------------------
    Save to anyonwiki 
----------------------------------------------------------=#

function save_fusion_category(C::SixJCategory, path::String, name::String)
    cat_path = joinpath(path, name)

    mkdir(cat_path)

    save_fusion_category_meta_data(C, joinpath(cat_path, "$(name)_meta"))

    save_symbols(F_symbols(C), joinpath(cat_path, "$(name)_F_symbols"), 4)

    save_symbols(P_symbols(C), joinpath(cat_path, "$(name)_P_symbols"))
    
    if is_braided(C) 
        save_symbols(R_symbols(C), joinpath(cat_path, "$(name)_R_symbols"))
    end
end

function anyonwiki_center_meta(i,j,k,l,m,n,o)
    p = anyonwiki_center_artifact_path(i,j,k,l,m,n,o)
    name = splitpath(p)[end]

    meta = include(joinpath(p, "$(name)_meta"))
end

function load_fusion_category(file::String)
    
    name = splitpath(file)[end]

    meta = include(joinpath(file, "$(name)_meta"))

    K = meta["field"]
    rank = meta["rank"]
    description = meta["name"]
    simples_names = meta["simples_names"]
    one = meta["one"]

    # include F/P/R-symbols as coefficient vectors, convert to number field elements and then to matrices
    F_symbols = load_F_symbols(rank,K,joinpath(file, "$(name)_F_symbols"))

    P_symbols = include(joinpath(file, "$(name)_P_symbols"))
    P_symbols = [K == QQ ? K(P_symbols[k]...) : K(P_symbols[k]) for k ∈ sort(collect(keys(P_symbols)))]

    C = six_j_category(K,  multiplication_table_from_F_symbols(F_symbols))
    set_associator!(C, F_symbols)
    set_pivotal!(C, P_symbols)

    if haskey(meta, "embedding")
        r = meta["embedding"]
        if K == QQ
            setfield!(C, :embedding, complex_embedding(rationals_as_number_field()[1], r))
        else
            setfield!(C, :embedding, complex_embedding(K, r))
        end
    end
    
    if isfile(joinpath(file, "$(name)_R_symbols"))
        R_symbols = load_R_symbols(rank,K,joinpath(file, "$(name)_R_symbols"))
        set_braiding!(C, R_symbols)
    end

    set_name!(C, description)
    set_simples_names!(C, simples_names)
    set_one!(C, one)

    C
end

function anyonwiki_center_multiplication_table(i,j,k,l,m,n,o)
    p = anyonwiki_center_artifact_path(i,j,k,l,m,n,o)
    name = splitpath(p)[end]

    meta = include(joinpath(p, "$(name)_meta"))

    rank = meta["rank"]
   
    p2 = joinpath(p, "$(name)_F_symbols")

    dir = filter(e -> e[1:3] == "[1,", readdir(p2))

    multiplicities = Dict(eval(Meta.parse(q))[2:4] => length(include(joinpath(p2,q))) for q in dir)

    [Int(sqrt(get(multiplicities, [i,j,k], 0))) for i in 1:rank, j in 1:rank, k in 1:rank]
end

function anyonwiki_center_grothendieck_ring(i,j,k,l,m,n,o)
    meta = anyonwiki_center_meta(i,j,k,l,m,n,o)
    names = meta["simples_names"]
    m = anyonwiki_center_multiplication_table(i,j,k,l,m,n,o)
    ℕRing(names, m, [1; zeros(Int, length(names)-1)])
end


function load_F_symbols(rank::Int, K::Field, path::String)
    ass = Array{MatElem,4}(undef, rank,rank,rank,rank)

    for i ∈ 1:rank, j ∈ 1:rank, k ∈ 1:rank, l ∈ 1:rank 
        _file = joinpath(path, "[$(i), $(j), $(k), $l]")

        if isfile(_file)
            symbols = include(_file)
            symbols_keys = collect(keys(symbols))
            if length(first(keys(symbols))) == 6 
                symbols_keys = sort(symbols_keys, by = v -> v[[6,5]])
            else
               symbols_keys = sort(symbols_keys, by = v -> v[[8,5,10,9,7,6]])
            end
            n = Int(sqrt(length(symbols_keys)))
            vals = [K == QQ ? K(symbols[v]...) : K(symbols[v]) for v ∈ symbols_keys]
            M = matrix(K,n,n, vals)
            ass[i,j,k,l] = transpose(M)
        else
            ass[i,j,k,l] = zero_matrix(K,0,0)
        end
    end
    ass 
end

function load_R_symbols(rank::Int, K::Field, path::String)
    braid = [zero_matrix(K,0,0) for _ ∈ 1:rank, _ ∈ 1:rank, _ ∈ 1:rank]
    symbols = include(path)
    chunks = group_dict_keys_by(e -> e[1:3], symbols)

    for ((i,j,k), D) ∈ chunks 
        
        symbols_keys = sort(collect(keys(D)))
        n = Int(sqrt(length(D)))
        vals = [K == QQ ? K(D[v]...) : K(D[v]) for v ∈ symbols_keys]

        M = matrix(K,n,n, vals)
        braid[i,j,k] = transpose(M)
    end
    braid 
end



function save_symbols(S::Dict, path::String, chunk::Int = 0)
    K = parent(first(S)[2])

    if chunk != 0
        chunks = group_dict_keys_by(e -> e[1:chunk], S)
        mkdir(path)

        for (k,ch) ∈ chunks
            open(joinpath(path, "$k"), "w") do io 
                write(io, "Dict(\n")
                
                if K == QQ
                    write(io, join(["\t$k => $([v])" for (k,v) ∈ ch], ",\n") )
                else
                    write(io, join(["\t$k => $(coefficients(v))" for (k,v) ∈ ch], ",\n"))
                end

                write(io, ")")
            end
        end
    else
        open(path, "w") do io 
            write(io, "Dict(\n")
            
            if K == QQ
                write(io, join(["\t$k => $([v])" for (k,v) ∈ S], ",\n") )
            else
                write(io, join(["\t$k => $(coefficients(v))" for (k,v) ∈ S], ",\n"))
            end

            write(io, ")")
        end
    end
end

function save_fusion_category_meta_data(C::SixJCategory, file::String)
    open(file, "w") do io 
        write(io, "# Meta data for $C\n\n")
        write(io, """Dict(\n
        \t\"name\" => \"$(C.name)\",\n""")
        if base_ring(C) == QQ 
            write(io, "\t\"field\" => QQ,\n")
        else
            write(io, "\t\"field\" => number_field(polynomial(QQ,$(collect(coefficients(base_ring(C).pol)))))[1],\n")
        end
        write(io, "
        \t\"rank\"=> $(rank(C)),\n
        \t\"multiplicity\" => $(multiplicity(C)),\n
        \t\"simples_names\" => $(simples_names(C)),\n
        \t\"one\" => $(C.one)
        ")

        if isdefined(C, :embedding)
            r = getfield(C, :embedding).r
            write(io, ",\n\t\"embedding\" => AcbField()(\"$(string(real(r)))\") + AcbField()(\"$(string(imag(r)))\")*AcbField()(im)\n")
        end
        write(io, ")")
    end
end


