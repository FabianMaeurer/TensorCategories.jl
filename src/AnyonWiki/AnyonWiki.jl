#=----------------------------------------------------------
    Load fusion categories from the Anyon Wiki by 
    Gert Vercleyen. 
----------------------------------------------------------=#

associator_path = joinpath(@__DIR__, "AnyonWikiData/6jSymbols/")
pivotal_path = joinpath(@__DIR__, "AnyonWikiData/PivotalStructures/")
anyon_path = joinpath(@__DIR__, "AnyonWikiData/")

@doc raw""" 

    anyonewiki(n::Int)

Load the n-th fusion category from the list of multiplicity free fusion categories of rank ≤ 7.
"""
function anyonwiki(n::Int, K::Union{Ring,Nothing} = nothing ; minimal = false)
    
    if ! minimal && K === nothing 
        if n ∈ ZZ.(readlines(joinpath(anyon_path, "cycloPos.txt")))
            return anyonwiki_cyclotomic(n)
        end
    end

    if K === nothing && minimal  
        K = load_anyonwiki_number_field(n::Int)
    elseif K === nothing
        K = QQBarField()
    end        


    data = []
    CC = AcbField(512)
    Q = QQBarField()
    deg = typeof(K) <: Union{<:NumField,QQField} ? degree(K) : 24
    if typeof(K) <: Union{<:NumField,QQField}
        es = complex_embeddings(K) 
        r = load_anyonwiki_associator_root(n)
        r2 = load_anyonwiki_pivotal_root(n) 
        for f ∈ es 
            try 
                preimage(f,r)
                preimage(f,r2)
                global e = f
                break 
            catch 
            end
        end
    end

    # Open File n and extract data
    open(joinpath(associator_path, "cat_$n")) do f 

        while ! eof(f)

            s = readline(f) 
            s = filter(!=('\"'), s)
            s = split(s, " ")

            indices = Int[eval(Meta.parse(a)) for a ∈ s[1:end-3]]

            found = false 

            while ! found
                x_real = CC(s[end-2] * "+/- 1e-510")
                x_imag = CC(s[end][1:end-2] * "+/- 1e-510") * CC(im) 
                try 
                global x = guess(Q, x_real + x_imag, deg)
                found = true
                catch  err
                    if precision(CC) > 5000 error(err) end
                    CC = AcbField(2*precision(CC)) 
                end
            end
            if typeof(K) <: Union{<:NumField,QQField}
                append!(data, [Any[indices; preimage(e,x)]])
            else 
                append!(data, [Any[indices; K(x)]])       
            end
        end
    end
    N = maximum([a[1] for a ∈ data])
    
    # Construct multiplication table from associators
    mult = zeros(Int,N,N,N)
    
    for (_,a,b,c) ∈ filter(e -> e[1] == 1, data) 
        
        mult[a,b,c] = 1 
    end

    C = six_j_category(K, mult)
    set_name!(C, "Fusion category number $n from AnyonWiki List")

    # Build matrices for 6j-Symbols
    for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N, d ∈ 1:N
        abc_d = filter(e -> e[[1,2,3,4]] == [a,b,c,d], data)
        l = Int(sqrt(length(abc_d)))
        abc_d = sort(abc_d, by = v -> v[6])
        abc_d = sort(abc_d, by = v -> v[5])
        
        M = matrix(K,l,l, [v[7] for v ∈ abc_d])
        set_associator!(C,a,b,c,d,M)  
    end
    set_one!(C, 1)
    
    if typeof(K) <: Union{<:NumField,QQField}
        set_pivotal!(C, load_anyonwiki_pivotal(n, K, e))
    else 
        set_pivotal!(C, load_anyonwiki_pivotal(n, K, nothing))     
    end
    

    C
end

function load_anyonwiki_number_field(n::Int)
    line = 1
    open(joinpath(associator_path, "Polynomials")) do f 

        while ! eof(f) 

            s = readline(f) 
            
            if line == n 
                global _,x = QQ[:x]
                min_pol = eval(Meta.parse(s))

                global K1 = min_pol == 1 ? QQ : number_field(min_pol)[1]
                break
            end
            
            line += 1
        end
    end
    
    line = 1
    open(joinpath(pivotal_path, "Polynomials")) do f 

        while ! eof(f) 

            s = readline(f) 
            
            if line == n 
                global _,x = QQ[:x]
                min_pol = eval(Meta.parse(s))

                global K2 =  min_pol == 1 ? QQ : number_field(min_pol)[1]
                break
            end
            
            line += 1
        end
    end
    
    K1 == K2 && return K1
    K1 == QQ && return K2
    K2 == QQ && return K1 

    simplify(splitting_field([defining_polynomial(K1), defining_polynomial(K2)]))[1]
end

function load_anyonwiki_associator_root(n::Int)
    line = 1
    CC = AcbField(1024)
    Q = QQBarField()
    deg = 24
    open(joinpath(associator_path, "Roots.dat")) do f 

        while ! eof(f) 
            s = readline(f) 
            if line == n
                s = filter(!=('\"'), s)
                s = split(s, " ")
                x_real = CC(s[1] * " +/- 1e-512")
                x_imag = CC(s[end][1:end-2] * "+/- 1e-512") * CC(im) 

                return guess(Q,x_real + x_imag,deg)
            end
            line += 1
        end
    end
end

function load_anyonwiki_pivotal_root(n::Int)
    line = 1
    CC = AcbField(1024)
    Q = QQBarField()
    deg = 24
    open(joinpath(pivotal_path, "Roots.dat")) do f 

        while ! eof(f) 
            s = readline(f) 
            if line == n
                s = filter(!=('\"'), s)
                s = split(s, " ")
                x_real = CC(s[1] * " +/- 1e-512")
                x_imag = CC(s[end][1:end-2] * "+/- 1e-512") * CC(im) 

                return guess(Q,x_real + x_imag,deg)
            end
            line += 1
        end
    end
end

function load_anyonwiki_pivotal(n::Int, K = load_anyonwiki_number_field(n), e = complex_embeddings(K)[1])

    CC = AcbField(512)
    Q = QQBarField()
    deg = typeof(K) <: Union{<:NumField,QQField} ? degree(K) : 24

    piv = elem_type(K)[]
    open(joinpath(pivotal_path, "pivots_cat_$n")) do f 

        while ! eof(f) 

            s = readline(f) 
            s = filter(!=('\"'), s)
            s = split(s, " ")

            x_real = CC(s[2] * " +/- 1e-512")
            x_imag = CC(s[4][1:end-2] * "+/- 1e-512") * CC(im) 

            x = guess(Q, x_real + x_imag, deg)

            if typeof(K) <: Union{<:NumField,QQField}
                push!(piv, preimage(e, x))
            else 
                push!(piv, K(x))    
            end
        end
    end
    return piv 
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

function anyonwiki_cyclotomic(n::Int, K::Field)
    
    r = anyonwiki_cyclotomic_root(n)
    root = root_of_unity(K,r)
    eval(Meta.parse("z$r = $(root)"))


    associator_array = anyonwiki_cyclotomic_associator(n,K)
   
    mult = multiplication_table_from_F_symbols(associator_array)


    C = six_j_category(K, mult)
    set_name!(C, "Fusion category number $n from AnyonWiki")
    set_associator!(C, associator_array)
    set_one!(C, [i == 1 for i ∈ 1:C.simples])

    pivotals = include(joinpath(anyon_path, "cyclic_P_Symbols/cat_$n.jl"))

    set_pivotal!(C, K.([pivotals[[i]] for i ∈ 1:length(pivotals)]))

    attrs = load_anyonwiki_attributes(n)

    if attrs[:braided] 
       
        braiding_array = anyonwiki_cyclotomic_braiding(n, K)
        set_braiding!(C, braiding_array)
    end

    return C
    
end

function anyonwiki_finite(n::Int)
    r = anyonwiki_cyclotomic_root(n)

    f_symbol_string = readlines(joinpath(anyon_path, "cyclic_6j_Symbols/cat_$n.jl"))[4]

    ind = findall(r"\//\d+", f_symbol_string)
    inverse = lcm([parse(Int, f_symbol_string[i[3:end]]) for i in ind])

    K = finite_prime_field_with_root_of_unity(r, inverse + 1)

    anyonwiki_cyclotomic(n,K)
end

function multiplication_table_from_F_symbols(ass::Array{MatElem,4})
    # Build multiplication_table
    N, = size(ass)

    mult = zeros(Int,N,N,N)
    
    for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N 
        mult[i,j,k] = size(ass[1,i,j,k])[1]
    end
    return mult
end

function anyonwiki_cyclotomic_braiding(n::Int, K::Field)
    r = anyonwiki_cyclotomic_root(n)
    root = root_of_unity(K,r)
    eval(Meta.parse("z$r = $(root)"))


    braiding_dict = include(joinpath(anyon_path, "cyclic_R_symbols/cat_$n.jl"))

    N = maximum([a[1] for a ∈ keys(braiding_dict)])

    braiding_array = Array{MatElem}(undef, N,N,N)

    for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N
        ab_c = filter(e -> e[[1,2,3]] == [a,b,c], collect(keys(braiding_dict)))
        l = Int(sqrt(length(ab_c)))

        M = matrix(K,l,l, [braiding_dict[v] for v ∈ ab_c])
        braiding_array[a,b,c] = M
    end

    return braiding_array
end

function anyonwiki_cyclotomic_root(n::Int)
    cyclopos = [parse(Int, s) for s ∈ readlines(joinpath(anyon_path, "cycloPos.txt"))]

    attrs = load_anyonwiki_attributes(n)

    if n ∉ cyclopos 
        error("Category number $n is not cyclotomic")
    end

    f_symbol_string = readlines(joinpath(anyon_path, "cyclic_6j_Symbols/cat_$n.jl"))[4]
    
    p_symbol_string = readlines(joinpath(anyon_path, "cyclic_P_Symbols/cat_$n.jl"))[4]

    r_symbol_string = if attrs[:braided]
        readlines(joinpath(anyon_path, "cyclic_R_Symbols/cat_$n.jl"))[4]
    else 
        ""
    end

    # read root
    i = findfirst(r"z\d*", f_symbol_string)
    j = findfirst(r"z\d*", p_symbol_string)
    k = if attrs[:braided]
        findfirst(r"z\d*", r_symbol_string)
    else 
        nothing
    end 

    if i === nothing && j === nothing && k == nothing
        return 1 
    else
        roots = [l === nothing ? 1 : parse(Int, s[l][2:end]) for (s,l) ∈ zip([f_symbol_string,p_symbol_string, r_symbol_string], [i,j,k])]

        return r = lcm(roots...)
    end
end

function anyonwiki_cyclotomic_associator(n::Int, K::Field)
    r = anyonwiki_cyclotomic_root(n)
        root = root_of_unity(K,r)
    eval(Meta.parse("z$r = $(root)"))

    # read 6j-symbols 
    associator_dict = include(joinpath(anyon_path, "cyclic_6j_Symbols/cat_$n.jl"))

    N = maximum([a[1] for a ∈ keys(associator_dict)])

    associator_array = Array{MatElem}(undef, N,N,N,N)
    
    # Build matrices for 6j-Symbols
    for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N, d ∈ 1:N
        abc_d = filter(e -> e[[1,2,3,4]] == [a,b,c,d], collect(keys(associator_dict)))
        l = Int(sqrt(length(abc_d)))
        abc_d = sort(abc_d, by = v -> v[6])
        abc_d = sort(abc_d, by = v -> v[5])
        
        M = matrix(K,l,l, [associator_dict[v] for v ∈ abc_d])
        associator_array[a,b,c,d] = M
    end
    return associator_array 
end

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
function load_anyonwiki_fusion_rules(n::Int)

    data = []

    open(joinpath(associator_path, "cat_$n")) do f 

        while ! eof(f)

            s = readline(f) 
            s = filter(!=('\"'), s)
            s = split(s, " ")

            push!(data, Int[eval(Meta.parse(a)) for a ∈ s[1:4]])
        end
    end

    N = maximum([a[1] for a ∈ data])
    
    # Construct multiplication table from associators
    mult = zeros(Int,N,N,N)
    
    for (_,a,b,c) ∈ filter(e -> e[1] == 1, data) 
        
        mult[a,b,c] = 1 
    end

    return mult 
end

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
            @show val
            push!(signatures, parse(ZZRingElem, val, base))
        end
    end
    return r,n, maximum(signatures)
end