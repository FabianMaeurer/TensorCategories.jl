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
function anyonwiki(n::Int, K::Ring = load_anyon_number_field(n::Int); cyclotomic = false)
    
    cyclotomic && return anyonwiki_cyclotomic(n)

    data = []
    CC = AcbField(512)
    Q = QQBarField()
    deg = typeof(K) <: Union{<:NumField,QQField} ? degree(K) : 24
    if typeof(K) <: Union{<:NumField,QQField}
        es = complex_embeddings(K) 
        r = load_anyon_associator_root(n)
        r2 = load_anyon_pivotal_root(n) 
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
        set_pivotal!(C, load_anyon_pivotal(n, K, e))
    else 
        set_pivotal!(C, load_anyon_pivotal(n, K, nothing))     
    end
    

    C
end

function load_anyon_number_field(n::Int)
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

function load_anyon_associator_root(n::Int)
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

function load_anyon_pivotal_root(n::Int)
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

function load_anyon_pivotal(n::Int, K = load_anyon_number_field(n), e = complex_embeddings(K)[1])

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


function load_anyon_attributes(n::Int) 
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
    cyclopos = [parse(Int, s) for s ∈ readlines(joinpath(anyon_path, "cycloPos.txt"))]
    attrs = load_anyon_attributes(n)

    if n ∉ cyclopos 
        error("Category number $n is not cyclotomic")
    end

    f_symbol_string = readlines(joinpath(anyon_path, "cyclic_6j_Symbols/cat_$n.jl"))[4]
    
    r_symbol_string = if attrs[:braided]
        readlines(joinpath(anyon_path, "cyclic_R_Symbols/cat_$n.jl"))[4]
    else 
        ""
    end

    # read root
    i = findfirst(r"z\d*", f_symbol_string)
    j = if attrs[:braided]
        findfirst(r"z\d*", r_symbol_string)
    else 
        nothing
    end 

    if i === nothing && j === nothing
        global K = QQ 
    else
        roots = [parse(Int, s[k][2:end]) for (s,k) ∈ zip([f_symbol_string, r_symbol_string], [i,j]) if k !== nothing]

        r = lcm([1 ; roots]...)

        eval(Meta.parse("global K, z$r = cyclotomic_field($r)"))
    end

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

    # Build multiplication_table
    mult = zeros(Int,N,N,N)
    
    for (_,a,b,c) ∈ filter(e -> e[1] == 1, keys(associator_dict)) 
        
        mult[a,b,c] = 1 
    end

    C = six_j_category(K, mult)
    set_name!(C, "Fusion category number $n from AnyonWiki")
    set_associator!(C, associator_array)
    set_one!(C, [i == 1 for i ∈ 1:N])

    pivotals = include(joinpath(anyon_path, "cyclic_P_Symbols/cat_$n.jl"))

    set_pivotal!(C, K.([pivotals[[i]] for i ∈ 1:length(pivotals)]))


    if attrs[:braided] 
        braiding_dict = include(joinpath(anyon_path, "cyclic_R_symbols/cat_$n.jl"))
        braiding_array = Array{MatElem}(undef, N,N,N)

        for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N
            ab_c = filter(e -> e[[1,2,3]] == [a,b,c], collect(keys(braiding_dict)))
            l = Int(sqrt(length(ab_c)))
     
            M = matrix(K,l,l, [braiding_dict[v] for v ∈ ab_c])
            braiding_array[a,b,c] = M
        end

        set_braiding!(C, braiding_array)
    end

    return C
    
end
