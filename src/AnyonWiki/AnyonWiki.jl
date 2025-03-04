#=----------------------------------------------------------
    Load fusion categories from the Anyon Wiki by 
    Gert Vercleyen. 
----------------------------------------------------------=#

associator_path = joinpath(@__DIR__, "MultFreeCenters/6jSymbols/")
pivotal_path = joinpath(@__DIR__, "MultFreeCenters/PivotalStructures/")

@doc raw""" 

    anyonewiki(n::Int)

Load the n-th fusion category from the list of multiplicity free fusion categories of rank ≤ 7.
"""
function anyonwiki(n::Int, K::Ring = load_anyon_number_field(n::Int))
    
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