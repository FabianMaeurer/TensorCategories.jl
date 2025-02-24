#=----------------------------------------------------------
    Load fusion categories from the Anyon Wiki by 
    Gert Vercleyen. 
----------------------------------------------------------=#

associator_path = joinpath(@__DIR__, "MultFreeCenters/6jSymbols/")
pivotal_path = joinpath(@__DIR__, "MultFreeCenters/PivotalStructures/")

@doc raw""" 

    anyonewiki(n::Int64)

Load the n-th fusion category from the list of multiplicity free fusion categories of rank ≤ 7.
"""
function anyonwiki(n::Int)
    
    data = []
    CC = AcbField(512)
    Q = QQBarField()
    K = load_anyon_number_field(n::Int)
    deg = degree(K)
    e = complex_embeddings(K)[1]

    # Open File n and extract data
    open(joinpath(associator_path, "cat_$n")) do f 

        while ! eof(f)

            s = readline(f) 
            s = filter(!=('\"'), s)
            s = split(s, " ")

            indices = Int[eval(Meta.parse(a)) for a ∈ s[1:end-3]]
            x_real = CC(s[end-2] * " +/- 1e-512")
            x_imag = CC(s[end][1:end-2] * "+/- 1e-512") * CC(im) 

            x = guess(Q, x_real + x_imag, deg)

            append!(data, [Any[indices; preimage(e,x)]])
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
        M = matrix(K,l,l, [v[7] for v ∈ abc_d])
        set_associator!(C,a,b,c,d,M)  
    end
    set_one!(C, 1)
    
    set_pivotal!(C, load_anyon_pivotal(n, K, e))

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

function load_anyon_pivotal(n::Int, K = load_anyon_number_field(n), e = complex_embeddings(K)[1])

    CC = AcbField(512)
    Q = QQBarField()
    deg = degree(K) 

    piv = elem_type(K)[]
    open(joinpath(pivotal_path, "pivots_cat_$n")) do f 

        while ! eof(f) 

            s = readline(f) 
            s = filter(!=('\"'), s)
            s = split(s, " ")

            x_real = CC(s[2] * " +/- 1e-512")
            x_imag = CC(s[4][1:end-2] * "+/- 1e-512") * CC(im) 

            x = guess(Q, x_real + x_imag, deg)

            append!(piv, [preimage(e, x)])
        end
    end
    return piv 
end