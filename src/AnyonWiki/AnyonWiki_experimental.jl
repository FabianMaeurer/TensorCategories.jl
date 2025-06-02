@doc raw""" 

    anyonwiki(r::Int, n::Int, m::Int)
    anyonwiki(r::Int, n::Int, m::Int)
    anyonwiki(r::Int, n::Int, m::Int; associator::Int, braiding::Int, pivotal::Int)
    anyonwiki(r::Int, n::Int, m::Int; associator::Int, braiding::Int, pivotal::Int)

Load a categorification of the ring FR_{r,n}^i with corresponding associator, braiding and pivotal structure from the list.
"""
function anyonwiki(r::Int, n::Int, m::Int; associator::Int = 1, braiding::Int = 1, pivotal::Int = 1)

    # The corresponding directory
    dir = joinpath(@__DIR__, "NumericCategories/FR_$(r)_1_$(n)_$m")

    ! isdir(dir) && error("There is no rank $(m)th $r fusion category with $n non self dual simples")

    ! isfile(joinpath(dir, "pentsol_$(associator).txt")) && error("There is no $(associator)th associator")

    data = []

    # Open the corresponding file and extract data
    open(joinpath(dir, "pentsol_$(associator).txt")) do f 

        CC = AcbField(256)
        RR = ArbField(256)
        while ! eof(f)

            s = readline(f) 
            #s = filter(!=(''), s)
            s = split(s, "\t")

            indices = Int[eval(Meta.parse(a)) for a ∈ s[1:end-2]]

            found = false 

            while ! found
                x_real = RR(s[end-1] * "+/- 2e-128")
                x_imag = s[end][1:end] == "0" ? RR(0) : CC(s[end][1:end] * "+/- 2e-128") * CC(im) 
                try 
                global x = guess(QQBarField(), x_real + x_imag, 24)
                found = true
                catch  err
                    return @show x_real + x_imag
                    # if precision(CC) > 5000 error(err) end
                    # CC = AcbField(2*precision(CC)) 
                end
            end
            append!(data, [Any[indices; x]]) 
        end
    end
    N = maximum([a[1] for a ∈ data])
    
    # Construct multiplication table from associators
    mult = zeros(Int,N,N,N)
    
    for (_,a,b,c) ∈ filter(e -> e[1] == 1, data) 
        
        mult[a,b,c] = 1 
    end

    C = six_j_category(QQBarField(), mult)
    set_name!(C, "Fusion category of rank $r with $n non self dual objects")

    # Build matrices for 6j-Symbols
    for a ∈ 1:N, b ∈ 1:N, c ∈ 1:N, d ∈ 1:N
        abc_d = filter(e -> e[[1,2,3,4]] == [a,b,c,d], data)
        l = Int(sqrt(length(abc_d)))
        abc_d = sort(abc_d, by = v -> v[6])
        abc_d = sort(abc_d, by = v -> v[5])
        
        M = matrix(QQBarField(),l,l, [v[7] for v ∈ abc_d])
        set_associator!(C,a,b,c,d,M)  
    end
    set_one!(C, 1)
    
    # if typeof(K) <: Union{<:NumField,QQField}
    #     set_pivotal!(C, load_anyonwiki_pivotal(n, K, e))
    # else 
    #     set_pivotal!(C, load_anyonwiki_pivotal(n, K, nothing))     
    # end
    

    C
end
