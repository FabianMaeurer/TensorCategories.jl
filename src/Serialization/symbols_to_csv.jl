#=----------------------------------------------------------
    Save F- and R-symbols as a readable csv file
----------------------------------------------------------=#

function numeric_symbols_to_csv(destination::String, F::Dict{Vector{Int}, AcbFieldElem}; delimiter = ", ")
    open(destination, "w") do io 
        sorted = collect(sort(F))

        for (k,v) ∈ sorted 
            s = join([string.(k); string(BigComplex(v))], delimiter) * "\n"
            write(io, s)
        end
    end
end

function numeric_symbols_from_csv(file::String, K::Field = AcbField(64); delimiter = ", ")
    F = Dict{Vector{Int}, elem_type(K)}()

    lines = readlines(file) 

    # check wether complex numbers are given as one or two columns
    split_number = length(split(first(lines), delimiter)) ∈ [8,12]
    
    for l ∈ lines 
        chunks = split(l, delimiter)
        index = parse.(Int, chunks[1:end-1 - split_number])

        real,imag = if split_number 
            K.(chunks[end-1:end])
        else
            K.(split(chunks[end][1:end-2], "+"))
        end
        
        push!(F, index => real + imag*K(im))
    end
    return F 
end

function load_numeric_fusion_category(F::String, K::AcbField = AcbField(64); delimiter = ", ", transpose = false)
    F = numeric_symbols_from_csv(F, K, delimiter = delimiter)

    ass = dict_to_associator(F)

    mult = multiplication_table_from_F_symbols(ass)
    N, = size(mult)

    one = findfirst([mult[i,:,:] == [Int(j == k) for k ∈ 1:N, j ∈ 1:N] for i ∈ 1:N])

    C = six_j_category(K, mult)

    transpose && (ass = Oscar.transpose.(ass))

    set_associator!(C, ass)

    set_one!(C, [Int(i == 1) for i in 1:N])


    set_canonical_spherical!(C)

    C
end

function load_numeric_fusion_category(F::String, R::String, K::AcbField = AcbField(64), delimiter = ", ", transpose = false)
    F = numeric_symbols_from_csv(F, K, delimiter = delimiter)
    R = numeric_symbols_from_csv(R, K, delimiter = delimiter)

    ass = dict_to_associator(F)
    braid = dict_to_braiding(R)

    mult = multiplication_table_from_F_symbols(ass)
    N, = size(mult)

    one = findfirst([mult[i,:,:] == [Int(j == k) for k ∈ 1:N, j ∈ 1:N] for i ∈ 1:N])

    C = six_j_category(K, mult)

    set_associator!(C, ass)
    set_braiding!(C, braid)

    set_one!(C, [Int(i == 1) for i in 1:N])

    set_canonical_spherical!(C)

    C
end