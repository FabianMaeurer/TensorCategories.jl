CC = AcbField(2048)

dict = sort(collect(include(joinpath(@__DIR__, "AnyonWikiData/base_field_generators.jl"))), by = x -> x[1])

reg = r"\"[^\"]*\""

error_log = open(joinpath(@__DIR__, "export.log"), "a")

Threads.@threads for (k, data) ∈ dict[936:end]
    
    cat_name = "cat_$(k[1])_$(k[2])_$(k[3])_$(k[4])_$(k[5])_$(k[6])_$(k[7])"
    println(k)
    # if isfile(joinpath(@__DIR__, "AnyonWikiData/algebraic_F_symbols/$cat_name.jl")) 
        
    #     # if k[6] != 0 file_R = joinpath(@__DIR__, "AnyonWikiData/new_R_Symbols/$cat_name.jl") end
    #     # file_P = joinpath(@__DIR__, "AnyonWikiData/new_P_Symbols/$cat_name.jl")

       
    #     # if k[6] != 0 str_R = read(file_R, String) end
    #     # str_P = read(file_P, String)


    #     # matches = unique([m.match for m in eachmatch(reg, k[6] == 0 ? str_P : str_R*str_P)])

    #     # K,e = TensorCategories.load_anyonwiki_number_field(k...)

    #     # subs = [String(m) => coefficients(preimage(e, TensorCategories.string_to_acb(CC, String(m)[2:end-1]), degree(K))) for m in matches]

    #     # if k[6] != 0 str_R = replace(str_R, subs...) end
    #     # str_P = replace(str_P, subs...)

    #     # if k[6] != 0 write(joinpath(@__DIR__, "AnyonWikiData/algebraic_R_symbols/$cat_name.jl"), str_R) end
    #     # write(joinpath(@__DIR__, "AnyonWikiData/algebraic_P_symbols/$cat_name.jl"), str_P)
    #     continue 
    # end
  

    file = joinpath(@__DIR__, "AnyonWikiData/new_6j_Symbols/$cat_name.jl")
    if k[6] != 0 file_R = joinpath(@__DIR__, "AnyonWikiData/new_R_Symbols/$cat_name.jl") end
    file_P = joinpath(@__DIR__, "AnyonWikiData/new_P_Symbols/$cat_name.jl")

    str = read(file, String)
    if k[6] != 0 str_R = read(file_R, String) end
    str_P = read(file_P, String)


    matches = unique([m.match for m in eachmatch(reg, k[6] == 0 ? str*str_P : str*str_R*str_P)])

    length(matches)
    K,e = TensorCategories.load_anyonwiki_number_field(k...)

    subs =[String(m) => coefficients(preimage(e, TensorCategories.string_to_acb(CC, String(m)[2:end-1]), degree(K))) for m in matches]

    str = replace(str, subs...)
    if k[6] != 0 str_R = replace(str_R, subs...) end
    str_P = replace(str_P, subs...)

    write(joinpath(@__DIR__, "AnyonWikiData/algebraic_F_symbols/$cat_name.jl"), str)
    if k[6] != 0 write(joinpath(@__DIR__, "AnyonWikiData/algebraic_R_symbols/$cat_name.jl"), str_R) end
    write(joinpath(@__DIR__, "AnyonWikiData/algebraic_P_symbols/$cat_name.jl"), str_P)

    C = anyonwiki(k...)
    write(error_log, "$k: $(pentagon_axiom(C))\n")
    flush(error_log)

end