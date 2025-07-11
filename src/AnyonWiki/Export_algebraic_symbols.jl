CC = AcbField(512)

dict = sort(include(joinpath(@__DIR__, "AnyonWikiData/base_field_generators.jl")))

reg = r"\"[^\"]*\""


for (k, data) ∈ dict
    println(k)
    cat_name = "cat_$(k[1])_$(k[2])_$(k[3])_$(k[4])_$(k[5])_$(k[6])_$(k[7])"
    file = joinpath(@__DIR__, "AnyonWikiData/new_6j_symbols/$cat_name.jl")
    if k[6] != 0 file_R = joinpath(@__DIR__, "AnyonWikiData/new_R_symbols/$cat_name.jl") end
    file_P = joinpath(@__DIR__, "AnyonWikiData/new_P_symbols/$cat_name.jl")

    str = read(file, String)
    if k[6] != 0 str_R = read(file_R, String) end
    str_P = read(file_P, String)


    matches = unique([m.match for m in eachmatch(reg, k[6] == 0 ? str*str_P : str*str_R*str_P)])

    K,e = TensorCategories.load_anyonwiki_number_field(k...)

    subs =[String(m) => coefficients(preimage(e, TensorCategories.string_to_acb(CC, String(m)[2:end-1]))) for m in matches]

    str = replace(str, subs...)
    if k[6] != 0 str_R = replace(str, subs...) end
    str_P = replace(str, subs...)

    write(joinpath(@__DIR__, "AnyonWikiData/algebraic_F_symbols/$cat_name.jl"), str)
    if k[6] != 0 write(joinpath(@__DIR__, "AnyonWikiData/algebraic_R_symbols/$cat_name.jl"), str) end
    write(joinpath(@__DIR__, "AnyonWikiData/algebraic_P_symbols/$cat_name.jl"), str)


end