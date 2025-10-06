function skiplines(io::IO, n)
    i = 1
    while i <= n
        eof(io) && error("File contains less than $n lines")
        i += read(io, Char) === '\n'
    end
end

function insertline(file, string, no)
    f = open(file, "r+");
    skiplines(f, no - 1);
    mark(f)
    buf = IOBuffer()
    write(buf, f)
    seekstart(buf)
    reset(f)
    println(f, string)
    write(f, buf)
    close(f)
end

function anyonwiki_center_simple_name_to_vec(s::String, simpls::Vector{String})
    s = split(s, ",")[1]
    s = replace(s, "(" => "", "⊕" => "", ")" => "")

    ret = []

    for S ∈ simpls 
        if occursin(Regex("(\\d+)⋅$S"), s)
            m = collect(eachmatch(Regex("(\\d+)⋅$S"),s))[1]
            push!(ret, parse(Int, match(r"\d+", m.match).match))
        elseif occursin(S, s)
            push!(ret, 1)
        else
            push!(ret, 0)
        end
    end
    ret
end

function anyonwiki_center_forgetful(i,j,k,l,m,n,o)
    names = TensorCategories.anyonwiki_center_meta(i,j,k,l,m,n,o)["simples_names"]
    vecs = transpose(hcat([anyonwiki_center_simple_name_to_vec(s, ["𝟙"; ["X$q" for q ∈ 2:i]]) for s ∈ names]...))
end

k = anyonwiki_keys(5)

p = "/tmpbig/maeurer/Data/TensorCategoriesDatabase/AnyonWikiCenters"

for i ∈ k 
    name = "center_$(i[1])_$(i[2])_$(i[3])_$(i[4])_$(i[5])_$(i[6])_$(i[7])"

    meta = include(joinpath(p, "$name", "$(name)_meta")) 

    if haskey(meta, "embedding")
        continue
    end
    println(i)
    # C = anyonwiki_center(i...)
    # Z = split(center(anyonwiki(i...)))
    C = anyonwiki(i...)
   
    L = TensorCategories.anyonwiki_center_meta(i...)["field"]
    K,_ = TensorCategories.load_anyonwiki_number_field(i...)

    f = K == QQ ? L : is_subfield(K, L)[2]

    e = getfield(C, :embedding)
    embL = complex_embeddings(L)
    j = argmin([abs(e(gen(K)) - q(f(gen(K)))) for q in embL]) 
    r = AcbField()(embL[j](gen(L)))

    insertline(joinpath(p, name, "$(name)_meta"), "\t\"embedding\" => AcbField()(\"$(string(real(r)))\") + AcbField()(\"$(string(imag(r)))\")*AcbField()(im),", 9)
    # sp = Dict(i => f(dim(S[i])) * inv(dim(Z[i])) for i ∈ 1:length(S))
   
    # TensorCategories.save_symbols(sp, joinpath(p, name, "$(name)_P_symbols"))

    # m = multiplicity(C)

    # insertline(joinpath(p, name, "$(name)_meta"), "\t\"multiplicity\" => $m,", 8)

    s = read(joinpath(p, "$name", "$(name)_meta"), String)
    s = replace(s, r"Skeletization" => "Skeletonization")
    write(joinpath(p, "$name", "$(name)_meta"), s)
    
    println("Done $name")
end