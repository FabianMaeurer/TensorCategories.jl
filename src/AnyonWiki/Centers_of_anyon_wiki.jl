#=----------------------------------------------------------
    Script to compute the centers of the anyon wiki 
----------------------------------------------------------=#
using Revise, TensorCategories, Oscar

log = open(joinpath(@__DIR__, "Logs/Centers_of_anyon_wiki.log"), "a")
time_log = open(joinpath(@__DIR__, "Logs/Centers_of_anyon_wiki_time.log"), "a")

code_to_name(a,b,c,d,e,f,g) = "center_$(a)_$(b)_$(c)_$(d)_$(e)_$(f)_$g"

codes = sort(collect(keys(include(joinpath(@__DIR__, "base_field_generators.jl")))))

for c ∈ filter(e -> e[1] == 6, codes)[1:5:end]
    file_name = code_to_name(c...)

    if isdir(joinpath(@__DIR__, "AnyonWikiData/Centers/$file_name"))
        continue 
    end

    try 

        C = anyonwiki(c...)

        t1 = @elapsed Z = center(C)

        t2 = @elapsed Z,em = split(Z)
        
        t3 = @elapsed Z2 = six_j_category(Z)

        K = base_ring(C)
        L = base_ring(Z2)

        embs = embeddings(L)

        j = argmin([abs(e(em(gen(K))) - getfield(C, :embedding)(gen(K))) for e ∈ embs])

        setfield!(Z2, :embedding, embs[j])

        t4 = @elapsed save_fusion_category(Z2, joinpath(@__DIR__, "AnyonWikiData/Centers"), file_name)


        write(time_log, "$file_name: $t1, $t2, $t3, $t4\n")
        flush(time_log)
    catch e 
        io = IOBuffer() 
        showerror(io,e)
        msg = String(take!(io))
        write(log, "$file_name:\n$(msg)\n\n")
        flush(log)
    end
end
