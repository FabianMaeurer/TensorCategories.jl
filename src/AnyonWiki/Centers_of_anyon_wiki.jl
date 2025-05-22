#=----------------------------------------------------------
    Script to compute the centers of the anyon wiki 
----------------------------------------------------------=#

log = open(joinpath(@__DIR__, "Logs/Centers_of_anyon_wiki.log"), "a")
time_log = open(joinpath(@__DIR__, "Logs/Centers_of_anyon_wiki_time.log"), "a")

for n ∈ 201:300  
     
    attributtes = TensorCategories.load_anyonwiki_attributes(n)

    # Works only for spherical categories
    if attributtes[3] != 1 continue end

    # Skip modular ones for now
    if attributtes[5] != 0 continue end
    
    C = anyonwiki(n, QQBarField())
    
    # if ! pentagon_axiom(C) 
    #     write(log, "Pentagon axiom failed for n = $n\n")
    #     println("Pentagon axiom failed for n = $n")
    #     continue
    # end

    # if ! is_spherical(C)
    #     write(log, "Spherical axiom failed for n = $n\n")
    #     println("Spherical axiom failed for n = $n")
    #     continue
    # end
    
    print("n = $n: ")

    try 
        t1 = @elapsed begin 
            Z = center(C) 
            simples(Z)
        end

        print("Center computed in $t1 seconds, ")

        t2 = @elapsed skel_Z = six_j_category(Z)

        println("Skeleton computed in $t2 seconds")

        write(time_log, "$n, $t1, $t2\n")
        flush(time_log)

        add_to_local_database(skel_Z, joinpath(@__DIR__, "AnyonWikiData/Centers/Center_$n"))

    catch e 
        println(" failed")
        
        str = IOBuffer()
        showerror(str, e) 
        str = String(take!(str))
        write(log, """n = $n: 
        $(str)\n\n""")
        flush(log)
    end
end