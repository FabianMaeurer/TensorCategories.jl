#=----------------------------------------------------------
    Script to compute the centers of the anyon wiki 
----------------------------------------------------------=#

log = open(joinpath(@__DIR__, "Cyclotomic Centers.log"), "a")
time_log = open(joinpath(@__DIR__, "Cyclotomic Centers Time.log"), "a")

for n ∈ 86:200
     
    attributes = TensorCategories.load_anyon_attributes(n)

    # Works only for spherical categories
    if !attributes[:spherical] continue end

    C = TensorCategories.anyonwiki_cyclotomic(n)
    
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

        if attributes[:modular]
            global t1 = @elapsed begin 
                Z = C ⊠ reverse_braiding(C)
            end
            set_name!(Z, "Center of fusion category number $n of AnyonWiki")

            println("Skeletal Center computed in $t1 seconds")
            add_to_local_database(Z, joinpath(@__DIR__, "AnyonWikiData/CyclotomicCenters/CyclotomicCenter_$n"))

            write(time_log, "$n, $t1\n")
            flush(time_log)
            continue
        end


        t1 = @elapsed begin 
            Z = center(C) 
            simples(Z)
            Z = split_cyclotomic(Z)
            sort_simples_by_dimension!(Z)
        end

        print("Center computed in $t1 seconds, ")

        t2 = @elapsed skel_Z = six_j_category(Z)

        println("Skeleton computed in $t2 seconds")

        write(time_log, "$n, $t1, $t2\n")
        flush(time_log)

        add_to_local_database(skel_Z, joinpath(@__DIR__, "AnyonWikiData/CyclotomicCenters/CyclotomicCenter_$n"))

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