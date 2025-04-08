#=----------------------------------------------------------
    Script to compute the centers of the anyon wiki 
----------------------------------------------------------=#

log = open("Centers_of_anyon_wiki.log", "a")

for n ∈ 48:200  
     
    C = anyonwiki(n, QQBarField())
    
    if ! pentagon_axiom(C) 
        write(log, "Pentagon axiom failed for n = $n\n")
        println("Pentagon axiom failed for n = $n")
        continue
    end

    if ! is_spherical(C)
        write(log, "Spherical axiom failed for n = $n\n")
        println("Spherical axiom failed for n = $n")
        continue
    end
    
    Z = center(C) 

    skel_Z = six_j_category(Z)

    add_to_local_database(skel_Z, joinpath(@__DIR__, "MultFreeCenters/Centers/Center_$n"))
end