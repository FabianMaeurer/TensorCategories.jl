#=----------------------------------------------------------
    Access Fusion Category Database  
----------------------------------------------------------=#

function load_fusion_category(name::String)
    load(joinpath(@__DIR__, name * ".mrdi"))
end

function add_to_local_database(C::SixJCategory, name::String)
    files = readdir(joinpath(@__DIR__))

    if name * ".mrdi" ∈ files 
        println("File already existent. Do you want to overwrite? [y/n]")
        ans = ""
        while ans ∉ ["y","n"]
            ans = readline()
        end

        if ans == "n" 
            return 
        end
    end
    save(joinpath(@__DIR__, name * ".mrdi"), C)
end