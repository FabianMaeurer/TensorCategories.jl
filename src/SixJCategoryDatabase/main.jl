#=----------------------------------------------------------
    Access Fusion Category Database  
----------------------------------------------------------=#

function load_fusion_category(name::String)
    load(joinpath(@__DIR__, name * ".mrdi"))
end

function add_to_local_database(C::SixJCategory, name::String, path::String = @__DIR__)

    files = readdir(joinpath(path))

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
    save(joinpath(path, name * ".mrdi"), C)
end