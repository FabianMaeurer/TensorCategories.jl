function save(path::String, C::Category)
    try 
        mkdir(path)
    catch 
        rm(path, recursive = true)
        mkdir(path)
    end

    for field ∈ fieldnames(typeof(C))
        if isdefined(C,field)
            save(joinpath(path, "$field"), (getfield(C,field)))
            serialize(joinpath(path, "serialized_$field"), (getfield(C,field)))
        end
    end
end

function load(path::String, T::Type{<:Category})
    if ismutabletype(T)
        C = T()
        for field ∈ fieldnames(T)
            try
                @show fieldtype(T,field)
                val = load(joinpath(path, "$field"), fieldtype(T,field))
                setfield!(C, field, val)
            catch 
                try 
                    val = load(joinpath(path, "$field"))
                    setfield!(C, field, val)
                catch 
                    try
                        val = deserialize(joinpath(path, "serialized_$field"))
                        setfield!(C, field, val)
                    catch
                        continue
                    end
                end
            end
        end
        return C
    else
        values = [load(joinpath(path, "$field"), fieldtype(T,field)) for field ∈ fieldnames(T)]
        return T(values...)
    end
end

#=-------------------------------------------------
    Arrays 
-------------------------------------------------=#

function save(path::String, A::Array{<:Any,N}) where N
    try 
        mkdir(path)
    catch 
        rm(path, recursive = true)
        mkdir(path)
    end
    dims = size(A)
    iter = Base.product([1:k for k ∈ dims]...)
    for i ∈ iter
        save(joinpath(path, "$i"), A[i...])
        serialize(joinpath(path, "serialized_$i"), A[i...])
    end
end

function load(path::String, ::Type{Array{MatElem,N}}) where {N}
    string_indicies = readdir(path)
    indicies = [eval(Meta.parse(s)) for s ∈ string_indicies]

    n = length(indicies[1])
    if n != N throw(ErrorException("Type does not match")) end

    dims = [maximum(map(e->e[i], indicies)) for i ∈ 1:n]

    A = Array{MatElem,N}(undef,dims...)

    for ind ∈ indicies
        try 
            A[ind...] = load(joinpath(path, "$ind"), MatElem)
        catch
            try
                A[ind...] = load(joinpath(path, "$ind"))
            catch
                A[ind...] = deserialize(joinpath(path, """serialized_$ind"""))
            end
        end
    end
    return A
end