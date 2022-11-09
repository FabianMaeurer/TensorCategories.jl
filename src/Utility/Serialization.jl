function save(path::String, C::Category)
    try 
        mkdir(path)
    catch 
        rm(path, recursive = true)
        mkdir(path)
    end

    for field ∈ fieldnames(typeof(C))
        save(joinpath(path, "$field"), (getfield(C,field)))
    end
end

function load(path::AbstractString, T::Type{<:Category})
    if ismutabletype(T)
        C = T()
        for field ∈ fieldnames(T)
            try
                val = load(joinpath(path, "$field"), fieldtype(T,field))
                setfield!(C, field, val)
            catch 
                continue
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
    end
end

function load(path::String, T::Array{<:Any,N}) where N
    indecies = [readdir(path)
    A = [load(joinpath(path, "$(i)")) ]
    