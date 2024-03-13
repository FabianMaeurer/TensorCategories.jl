const HashTypes = Union{Morphism, Category, Object}

function Base.hash(C::T, h::UInt) where T <: HashTypes
    content = (getfield(C, s) for s ∈ fieldnames(typeof(C)) if isdefined(C, s))
    hash(content, h)
end

function ==(X::T,Y::T) where T <: Union{Morphism, Category, Object}
    for s ∈ fieldnames(typeof(X)) 
        if (isdefined(X, s) ⊻ isdefined(Y, s)) 
            return false
        elseif isdefined(X,s) && isdefined(Y,s)
            if getfield(X,s) != getfield(Y,s)
                return false
            end
        end
    end
    return true
end

function unique_without_hash(A::AbstractArray{T,N}) where {T,N}
    if length(A) ≤ 1 
        return A
    end
    B = T[]
    A2 = deepcopy(A)
    while length(A2) > 0
        f = popfirst!(A2)
        B = [B; f]
        A2 = [g for g ∈ A2 if g != f]
    end
    return B
end

            
#=----------------------------------------------------------
    Replace terms in Expression 
----------------------------------------------------------=#

function replace!(e, old, new)
    for (i,a) in enumerate(e.args)
        if a==old
            e.args[i] = new
        elseif a isa Expr
            replace!(a, old, new)
        end
        ## otherwise do nothing
    end
    e
end


#=----------------------------------------------------------
    Dicts    
----------------------------------------------------------=#

function Base.getindex(D::Dict{<:Object,<:Any}, X::Object)
    i = findfirst(e -> e == X, collect(keys(D)))

    i === nothing && KeyError("key $X not found")

    return collect(values(D))[i]
end