function Base.hash(C::T, h::UInt) where T <: Union{Morphism, Category, Object}
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
            
