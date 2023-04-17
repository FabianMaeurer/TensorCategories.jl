const HashTypes = Union{CategoryMorphism, Category, CategoryObject, GrothendieckGroup, GrothendieckGroupElem}

function Base.hash(C::T, h::UInt) where T <: HashTypes
    content = (getfield(C, s) for s ∈ fieldnames(typeof(C)) if isdefined(C, s))
    hash(content, h)
end

function ==(X::T,Y::T) where T <: Union{CategoryMorphism, Category, CategoryObject}
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
