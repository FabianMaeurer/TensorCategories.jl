"""
    pentagon_axiom(X::T, Y::T, Z::T, W::T) where T <: Object

Check the pentagon axiom for ```X, Y, Z, W```.
"""
function pentagon_axiom(X::T, Y::T, Z::T, W::T) where T <: Object
    f = (id(X)⊗associator(Y,Z,W)) ∘ associator(X,Y⊗Z,W) ∘ (associator(X,Y,Z)⊗id(W))
    g = associator(X,Y,Z⊗W) ∘ associator(X⊗Y,Z,W)
    return f == g
end

"""
    pentagon_axiom(objects::Vector{<:Object}, log::Bool = false)

Check the pentagon axiom for all combinations of objects in ```objects```. If
```log = true``` an array with the failing combinations is returned
"""
function pentagon_axiom(objects::Vector{<:Object}, log::Bool = false)
    failed = []
    for x ∈ objects, y ∈ objects,
        z ∈ objects, w ∈ objects

        pentagon_axiom(x,y,z,w) ? nothing : push!(failed, (x,y,z,w))
        if !log && length(failed) > 0 
            return false
        end
    end

    if log
        return length(failed) == 0, failed
    else
        return length(failed) == 0
    end
end

"""
    pentagon_axiom(C::Category, log::Bool = false)

Check the pentagon axiom for all combinations of  simple objects of ```C```. If 
```log = true``` an array with the failing combinations is returned
"""
pentagon_axiom(C::Category, log::Bool = false) = pentagon_axiom(simples(C), log)
