"""
    hexagon_axiom(X::T, Y::T, Z::T) where T <: Object

Check the hexagon axiom for ```X, Y, Z```.
"""
function hexagon_axiom(X::T, Y::T, Z::T) where T <: Object
    f = associator(Y,Z,X) ∘ braiding(X,Y⊗Z) ∘ associator(X,Y,Z)
    g = (id(Y)⊗braiding(X,Z)) ∘ associator(Y,X,Z) ∘ (braiding(X,Y)⊗id(Z))
    return f == g
end

"""
    hexagon_axiom(objects::Vector{<:Object}, log::Bool = false)

Check the hexagon axiom for all combinations of objects in ```objects```. If
```log = true``` an array with the failing combinations is returned
"""
function hexagon_axiom(objects::Vector{<:Object}, log::Bool = false)
    failed = []
    for x ∈ objects, y ∈ objects, z ∈ objects

        hexagon_axiom(x,y,z) ? nothing : push!(failed, (x,y,z))
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
    hexagon_axiom(C::Category, log::Bool = false)

Check the hexagon axiom for all combinations of  simple objects of ```C```. If 
```log = true``` an array with the failing combinations is returned
"""
hexagon_axiom(C::Category, log::Bool = false) = hexagon_axiom(simples(C), log)
