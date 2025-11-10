"""
    pentagon_axiom(X::T, Y::T, Z::T, W::T) where T <: Object

Check the pentagon axiom for ```X, Y, Z, W```.
"""
function pentagon_axiom(X::T, Y::T, Z::T, W::T) where T <: Object
    if typeof(base_ring(X)) <: Union{ArbField, ComplexField, AcbField}
        return pentagon_axiom_numeric(X, Y, Z, W)
    end 

    f = (id(X)⊗associator(Y,Z,W)) ∘ associator(X,Y⊗Z,W) ∘ (associator(X,Y,Z)⊗id(W))
    g = associator(X,Y,Z⊗W) ∘ associator(X⊗Y,Z,W)
    return f == g
end

function pentagon_axiom_numeric(X::T, Y::T, Z::T, W::T) where T <: Object
    f = (id(X)⊗associator(Y,Z,W)) ∘ associator(X,Y⊗Z,W) ∘ (associator(X,Y,Z)⊗id(W))
    g = associator(X,Y,Z⊗W) ∘ associator(X⊗Y,Z,W)
    return overlaps(matrix(f), matrix(g))
end

"""
    pentagon_axiom(objects::Vector{<:Object}, log::Bool = false)

Check the pentagon axiom for all combinations of objects in ```objects```. If
```log = true``` an array with the failing combinations is returned
"""
function pentagon_axiom(objects::Vector{<:Object}, log::Bool = false; show_progress = false)
    failed = []
    N = (length(objects)-1)^4
    checked = 0

    Threads.@threads for x ∈ objects

        x == one(parent(x)) && continue
        
        for y ∈ objects,  z ∈ objects, w ∈ objects
            one(parent(y)) ∈ [y,z,w] && continue

            pentagon_axiom(x,y,z,w) ? nothing : push!(failed, (x,y,z,w))
                
            if !log && length(failed) > 0 
                show_progress && print("\n")
                return false
            end

            if show_progress
                checked += 1
                print("\e[2K\rChecked $(checked) / $N combinations ($(round(checked / N * 100, digits = 2))%)")
            end
        end
    end
    show_progress && print("\n")

    if log
        return length(failed) == 0, failed
    else
        return length(failed) == 0
    end
end

function randomized_pentagon_axiom(C::Category, n::Int = 0)
    S = simples(C)
    m = length(S)
    if n == 0
        n = m^2
    end

    for _ ∈ 1:n
        if !pentagon_axiom(C[rand(2:m, 4)]...)
            return false
        end
    end

    true
end

"""
    pentagon_axiom(C::Category, log::Bool = false)

Check the pentagon axiom for all combinations of  simple objects of ```C```. If 
```log = true``` an array with the failing combinations is returned
"""
pentagon_axiom(C::Category, log::Bool = false; show_progress = false) = pentagon_axiom(simples(C), log, show_progress = show_progress)
