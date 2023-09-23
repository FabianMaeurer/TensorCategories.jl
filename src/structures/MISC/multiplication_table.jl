function multiplication_table(C::Category, indecomposables::Vector{<:Object} = indecomposables(C))
    @assert is_multitensor(C) "Category needs to be multitensor"
    #m = [s⊗t for s ∈ indecomposables, t ∈ indecomposables]
    coeffs = [coefficients(s⊗t, indecomposables) for s ∈ indecomposables, t ∈ indecomposables]
    return [c[k] for c ∈ coeffs, k ∈ 1:length(indecomposables)]
end

function multiplication_table(indecomposables::Vector{<:Object})
    @assert is_semisimple(parent(indecomposables[1])) "Category needs to be semi-simple"

    return multiplication_table(parent(indecomposables[1]), indecomposables)
end

function print_multiplication_table(indecomposables::Vector{<:Object}, names::Vector{String} = ["v$i" for i ∈ 1:length(indecomposables)])
    @assert length(indecomposables) == length(names) "Invalid input"
    mult_table = multiplication_table(parent(indecomposables[1]), indecomposables)

    return [pretty_print_semisimple(s⊗t,indecomposables,names) for s ∈ indecomposables, t ∈ indecomposables]
end

print_multiplication_table(C::Category) = print_multiplication_table(indecomposables(C), indecomposables_names(C))

function pretty_print_decomposable(m::Object,indecomposables::Vector{<:Object},names::Vector{String})
    facs = decompose(m, indecomposables)

    if length(facs) == 0 return "0" end

    str = ""
    for (o,k) ∈ facs
        i = findfirst(x -> x == o, indecomposables)
        if i == nothing i = findfirst(x -> is_isomorphic(x,o)[1], indecomposables) end

        str = length(str) > 0 ? str*" ⊕ "* (k > 1 ? "$(names[i])^$k" : "$(names[i])") : (k > 1 ? str*"$(names[i])^$k" : str*"$(names[i])")
    end
    return str
end


function coefficients(X::T, indecomposables::Vector{T} = indecomposables(parent(X))) where {T <: Object}
    facs = decompose(X)
    coeffs = [0 for i ∈ 1:length(indecomposables)]
    for (x,k) ∈ facs
        i = findfirst(y -> y == x, indecomposables)
        if i === nothing 
            i = findfirst(y -> is_isomorphic(y,x)[1], indecomposables) 
        end
        coeffs[i] = k
    end
    return coeffs
end
