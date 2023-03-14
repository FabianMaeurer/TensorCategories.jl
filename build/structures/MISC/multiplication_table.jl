function multiplication_table(C::Category, simples::Vector{<:CategoryObject} = simples(C))
    @assert is_semisimple(C) "Category needs to be semi-simple"
    m = [s⊗t for s ∈ simples, t ∈ simples]
    coeffs = [coefficients(m,simples) for m ∈ m]
    return [c[k] for c ∈ coeffs, k ∈ 1:length(simples)]
end

function multiplication_table(simples::Vector{<:CategoryObject})
    @assert is_semisimple(parent(simples[1])) "Category needs to be semi-simple"

    return multiplication_table(parent(simples[1]), simples)
end

function print_multiplication_table(simples::Vector{<:CategoryObject}, names::Vector{String} = ["v$i" for i ∈ 1:length(simples)])
    @assert length(simples) == length(names) "Invalid input"
    mult_table = multiplication_table(parent(simples[1]), simples)

    return [pretty_print_semisimple(s⊗t,simples,names) for s ∈ simples, t ∈ simples]
end

print_multiplication_table(C::Category) = print_multiplication_table(simples(C), simples_names(C))

function pretty_print_semisimple(m::CategoryObject,simples::Vector{<:CategoryObject},names::Vector{String})
    facs = decompose(m, simples)

    if length(facs) == 0 return "0" end

    str = ""
    for (o,k) ∈ facs
        i = findfirst(x -> x == o, simples)
        if i == nothing i = findfirst(x -> is_isomorphic(x,o)[1], simples) end

        str = length(str) > 0 ? str*" ⊕ "* (k > 1 ? "$(names[i])^$k" : "$(names[i])") : (k > 1 ? str*"$(names[i])^$k" : str*"$(names[i])")
    end
    return str
end


function coefficients(X::T, simples::Vector{T} = simples(parent(X))) where {T <: CategoryObject}
    facs = decompose(X)
    coeffs = [0 for i ∈ 1:length(simples)]
    for (x,k) ∈ facs
        i = findfirst(y -> y == x, simples)
        if i == nothing i = findfirst(y -> is_isomorphic(y,x)[1], simples) end
        coeffs[i] = k
    end
    return coeffs
end
