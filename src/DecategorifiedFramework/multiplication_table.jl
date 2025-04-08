function multiplication_table(C::Category, indecomposables::Vector{<:Object} = indecomposables(C))
    #@assert is_multitensor(C) "Category needs to be multitensor"
    #m = [s⊗t for s ∈ indecomposables, t ∈ indecomposables]
    
    if hasfield(typeof(C), :__attrs) && has_attribute(C, :multiplication_table)
        return get_attribute(C, :multiplication_table)
    end
    
    coeffs = [coefficients(s⊗t, indecomposables) for s ∈ indecomposables, t ∈ indecomposables]

        
    mult = [c[k] for c ∈ coeffs, k ∈ 1:length(indecomposables)]
    if hasfield(typeof(C), :__attrs)
        set_attribute!(C, :multiplication_table, mult)
        return mult
    else 
        return mult
    end
end

function multiplication_table(indecomposables::Vector{<:Object})
    #@assert is_semisimple(parent(indecomposables[1])) "Category needs to be semi-simple"

    return multiplication_table(parent(indecomposables[1]), indecomposables)
end

function print_multiplication_table(indecomposables::Vector{<:Object}, names::Vector{String} = ["v$i" for i ∈ 1:length(indecomposables)])
    @assert length(indecomposables) == length(names) "Invalid input"
    #mult_table = multiplication_table(parent(indecomposables[1]), indecomposables)

    return [pretty_print_decomposable(s⊗t,indecomposables,names) for s ∈ indecomposables, t ∈ indecomposables]
end

function print_multiplication_table(C::Category)
    if is_semisimple(C) 
        print_multiplication_table(multiplication_table(C), simples_names(C))
    else
        print_multiplication_table(indecomposables(C), indecomposables_names(C))
    end
end

function pretty_print_decomposable(m::Object,indecomposables::Vector{<:Object},names::Vector{String})
    facs = is_semisimple(parent(m)) ? decompose(m, indecomposables) : decompose(m)

    if length(facs) == 0 return "0" end

    str = ""
    for (o,k) ∈ facs
        i = findfirst(x -> x == o, indecomposables)
        if i == nothing i = findfirst(x -> is_isomorphic(x,o)[1], indecomposables) end

        str = length(str) > 0 ? str*" ⊕ "* (k > 1 ? "$k⋅$(names[i])" : "$(names[i])") : (k > 1 ? str*"$k⋅$(names[i])" : str*"$(names[i])")
    end
    return str
end

function pretty_print_decomposable(facs::Vector{Int}, names::Vector{String})
    
    str = ""
    
    for (k,s) ∈ zip(facs, names)
        
        k == 0 && continue
        str = length(str) > 0 ? str*" ⊕ "* (k > 1 ? "$k⋅$(s)" : "$(s)") : (k > 1 ? str*"$k⋅$(s)" : str*"$(s)")
    end
    return str
end

function coefficients(X::T, indecomposables::Vector{T} = indecomposables(parent(X))) where {T <: Object}
    if is_semisimple(parent(X))
        return fusion_coefficients(X, indecomposables)
    else
        return _coefficients(X, indecomposables)
    end
end

function fusion_coefficients(X::T, simpls::Vector{T} = simples(parent(X))) where {T <: Object}
    [int_dim(Hom(s,X)) for s ∈ simpls]
end

function _coefficients(X::T, indecomposables::Vector{T} = indecomposables(parent(X))) where {T <: Object}
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


function print_multiplication_table(A::Array{T,3}, names::Vector{String} =    ["X$i" for i ∈ 1:size(A)[1]]) where T <: Union{Int, RingElem}
    n = size(A)[1]
    [print_sum(A[i,j,:], names) for i ∈ 1:n, j ∈ 1:n]
end
        
    
function print_sum(facs::Vector{T}, names::Vector{String}) where T <: Union{Int, RingElem}
    if sum(facs) == 0 return "0" end

    str = ""
    for (k,name) ∈ zip(facs, names)
        if k == 0 
            continue
        end

        str = length(str) > 0 ? str*" ⊕ "* (k ∉ [0,1] ? "$k⋅$(name)" : "$(name)") : (k ∉ [0,1] ? str*"$k⋅$(name)" : str*"$(name)")
    end
    return str
end


#=----------------------------------------------------------
    print module action 
----------------------------------------------------------=#

function print_module_action(M::T, names = nothing) where T <: Union{RightModuleCategory, LeftModuleCategory}
    @assert is_semisimple(M)

    simples_C = simples(category(M))
    simples_M = simples(M)
    names = names === nothing ? ["m$i" for i ∈ 1:length(simples_M)] : names

    reshape(vcat([[print_sum(collect(r), names) for r ∈ eachrow(action_matrix(X,M))] for X ∈ simples_C]...), length(simples_C), length(simples_M))
end

#=----------------------------------------------------------
    Live progress 
----------------------------------------------------------=#

function multiplication_table_with_progress(C::Category, indecs::Vector{<:Object} = indecomposables(C); symmetric::Bool = false, names = simples_names(C), one = nothing)

    n = length(indecs)
    m = Array{Int}(undef,n,n,n)

    displ = ["⋅" for _ ∈ 1:n, _ ∈ 1:n]

    for i ∈ 1:n, j ∈ (symmetric ? i : 1):n

        if one != nothing && (i ∈ one || j in one)
            k = i ∈ one ? j : i
            m[i,j,:] = m[j,i,:] = [l == k for l ∈ 1:n]
        else
            m[i,j,:] = coefficients(indecs[i] ⊗ indecs[j], indecs)
            if symmetric m[j,i] = m[i,j] end
        end

        displ[i,j] = pretty_print_decomposable(m[i,j,:], names)
        if symmetric displ[j,i] = displ[i,j] end

        if i*j != 1 
            for _ ∈ 1:n+1
                print("\e[A\e[2K")
            end
        end

        display(displ)
    end

    println("")
    return m
end