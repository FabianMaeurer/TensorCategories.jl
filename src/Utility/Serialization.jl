encodeType(::Type{CenterCategory}) = "TensorCategories.CenterCategory"

function save_internal(s::SerializerState, C::CenterCategory)
    if isdefined(C,:simples)
        return Dict(
            :base_ring => save_type_dispatch(s, C.base_ring),
            :category => save_type_dispatch(s, C.category),
            :simples => save_type_dispatch(s, C.simples)
        )
    else
        return Dict(
            :base_ring => save_type_dispatch(s, C.base_ring),
            :category => save_type_dispatch(s, C.category)
        )
    end
end

function load_internal(s::DeserializerState, ::Type{CenterCategory}, dict::Dict)
    C = CenterCategory(dict[:base_ring], dict[:category])
    if :simples ∈ keys(dict)
        C.simples = dict[:simples]
    end
    return C
end

encodeType(T::Type{<:Category}) = "TensorCategories.$T"

function save_internal(s::SerializerState, C::Category)
    fields = [s for s ∈ fieldnames(typeof(C)) if isdefined(C,s)]
    return Dict(v => save_type_dispatch(s, getfield(C,v)) for v ∈ fields)
end

