#=----------------------------------------------------------
    Generic checks for categories 
----------------------------------------------------------=#

function is_fusion(C::Category) 
    try 
        return get_attribute(C, :fusion)
    catch end
    is_multifusion(C) && int_dim(End(one(C))) == 1
end

function is_multifusion(C::Category) 
    try
        return get_attribute(C, :fusion)
    catch end
    is_multitensor(C) && is_semisimple(C) && is_finite(C)
end

function is_tensor(C::Category) 
    try 
        return get_attribute(C, :tensor)
    catch end
    is_multitensor(C) && int_dim(End(one(C))) == 1
end

function is_multitensor(C::Category) 
    try 
        return get_attribute(C, :multitensor)
    catch end
    is_multiring(C) && is_rigid(C)
end

function is_ring(C::Category) 
    try 
        get_attribute(C, :ring)
    catch end
    is_multiring(C) && int_dim(End(one(C))) == 1
end

function is_multiring(C::Category) 
    try 
        get_attribute(C, :multiring)
    catch end
    is_abelian(C) && is_linear(C) && is_monoidal(C)
end

function is_finite(C::Category) 
    try 
        return length(simples(C)) ≥ 0 
    catch end
    return false
end

function is_monoidal(C::Category) 
    try 
        return get_attribute(C, :monoidal)
    catch end
    T = object_type(C)
    hasmethod(one, Tuple{typeof(C)}) &&
    hasmethod(tensor_product, Tuple{T,T}) 
end

function is_abelian(C::Category) 
    try 
        return get_attribute(C, :abelian)
    catch end
    if is_additive(C) && is_linear(C)
        T = morphism_type(C)
        return hasmethod(kernel, Tuple{T}) && hasmethod(cokernel, Tuple{T})
    end
    return false
end

function is_additive(C::Category) 
    try 
        return get_attribute(C, :additive)
    catch end
    T = object_type(C)
    hasmethod(direct_sum, Tuple{T,T}) && hasmethod(zero, Tuple{typeof(C)})
end


function is_linear(C::Category) 
    try 
        return get_attribute(C, :linear)
    catch end
    hasmethod(base_ring, Tuple{typeof(C)})
end

is_semisimple(C::Category) = is_multitensor(C)

function is_modular(C::Category) 
    try 
        return get_attribute(C, :modular)
    catch end

    try
        return det(smatrix(C)) != 0 
    catch 
        return false 
    end
end

function is_spherical(C::Category)
    try 
        return get_attribute(C, :spherical)
    catch end

    @assert is_multifusion(C) "Generic checking only available for multifusion categories"

    obj_type = typeof(one(C))
    if  !hasmethod(spherical, Tuple{obj_type})
        return false
    end
    try 
        for x ∈ simples(C)
            spherical(x)
        end
        return true
    catch
        return false
    end
end

function is_rigid(C::Category)
    try 
        return get_attribute(C, :rigid)
    catch end
    T = object_type(C)
    is_monoidal(C) && hasmethod(dual, Tuple{T}) && hasmethod(ev, Tuple{T}) && hasmethod(coev, Tuple{T})
end

function is_braided(C::Category)
    try 
        return get_attribute(C, :braided)
    catch end
    T = object_type(C)
    is_monoidal(C) && hasmethod(braiding, Tuple{T,T})
end

function is_krull_schmidt(C::Category)
    try 
        return get_attribute(C, :krull_schmidt)
    catch end
    false
end


#=----------------------------------------------------------
    Helpers 
----------------------------------------------------------=#

function all_subtypes(T::Type)
    sub_types = subtypes(T)
    
    is_abstract = isabstracttype.(sub_types)

    concrete_types = sub_types[true .⊻ (is_abstract)]
    abstract_types = sub_types[is_abstract]

    return [concrete_types; vcat(all_subtypes.(abstract_types))...]
end


function object_type(C::Category)
    object_types = all_subtypes(Object)

    for T ∈ object_types
        if hasfield(T, :parent)
            if typeof(C) <: fieldtype(T,:parent)
                return T
            end
        end
    end
end 

function morphism_type(C::Category)
    morphism_types = all_subtypes(Morphism)

    for T ∈ morphism_types
        if hasfield(T, :domain)
            if object_type(C) <: fieldtype(T,:domain)
                return T
            end
        end
    end
end 