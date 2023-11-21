#=----------------------------------------------------------
    Generic checks for categories 
----------------------------------------------------------=#

function is_fusion(C::Category) 
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :fusion) do 
            is_multifusion(C) && int_dim(End(one(C))) == 1
        end
    end
    is_multifusion(C) && int_dim(End(one(C))) == 1
end

function is_multifusion(C::Category) 
    if is_fusion(C) 
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :multifusion) do 
            is_multitensor(C) && is_semisimple(C) && is_finite(C)
        end
    end
    is_multitensor(C) && is_semisimple(C) && is_finite(C)
end

function is_tensor(C::Category) 
if is_fusion(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :tensor) do 
            is_multitensor(C) && int_dim(End(one(C))) == 1
        end
    end
    is_multitensor(C) && int_dim(End(one(C))) == 1
end

function is_multitensor(C::Category) 
if is_tensor(C) || is_multifusion(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :multitensor) do 
            is_multiring(C) && is_rigid(C)
        end
    end
    is_multiring(C) && is_rigid(C)
end

function is_ring(C::Category) 
if is_tensor(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :ring) do 
            is_multiring(C) && int_dim(End(one(C))) == 1
        end
    end
    is_multiring(C) && int_dim(End(one(C))) == 1
end

function is_multiring(C::Category) 
    if is_multitensor(C) || is_ring(C)
        return true
    end 
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :multiring) do 
            is_abelian(C) && is_linear(C) && is_monoidal(C)
        end
    end
    is_abelian(C) && is_linear(C) && is_monoidal(C)
end

function is_finite(C::Category) 
    if is_fusion(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :finite) do 
            try 
                return length(simples(C)) ≥ 0 
            catch end
            return false
        end
    end
    try 
        return length(simples(C)) ≥ 0 
    catch end
    return false
end

function is_monoidal(C::Category) 
    if is_multiring(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :monoidal) do 
            T = object_type(C)
            hasmethod(one, Tuple{typeof(C)}) &&
            hasmethod(tensor_product, Tuple{T,T}) 
        end
    end
    T = object_type(C)
    hasmethod(one, Tuple{typeof(C)}) &&
    hasmethod(tensor_product, Tuple{T,T}) 
end

function is_abelian(C::Category) 
    if is_multiring(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :abelian) do 
            if is_additive(C) && is_linear(C)
                T = morphism_type(C)
                return hasmethod(kernel, Tuple{T}) && hasmethod(cokernel, Tuple{T})
            end
            false
        end
    end
    if is_additive(C) && is_linear(C)
        T = morphism_type(C)
        return hasmethod(kernel, Tuple{T}) && hasmethod(cokernel, Tuple{T})
    end
    return false
end

function is_additive(C::Category) 
if is_abelian(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :additive) do 
            T = object_type(C)
            hasmethod(direct_sum, Tuple{T,T}) && hasmethod(zero, Tuple{typeof(C)})
        end
    end
    T = object_type(C)
    hasmethod(direct_sum, Tuple{T,T}) && hasmethod(zero, Tuple{typeof(C)})
end


function is_linear(C::Category) 
if is_multiring(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute(C, :additive) do 
            hasmethod(base_ring, Tuple{typeof(C)})
        end
    end
    hasmethod(base_ring, Tuple{typeof(C)})
end

is_semisimple(C::Category) = is_multifusion(C)

function is_modular(C::Category) 
    if hasfield(typeof(C), :__attrs) 
        return get_attribute!(C, :modular) do
            _is_modular(C)
        end
    end

    _is_modular(C)
end

function _is_modular(C::Category) 
    try
        return det(smatrix(C)) != 0 
    catch 
        return false 
    end
end

function is_spherical(C::Category)
    if hasfield(typeof(C), :__attrs) 
        return get_attribute!(C, :spherical) do
            _is_spherical(C)
        end
    end
    _is_spherical(C)
end

function _is_spherical(C::Category)
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
if is_multitensor(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute!(C, :spherical) do
            T = object_type(C)
            is_monoidal(C) && hasmethod(dual, Tuple{T}) && hasmethod(ev, Tuple{T}) && hasmethod(coev, Tuple{T})
        end
    end
    T = object_type(C)
    is_monoidal(C) && hasmethod(dual, Tuple{T}) && hasmethod(ev, Tuple{T}) && hasmethod(coev, Tuple{T})
end

function is_braided(C::Category)
    if hasfield(typeof(C), :__attrs) 
        return get_attribute!(C, :spherical) do
            T = object_type(C)
            is_monoidal(C) && hasmethod(braiding, Tuple{T,T})
        end
    end
    T = object_type(C)
    is_monoidal(C) && hasmethod(braiding, Tuple{T,T})
end

function is_krull_schmidt(C::Category)
if is_multiring(C)
        return true
    end
    if hasfield(typeof(C), :__attrs) 
        return get_attribute!(C, :spherical) do
            false
        end
    end
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