#=----------------------------------------------------------
    Generic checks for categories 
----------------------------------------------------------=#

# Structural predicates are conservative capability queries. A true result
# records declared or category-specific knowledge and includes the logical
# consequences below. A false result may also mean that a property has not been
# declared; generic methods cannot prove axioms from the existence of methods.
_declared_structure(C::Category, key::Symbol) =
    hasfield(typeof(C), :__attrs) && get_attribute(C, key, false) === true

is_fusion(C::Category) = _declared_structure(C, :fusion)
is_multifusion(C::Category) = is_fusion(C) || _declared_structure(C, :multifusion)
is_weak_fusion(C::Category) = is_fusion(C) || _declared_structure(C, :weak_fusion)
is_weak_multifusion(C::Category) = is_multifusion(C) || is_weak_fusion(C) ||
    _declared_structure(C, :weak_multifusion)

function is_split_semisimple(C::Category)
    is_multifusion(C) && return true
    is_semisimple(C) && all(s -> int_dim(End(s)) == 1, simples(C))
end

is_tensor(C::Category) = is_weak_fusion(C) || _declared_structure(C, :tensor)
is_multitensor(C::Category) = is_tensor(C) || is_weak_multifusion(C) ||
    _declared_structure(C, :multitensor)
is_ring(C::Category) = is_tensor(C) || _declared_structure(C, :ring)
is_multiring(C::Category) = is_multitensor(C) || is_ring(C) ||
    _declared_structure(C, :multiring)

"""
    is_finite(C::Category)

Return whether `C` is known to be a finite abelian category.
"""
is_finite(C::Category) = is_weak_multifusion(C) || _declared_structure(C, :finite)

"""
    is_locally_finite(C::Category)

Return whether `C` is known to be a locally finite linear abelian category.
"""
is_locally_finite(C::Category) = is_finite(C) || is_multiring(C) ||
    _declared_structure(C, :locally_finite)

is_monoidal(C::Category) = is_multiring(C) || any(
    key -> _declared_structure(C, key), (:monoidal, :rigid, :spherical, :is_braided))
is_semisimple(C::Category) = is_weak_multifusion(C) ||
    _declared_structure(C, :semisimple)
is_abelian(C::Category) = is_locally_finite(C) || is_semisimple(C) ||
    _declared_structure(C, :abelian)
is_krull_schmidt(C::Category) = is_locally_finite(C) ||
    _declared_structure(C, :krull_schmidt)
is_additive(C::Category) = is_abelian(C) || is_krull_schmidt(C) ||
    _declared_structure(C, :additive)
is_linear(C::Category) = is_locally_finite(C) || _declared_structure(C, :linear)

function is_modular(C::Category) 
    if hasfield(typeof(C), :__attrs) 
        return get_attribute!(C, :modular) do
            _is_modular(C)
        end
    end

    _is_modular(C)
end

function _is_modular(C::Category) 
    is_fusion(C) && is_braided(C) && is_spherical(C) || return false
    d = det(smatrix(C))
    if base_ring(C) isa Union{ArbField,AcbField,ComplexField}
        return !Oscar.contains_zero(d)
    end
    !iszero(d)
end

"""
    is_spherical(C::Category; check=false)

Return the declared or cached spherical status. Pass `check=true` to verify the
chosen pivotal structure and equality of left and right dimensions.
"""
function is_spherical(C::Category; check::Bool=false)
    check && return _is_spherical(C)
    if hasfield(typeof(C), :__attrs) 
        return get_attribute!(C, :spherical) do
            _is_spherical(C)
        end
    end
    _is_spherical(C)
end

function _is_spherical(C::Category)
    # EGNO, Section 4.7: for a split semisimple pivotal category it suffices
    # to compare the left and right dimensions on simple objects. A method
    # producing components does not by itself prove pivotal monoidality.
    is_split_semisimple(C) || return false
    S = simples(C)
    all(X -> applicable(spherical,X),S) || return false
    is_pivotal(C;check=true) || return false
    if base_ring(C) isa Union{ArbField,AcbField,ComplexField}
        return all(overlaps(dim(X),dim(dual(X))) for X in S)
    end
    all(dim(X) == dim(dual(X)) for X in S)
end

function is_rigid(C::Category)
    is_multitensor(C) || _declared_structure(C, :rigid) ||
        _declared_structure(C, :spherical)
end

is_braided(C::Category) = _declared_structure(C, :is_braided)

is_unitary(C::Category) = false
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
            if typeof(C) == fieldtype(T,:parent)
                return T
            end
        end
    end
end

function morphism_type(C::Category)
    morphism_types = all_subtypes(Morphism)

    for T ∈ morphism_types
        if hasfield(T, :domain)
            if object_type(C) == fieldtype(T,:domain)
                return T
            end
        end
    end
end
