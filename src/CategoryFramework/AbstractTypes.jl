#------------------------------------------------------------------------
#   Structs for categories
#------------------------------------------------------------------------

abstract type Category end

abstract type Object end

abstract type Morphism end


"""
    VectorSpaceObject

An object in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceObject <: Object end

"""
    VectorSpaceMorphism

A morphism in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceMorphism <: Morphism end

abstract type HomSet end

abstract type AbstractHomSpace <: VectorSpaceObject end

struct HomSpace <: AbstractHomSpace
    X::Object
    Y::Object
    basis::Vector{<:Morphism}
    parent
end

#=----------------------------------------------------------
    Getter & Setter 
----------------------------------------------------------=#

domain(m::Morphism) = m.domain
codomain(m::Morphism) = m.codomain

"""
    parent(X::Object)

Return the parent category of the object X.
"""
parent(X::Object) = X.parent

"""
    function parent(f::Morphism)

Return the parent category of ``f``.
"""
parent(f::Morphism) = parent(domain(f))

"""
    base_ring(X::Object)

Return the base ring ```k``` of the ```k```-linear parent category of ```X```.
"""
base_ring(X::Object) = base_ring(parent(X))
base_ring(X::Morphism) = base_ring(parent(domain(X)))

"""
    base_ring(C::Category)

Return the base ring ```k```of the ```k```-linear category ```C```.
"""
function base_ring(C::Category) 
    if hasfield(typeof(C), :base_ring) 
        return C.base_ring
    elseif hasfield(typeof(C), :category)
        return base_ring(category(C))
    else
        error("Category is not linear")
    end
end

base_group(C::Category) = C.base_group
base_group(X::Object) = parent(X).base_group

category(C::Category) = C.category
object(X::Object) = X.object
morphism(f::Morphism) = f.morphism

#=----------------------------------------------------------
    Functors
----------------------------------------------------------=#

abstract type AbstractFunctor  <: Object end


domain(F::AbstractFunctor) = F.domain
codomain(F::AbstractFunctor) = F.codomain


struct Functor <: AbstractFunctor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

#-------------------------------------------------------------------------------
#   Forgetful Functors
#-------------------------------------------------------------------------------

struct Forgetful <: AbstractFunctor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end



(F::AbstractFunctor)(x::T) where {T <: Object} = F.obj_map(x)
(F::AbstractFunctor)(x::T) where {T <: Morphism} = F.mor_map(x)

function show(io::IO, F::Forgetful)
    print(io, "Forgetful functor from $(domain(F)) to $(codomain(F))")
end
#-------------------------------------------------------------------------------
#   Hom Functors
#-------------------------------------------------------------------------------

struct HomFunctor <: AbstractFunctor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

function Hom(X::Object,::Colon)
    K = base_ring(parent(X))
    C = VectorSpaces(K)
    obj_map = Y -> Hom(X,Y)
    mor_map = f -> g -> g ∘ f
    return HomFunctor(parent(X),C,obj_map,mor_map)
end

function Hom(::Colon,X::Object)
    K = base_ring(parent(X))
    C = VectorSpaces(K)
    obj_map = Y -> Hom(Y,X)
    mor_map = g -> f -> g ∘ f
    return HomFunctor(OppositeCategory(parent(X)),C,obj_map,mor_map)
end

# function Hom(X::SetObject,::Colon)
#     obj_map = Y -> Hom(X,Y)
#     mor_map = f -> g -> g ∘ f
#     return HomFunctor(parent(X),Sets(),obj_map,mor_map)
# end

# function Hom(::Colon,X::SetObject)
#     obj_map = Y -> Hom(Y,X)
#     mor_map = g -> f -> g ∘ f
#     return HomFunctor(parent(X),Sets(),obj_map,mor_map)
# end

function show(io::IO, H::HomFunctor)
    print(io,"$(typeof(H.domain) == OppositeCategory ? "Contravariant" : "Covariant") Hom-functor in $(H.domain)")
end