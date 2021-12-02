module JuCat

import Base: show,^,==, getindex, in, issubset, iterate, length

import Oscar: VectorSpace, Field, elem_type, QQ, FieldElem,
                dim, base_ring, MatrixSpace, GAPGroup, GroupElem,
                ModuleIsomorphism, parent, matrix, basis, MatElem, ∘,
                ⊕, compose, ⊗, tensor_product, Map, MatrixElem, kronecker_product,
                id, domain, one, zero, MatrixSpace, size, AbstractSet, inv, product,
                Ring, RingElem, base_field, MPolyQuo


export Category, TensorCategory, Morphism, Object, VectorSpaces, base_ring, hom,
        GradedVectorSpaces, GradedVectorSpaceObject, GradedVectorSpaceMorphism,
        VectorSpaceObject,
        VectorSpaceMorphism, parent, dsum,⊕, domain, codomain, compose, ∘, ^, ⊗,
        tensor_product,==, associator, basis, id, getindex, one, zero, Forgetful,
        Functor, Sets, SetObject, SetMorphism, inv, product, coproduct,
        features, issemisimple, isabelian, ismonoidal, ×, ∐, RepresentationCategory,
        GroupRepresentationCategory, AlgebraRepresentationCategory,
        HopfAlgebraRepresentationCategory, FusionCategory, VSHomSpace,HomSpace,
        Hom, GVSHomSpace, HomFunctor, VSObject, GVSObject, GVSMorphism, SetHomSet,
        HomSet

include("structures/abstracts.jl")
include("structures/VectorSpaces.jl")
include("structures/GradedVectorSpaces.jl")
include("structures/Functors.jl")
include("structures/set.jl")
include("structures/Representations/Representations.jl")
include("structures/Representations/AlgebraRepresentations.jl")



end
