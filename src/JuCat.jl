module JuCat

import Base: show,^,==, getindex, in, issubset, iterate, length,*,+,-, iterate,
                getproperty

import Oscar.AbstractAlgebra.Integers


import Oscar: VectorSpace, Field, elem_type, QQ, FieldElem,
                dim, base_ring, MatrixSpace, GAPGroup, GroupElem,
                ModuleIsomorphism, parent, matrix, basis, MatElem, ∘, gens,
                ⊕, compose, ⊗, tensor_product, Map, MatrixElem, kronecker_product,
                id, domain, one, zero, MatrixSpace, size, AbstractSet, inv, product,
                Ring, RingElem, base_field, MPolyQuo, iscommutative, isinvertible,
                MatrixGroup, hom, GAPGroupHomomorphism, GL, MatrixSpace, matrix,
                codomain, GAP, characteristic, degree, julia_to_gap, GSet, gset,
                FinField,gen, GSet, gset, orbits, stabilizer, orbit,
                isisomorphic, issubgroup, left_transversal, ComplexField, order,
                elements, index, symmetric_group, gap_to_julia, multiplication_table,
                issemisimple, AlgAss, AlgAssElem, FiniteField, abelian_closure,
                irreducible_modules, action




export Category, TensorCategory, Morphism, Object, VectorSpaces, base_ring, hom,
        GradedVectorSpaces, GradedVectorSpaceObject, GradedVectorSpaceMorphism,
        VectorSpaceObject, simples,
        VectorSpaceMorphism, parent, dsum,⊕, domain, codomain, compose, ∘, ^, ⊗,
        tensor_product,==, associator, basis, id, getindex, one, zero, Forgetful,
        Functor, Sets, SetObject, SetMorphism, inv, product, coproduct,
        features, issemisimple, isabelian, ismonoidal, ×, ∐, RepresentationCategory,
        GroupRepresentationCategory,
         FusionCategory, VSHomSpace,HomSpace,
        Hom, GVSHomSpace, HomFunctor, VSObject, GVSObject, GVSMorphism, SetHomSet,
        HomSet, Cocycle, trivial_3_cocylce,*, Algebra, FreeAlgebra, GroupAlgebra,
        AlgebraElem,+,-, GroupRepresentation, GroupRepresentationMorphism,
        isinvertible, Representation, isequivariant, matrix, GRHomSpace,
        OppositeCategory, OppositeMorphism, OppositeObject, ProductCategory,
        ProductObject, ProductMorphism, CohSheaves, CohSheaf, CohSheafMorphism,
        stalks, PullbackFunctor, Pullback, PushforwardFunctor, Pushforward,
        CohSfHomSpace, ConvolutionCategory, ConvolutionObject, ConvolutionMorphism,
        ConvHomSpace,stalk, induction, restriction, orbit_index, dsum_morphisms,
        decompose, multiplication_table, print_multiplication_table, groethendieck_ring


include("Utility/FFE_to_FinField.jl")
include("structures/abstracts.jl")
include("structures/MISC/ProductCategory.jl")
include("structures/MISC/OppositeCategory.jl")
include("structures/VectorSpaces/VectorSpaces.jl")
include("structures/VectorSpaces/Cocycles.jl")
include("structures/VectorSpaces/GradedVectorSpaces.jl")
include("structures/set.jl")
include("structures/Functors.jl")
include("structures/Representations/Representations.jl")
include("structures/Representations/AlgebraRepresentations.jl")
include("structures/Representations/GroupRepresentations.jl")
include("structures/AlgebraicStructures/Algebras.jl")
include("structures/AlgebraicStructures/AlgebraMorphisms.jl")
include("structures/ConvolutionCategory/CoherentSheaves.jl")
include("structures/ConvolutionCategory/ConvolutionCategory.jl")
include("structures/MISC/multiplication_table.jl")
include("structures/GroethendieckRing.jl")






end
