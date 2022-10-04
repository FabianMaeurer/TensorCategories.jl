module TensorCategories

import Base: show,^,==, getindex, in, issubset, iterate, length,*,+,-, iterate,
                getproperty

import Oscar.AbstractAlgebra.Integers

import Oscar:  Field, elem_type, QQ, FieldElem,
                dim, base_ring, MatrixSpace, GAPGroup, GroupElem,
                ModuleIsomorphism, parent, matrix, basis, MatElem, ∘, gens,
                ⊕, compose, ⊗, tensor_product, Map, MatrixElem, kronecker_product,
                id, domain, one, zero, MatrixSpace, size, AbstractSet,
                 inv, product,
                Ring, RingElem, base_field, MPolyQuo, iscommutative, isinvertible,
                MatrixGroup, hom, GAPGroupHomomorphism, GL, MatrixSpace, matrix,
                codomain, GAP, characteristic, degree, julia_to_gap, GSet, gset,
                FinField,gen, GSet, gset, orbits, stabilizer, orbit,
                isisomorphic, issubgroup, left_transversal, ComplexField, order,
                elements, index, symmetric_group, gap_to_julia, multiplication_table,
                issemisimple, AlgAss, AlgAssElem, FiniteField, abelian_closure,
                irreducible_modules, action, decompose,+, dual, tr, iscentral, rank,
                ZZ, solve_left, PolynomialRing, groebner_basis, ideal, roots,
                splitting_field, change_base_ring, isconstant, coeff, isindependent,
                coefficients, isabelian, leading_monomial, gcd, msolve, fmpz, fmpq,
                rref, NumberField, nf_elem, kernel, cokernel, primary_decomposition,
                Ideal, minpoly, image, solve, eigenspace, eigenspaces, diagonal_matrix,
                permutation_matrix, cyclotomic_field, jordan_normal_form

export Category, TensorCategory, Morphism, Object, VectorSpaces, base_ring, hom,
        GradedVectorSpaces, VectorSpaceObject, simples,
        VectorSpaceMorphism, parent, dsum,⊕, domain, codomain, compose, ∘, ^, ⊗,
        tensor_product,==, associator, basis, id, getindex, one, zero, Forgetful,
        Functor, Sets, SetObject, SetMorphism, inv, product, coproduct,
        features, issemisimple, isabelian, ismonoidal, ×, ∐, RepresentationCategory,
        GroupRepresentationCategory, ismultiring, ismultifusion, isring, ismultitensor,
        istensor, isfusion,
         FusionCategory, VSHomSpace,HomSpace,
        Hom, GVSHomSpace, HomFunctor, VSObject, GVSObject, GVSMorphism, SetHomSet,
        HomSet, Cocycle, trivial_3_cocylce,*,
        +,-, GroupRepresentation, GroupRepresentationMorphism,
        GroupRepresentationCategory,
        isinvertible, Representation, isequivariant, matrix, GRHomSpace,
        OppositeCategory, OppositeMorphism, OppositeObject, ProductCategory,
        ProductObject, ProductMorphism, CohSheaves, CohSheaf, CohSheafMorphism,
        stalks, PullbackFunctor, Pullback, PushforwardFunctor, Pushforward,
        CohSfHomSpace, ConvolutionCategory, ConvolutionObject, ConvolutionMorphism,
        ConvHomSpace,stalk, induction, restriction, orbit_index, dsum_with_morphisms,
        decompose, multiplication_table, print_multiplication_table, grothendieck_ring,
        dual, left_dual, right_dual, ev, coev, left_trace, right_trace, tr, braiding,
        drinfeld_morphism, smatrix, End, CenterCategory, CenterObject, CenterMorphism,
        spherical, iscentral, center_simples, RingCategory, set_tensor_product!,
        set_braiding!, set_one!, Ising, zero_morphism, express_in_basis, solve_groebner,
        Center, CenterCategory, CenterObject, CenterMorphism, ev_coev, matrices,
        orbit_stabilizers, GRepRestriction, GRepInduction, Restriction, Induction,
        print_multiplication_table, RingCatObject, RingCatMorphism, kernel, cokernel,
        image, isgraded, cyclic_group_3cocycle, decompose_morphism,
        central_objects, half_braiding, half_braidings, left_inverse, right_inverse,
        add_simple!, pentagon_axiom, set_associator!, dual_basis, pairing,
        eigenspace, eigenspaces, simple_subobjects, I2, I2subcategory,
        DeligneProdMorphism, DeligneProdObject, DeligneProduct, ⊠, op, AbstractHomSpace,
        is_half_braiding, object, distribute_right, distribute_left

include("Utility/FFE_to_FinField.jl")
include("Utility/SolveGroebner.jl")
include("Utility/Technicallities.jl")
include("structures/abstracts.jl")
include("structures/VectorSpaces/VectorSpaces.jl")
include("structures/MISC/ProductCategory.jl")
include("structures/MISC/OppositeCategory.jl")
include("structures/MISC/DeligneTensorProduct.jl")
include("structures/VectorSpaces/Cocycles.jl")
include("structures/VectorSpaces/GradedVectorSpaces.jl")
include("structures/set.jl")
include("structures/Representations/Representations.jl")
include("structures/Representations/GroupRepresentations.jl")
include("structures/Functors.jl")
include("structures/ConvolutionCategory/CoherentSheaves.jl")
include("structures/ConvolutionCategory/ConvolutionCategory.jl")
include("structures/MultiFusionCategories/FusionCategory.jl")
include("structures/MISC/multiplication_table.jl")
include("structures/MISC/PentagonAxiom.jl")
include("structures/GrothendieckRing.jl")
include("structures/Center/Induction.jl")
include("structures/Center/Center.jl")
include("structures/Center/CenterChecks.jl")
include("structures/Examples/I2-fusion.jl")


end
