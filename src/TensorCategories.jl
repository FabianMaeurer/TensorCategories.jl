module TensorCategories

import Base: *, +, -, ==, ^, getindex, getproperty, in, issubset, iterate, length, show

import Oscar: +, AbstractSet, AlgAss, AlgAssElem, ComplexField, Field, FieldElem, FinField,
    FiniteField, GAP, GAPGroup, GAPGroupHomomorphism, GL, GSet, GroupElem, Ideal, MPolyQuo,
    Map, MatElem, MatrixElem, MatrixGroup, MatrixSpace, ModuleIsomorphism, NumberField,
    PolynomialRing, QQ, Ring, RingElem, ZZ, abelian_closure, action, base_field, base_ring,
    basis, change_base_ring, characteristic, codomain, coeff, coefficients, cokernel,
    compose, cyclotomic_field, decompose, degree, diagonal_matrix, dim, domain, dual,
    eigenspace, eigenspaces, elem_type, elements, fmpq, fmpz, gap_to_julia, gcd, gen, gens,
    groebner_basis, gset, hom, id, ideal, image, index, inv, irreducible_modules, isabelian,
    is_central, iscommutative, isconstant, isindependent, isinvertible, isisomorphic,
    issemisimple, issubgroup, jordan_normal_form, julia_to_gap, kernel, kronecker_product,
    leading_monomial, left_transversal, matrix, minpoly, msolve, multiplication_table,
    nf_elem, one, orbit, orbits, order, parent, permutation_matrix, primary_decomposition,
    product, rank, roots, rref, size, solve, solve_left, splitting_field, stabilizer,
    symmetric_group, tensor_product, tr, zero, ∘, ⊕, ⊗, iso_oscar_gap, preimage, is_simple,
    CyclotomicField, absolute_simple_field, is_abelian, is_square, charpoly, det, load,save,
    factor, zero_matrix, identity_matrix, complex_embeddings, QQBar, eigenvalues, @alias


using Memoization

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
        drinfeld_morphism, smatrix, tmatrix, End, CenterCategory, CenterObject, CenterMorphism,
        spherical, is_central, center_simples, RingCategory, set_tensor_product!,
        set_braiding!, set_one!, Ising, zero_morphism, express_in_basis, solve_groebner,
        Center, CenterCategory, CenterObject, CenterMorphism, ev_coev, matrices,
        orbit_stabilizers, GRepRestriction, GRepInduction, Restriction, Induction,
        print_multiplication_table, RingCatObject, RingCatMorphism, kernel, cokernel,
        image, isgraded, cyclic_group_3cocycle, decompose_morphism,
        central_objects, half_braiding, half_braidings, left_inverse, right_inverse,
        add_simple!, pentagon_axiom, set_associator!, dual_basis, pairing,
        eigenspace, eigenspaces, simple_subobjects, I2, I2subcategory,
        DeligneProdMorphism, DeligneProdObject, DeligneProduct, ⊠, op, AbstractHomSpace,
        is_half_braiding, object, distribute_right, distribute_left, is_simple,
        decompose_morphism, TambaraYamagami, RingSubcategory, SubcategoryMorphism,
        SubcategoryObject,load,save, cat_fr_8122, E6subfactor, fpdim, set_cannonical_spherical!,
        normalized_smatrix, sort_simples_by_dimension!, set_spherical!, inv_associator

        GAP.Packages.install("HAP")
        GAP.Packages.load("HAP") 

include("structures/abstracts.jl")
include("Utility/Technicallities.jl")
include("structures/VectorSpaces/VectorSpaces.jl")
include("structures/MISC/ProductCategory.jl")
include("structures/MISC/OppositeCategory.jl")
include("structures/MISC/DeligneTensorProduct.jl")
include("structures/VectorSpaces/Cocycles.jl")
include("structures/VectorSpaces/Unitary-3-cocycle.jl")
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
include("structures/MISC/Subcategories.jl")
include("structures/GrothendieckRing.jl")
include("structures/Center/Center.jl")
include("structures/Center/Induction.jl")
include("structures/Center/CenterChecks.jl")
include("structures/Examples/I2-fusion.jl")
include("structures/Examples/TambaraYamagami.jl")
include("structures/Examples/VercleyenSingerland/vercleyen_singerland.jl")
include("structures/Examples/E6subfactor.jl")

#include("Utility/Serialization.jl")

@alias Centre Center

end
