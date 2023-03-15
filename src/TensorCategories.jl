module TensorCategories

import Base: *, +, -, ==, ^, getindex, getproperty, in, issubset, iterate, length, show

import Oscar: +, AbstractSet, AlgAss, AlgAssElem, ComplexField, Field, FieldElem, FinField,
    FiniteField, GAP, GAPGroup, GAPGroupHomomorphism, GL, GSet, GroupElem, Ideal, MPolyQuo,
    Map, MatElem, MatrixElem, MatrixGroup, MatrixSpace, ModuleIsomorphism, NumberField,
    PolynomialRing, QQ, Ring, RingElem, ZZ, abelian_closure, action, base_field, base_ring,
    basis, change_base_ring, characteristic, codomain, coeff, coefficients, cokernel,
    compose, cyclotomic_field, decompose, degree, diagonal_matrix, dim, domain, dual,
    eigenspace, eigenspaces, elem_type, elements, fmpq, fmpz, gap_to_julia, gcd, gen, gens,
    groebner_basis, gset, hom, id, ideal, image, index, inv, irreducible_modules, is_abelian,
    is_central, iscommutative, isconstant, isindependent, isinvertible, is_isomorphic,
    is_semisimple, issubgroup, jordan_normal_form, julia_to_gap, kernel, kronecker_product,
    leading_monomial, left_transversal, matrix, minpoly, real_solutions, multiplication_table,
    nf_elem, one, orbit, orbits, order, parent, permutation_matrix, primary_decomposition,
    product, rank, roots, rref, size, solve, solve_left, splitting_field, stabilizer,
    symmetric_group, tensor_product, tr, zero, ∘, ⊕, ⊗, iso_oscar_gap, preimage, is_simple,
    CyclotomicField, absolute_simple_field, is_abelian, is_square, charpoly, det, load,save,
    factor, zero_matrix, identity_matrix, complex_embeddings, QQBar, eigenvalues, @alias,
    abelian_group, PcGroup, is_modular, subgroup, nullspace, factor, qqbar,
    leading_coefficient, roots, is_rational, fmpq_mpoly, lex, Fac, root_of_unity, PolyElem, MPolyElem, monomials, fmpq_poly, MPolyIdeal,
    height_bits, lcm, change_base_ring, guess, direct_sum



using Memoization, InteractiveUtils

export Category, TensorCategory, CategoryMorphism, CategoryObject, VectorSpaces, base_ring, hom, Morphism,
        GradedVectorSpaces, VectorSpaceCategoryObject, simples,
        VectorSpaceCategoryMorphism, parent, direct_sum,⊕, domain, codomain, compose, ∘, ^, ⊗,
        tensor_product,==, associator, basis, id, getindex, one, zero, Forgetful,
        Functor, Sets, SetCategoryObject, SetCategoryMorphism, inv, product, coproduct,
        features, is_semisimple, is_abelian, is_monoidal, ×, ∐, RepresentationCategory, is_linear,
        GroupRepresentationCategory, is_multiring, is_multifusion, is_ring, is_multitensor, is_additive,
        is_tensor, is_fusion, is_multifusion, is_abelian, is_ring,
         FusionCategory, VSCategoryHomSpace,CategoryHomSpace,
        Hom, GVSCategoryHomSpace, HomFunctor, VSCategoryObject, GVSCategoryObject, GVSCategoryMorphism, SetCategoryHomSet,
        CategoryHomSet, Cocycle, trivial_3_cocylce,*,
        +,-, GroupRepresentation, GroupRepresentationCategoryMorphism,
        GroupRepresentationCategory,
        isinvertible, Representation, isequivariant, matrix, GRCategoryHomSpace,
        OppositeCategory, OppositeCategoryMorphism, OppositeCategoryObject, ProductCategory,
        ProductCategoryObject, ProductCategoryMorphism, CohSheaves, CohSheafCategoryObject, CohSheafCategoryMorphism,
        stalks, PullbackFunctor, Pullback, PushforwardFunctor, Pushforward,
        CohSfCategoryHomSpace, ConvolutionCategory, ConvolutionCategoryObject, ConvolutionCategoryMorphism,
        ConvCategoryHomSpace,stalk, induction, restriction, orbit_index, direct_sum,
        decompose, multiplication_table, print_multiplication_table, grothendieck_ring,
        dual, left_dual, right_dual, ev, coev, left_trace, right_trace, tr, braiding,
        drinfeld_morphism, smatrix, tmatrix, End, CenterCategory, CenterCategoryObject, CenterCategoryMorphism,
        spherical, is_central, center_simples, RingCategory, set_tensor_product!,
        set_braiding!, set_one!, Ising, zero_morphism, express_in_basis, solve_groebner,
        Center, CenterCategory, CenterCategoryObject, CenterCategoryMorphism, ev_coev, matrices,
        orbit_stabilizers, GRepRestriction, GRepInduction, Restriction, Induction,
        print_multiplication_table, RingCategoryObject, RingCatMorphism, kernel, cokernel,
        image, isgraded, cyclic_group_3cocycle, decompose_morphism,
        central_objects, half_braiding, half_braidings, left_inverse, right_inverse,
        add_simple!, pentagon_axiom, set_associator!, dual_basis, pairing,
        eigenspace, eigenspaces, simple_subobjects, I2, I2subcategory,
        DeligneProdCategoryMorphism, DeligneProdCategoryObject, DeligneProduct, ⊠, op, AbstractCategoryHomSpace,
        is_half_braiding, object, distribute_right, distribute_left, is_simple,
        decompose_morphism, TambaraYamagami, RingSubcategory, SubcategoryMorphism,
        SubcategoryCategoryObject,load,save, cat_fr_8122, E6subfactor, fpdim, set_cannonical_spherical!, cat_fr_9143,
        normalized_smatrix, sort_simples_by_dimension!, set_spherical!, inv_associator, 
        is_modular, is_spherical, TwistedGradedVectorSpaces, six_j_symbols, simples_names, HaagerupH3, factor, roots, is_subobject, rational_lift, Fibonacci, object_type



include("structures/abstracts.jl")
include("structures/FrameworkChecks.jl")
include("Utility/Technicallities.jl")
include("Utility/QQBar_Polynomials.jl")
include("Utility/SolveGroebner.jl")
include("structures/VectorSpaces/VectorSpaces.jl")
include("structures/MISC/ProductCategory.jl")
include("structures/MISC/OppositeCategory.jl")
include("structures/MISC/DeligneTensorProduct.jl")
include("structures/VectorSpaces/Cocycles.jl")
include("structures/VectorSpaces/GradedVectorSpaces.jl")
include("structures/VectorSpaces/Unitary-3-cocycle.jl")
include("structures/set.jl")
include("structures/Representations/Representations.jl")
include("structures/Representations/GroupRepresentations.jl")
include("structures/Functors.jl")
include("structures/ConvolutionCategory/CoherentSheaves.jl")
include("structures/ConvolutionCategory/ConvolutionCategory.jl")
include("structures/MultiFusionCategories/FusionCategory.jl")
include("structures/MultiFusionCategories/Skeletization.jl")
include("structures/MISC/multiplication_table.jl")
include("structures/MISC/PentagonAxiom.jl")
include("structures/MISC/Subcategories.jl")
include("structures/GrothendieckRing.jl")
include("structures/Center/Center.jl")
include("structures/Center/Induction.jl")
include("structures/Center/CenterChecks.jl")
include("structures/Examples/I2-fusion.jl")
include("structures/Examples/TambaraYamagami.jl")
include("structures/Examples/VercleyenSingerland/FR_8211/fr_8211.jl")
include("structures/Examples/VercleyenSingerland/FR_9143/fr_9143.jl")
include("structures/Examples/E6subfactor.jl")
include("structures/Examples/Haagerup/Haagerup.jl")

#include("Utility/Serialization.jl")

@alias Centre Center

end
