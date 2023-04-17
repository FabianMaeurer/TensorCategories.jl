module TensorCategories

import Base: *, +, -, ==, ^, getindex, getproperty, in, issubset, iterate, length, show,div

import Oscar: +, AbstractSet, AlgAss, AlgAssElem, ComplexField, Field, FieldElem, FinField, GroupsCore,
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
    height_bits, lcm, change_base_ring, guess, direct_sum, matrix_algebra,
    @attributes, Hecke.AbsAlgAss, Hecke.AbsAlgAssElem, has_attribute



using Memoization, InteractiveUtils

export - 
export * 
export ^ 
export ∐ 
export + 
export × 
export == 
export ∘ 
export ⊕ 
export ⊗ 
export ⊠ 
export AbstractCategoryHomSpace 
export add_simple! 
export associator 
export base_ring 
export basis 
export braiding 
export cat_fr_8122 
export cat_fr_9143 
export Category 
export CategoryHomSet 
export CategoryHomSpace 
export CategoryMorphism 
export CategoryObject 
export Center 
export center_simples 
export CenterCategory 
export CenterCategory 
export CenterCategoryMorphism 
export CenterCategoryMorphism 
export CenterCategoryObject 
export CenterCategoryObject 
export central_objects 
export Cocycle 
export codomain 
export coev 
export CohSfCategoryHomSpace 
export CohSheafCategoryMorphism 
export CohSheafCategoryObject 
export CohSheaves 
export cokernel 
export compose 
export ConvCategoryHomSpace 
export ConvolutionCategory 
export ConvolutionCategoryMorphism 
export ConvolutionCategoryObject 
export coproduct 
export cyclic_group_3cocycle 
export decompose 
export decompose_morphism 
export decompose_morphism 
export DeligneProdCategoryMorphism 
export DeligneProdCategoryObject 
export DeligneProduct 
export dim 
export direct_sum 
export direct_sum 
export distribute_left 
export distribute_right 
export domain 
export drinfeld_morphism 
export dual 
export dual_basis 
export E6subfactor 
export eigenspace 
export eigenspaces 
export End 
export endomorphism_ring 
export ev 
export ev_coev 
export express_in_basis 
export factor 
export Fibonacci 
export Forgetful 
export fpdim 
export Functor 
export FusionCategory 
export getindex 
export GradedVectorSpaces 
export GRCategoryHomSpace 
export GRepInduction 
export GRepRestriction 
export grothendieck_ring 
export GrothendieckRing
export GroupRepresentation 
export GroupRepresentationCategory 
export GroupRepresentationCategory 
export GroupRepresentationCategoryMorphism 
export GVSCategoryHomSpace 
export GVSCategoryMorphism 
export GVSCategoryObject 
export HaagerupH3 
export half_braiding 
export half_braidings 
export hom 
export Hom 
export HomFunctor 
export I2 
export I2subcategory 
export id 
export image 
export indecomposable_subobjects
export indecomposables
export induction 
export Induction 
export int_dim 
export inv 
export inv_associator 
export involution
export is_abelian 
export is_abelian 
export is_additive 
export is_braided
export is_central 
export is_fusion 
export is_half_braiding 
export is_linear 
export is_modular 
export is_monoidal 
export is_multifusion 
export is_multifusion 
export is_multiring 
export is_multitensor 
export is_ring 
export is_ring 
export is_semisimple 
export is_simple 
export is_spherical 
export is_subobject 
export is_tensor 
export isequivariant 
export isgraded 
export Ising 
export isinvertible 
export kernel 
export left_dual 
export left_inverse 
export left_trace 
export load 
export matrices 
export matrix 
export Morphism 
export multiplication_table 
export normalized_smatrix 
export object 
export object_type 
export one 
export op 
export OppositeCategory 
export OppositeCategoryMorphism 
export OppositeCategoryObject 
export orbit_index 
export orbit_stabilizers 
export pairing 
export parent 
export pentagon_axiom 
export print_multiplication_table 
export print_multiplication_table 
export product 
export ProductCategory 
export ProductCategoryMorphism 
export ProductCategoryObject 
export Pullback 
export PullbackFunctor 
export Pushforward 
export PushforwardFunctor 
export rational_lift 
export Representation 
export RepresentationCategory 
export restriction 
export Restriction 
export right_dual 
export right_inverse 
export right_trace 
export RingCategory 
export RingCategoryObject 
export RingCatMorphism 
export RingSubcategory 
export roots 
export save 
export set_associator! 
export set_braiding! 
export set_canonical_spherical! 
export set_one! 
export set_spherical! 
export set_tensor_product! 
export SetCategoryHomSet 
export SetCategoryMorphism 
export SetCategoryObject 
export Sets 
export simple_subobjects 
export simples 
export simples_names 
export six_j_symbols 
export smatrix 
export solve_groebner 
export sort_simples_by_dimension! 
export spherical 
export stalk 
export stalks 
export SubcategoryCategoryObject 
export SubcategoryMorphism 
export TambaraYamagami 
export tensor_product 
export TensorCategory 
export tmatrix 
export tr 
export trivial_3_cocylce 
export TwistedGradedVectorSpaces 
export VectorSpaceCategoryMorphism 
export VectorSpaceCategoryObject 
export VectorSpaces 
export VSCategoryHomSpace 
export VSCategoryObject 
export zero 
export zero_morphism 
export ZPlusRing, ℕRing, ℤ₊Ring
export ZPlusRingElem, ℕRingElem, ℤ₊RingElem



include("structures/abstracts.jl")
include("structures/FrameworkChecks.jl")
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
include("structures/MISC/Fusionrings.jl")
include("structures/GrothendieckRing.jl")
include("structures/Center/Center.jl")
include("structures/Center/Induction.jl")
include("structures/Center/CenterChecks.jl")
include("structures/Examples/I2-fusion.jl")
include("structures/Examples/TambaraYamagami.jl")
include("structures/Examples/FibonacciCategory.jl")
include("structures/Examples/VercleyenSingerland/FR_8211/fr_8211.jl")
include("structures/Examples/VercleyenSingerland/FR_9143/fr_9143.jl")
include("structures/Examples/E6subfactor.jl")
include("structures/Examples/Haagerup/Haagerup.jl")
include("Utility/Technicallities.jl")

#include("Utility/Serialization.jl")

@alias Centre Center

end
