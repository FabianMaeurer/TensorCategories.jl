module TensorCategories

import Base: *, +, -, ==, ^, getindex, getproperty, in, issubset, iterate, length, show,div, rand, split

import Oscar.AbstractAlgebra.Generic: Poly
import Oscar.Hecke: RelSimpleNumField, regular_module
import Oscar: +, @alias, @attributes, AbstractSet, AcbField, StructureConstantAlgebra, AssociativeAlgebraElem,
    cyclotomic_field, Fac, Field, FieldElem, FinField, GF, GAP, GAPGroup,
    GAPGroupHomomorphism, GL, GSet, GroupElem, Hecke.AbstractAssociativeAlgebra,
    Hecke.AbstractAssociativeAlgebraElem, Ideal, MPolyRingElem, MPolyIdeal, Map, MatElem, MatrixElem,
    MatrixGroup, matrix_space, ModuleIsomorphism, number_field, PcGroup, PolyRingElem,
    polynomial_ring, QQ, QQBar, QQField, QQFieldElem, QQMPolyRingElem, Ring, RingElem, ZZ, QQBarFieldElem,
    ZZRingElem, abelian_closure, abelian_group, absolute_simple_field, action, base_field,
    base_ring, basis, central_primitive_idempotents, change_base_ring, characteristic,
    charpoly, codomain, coeff, coefficients, cokernel, complex_embeddings, compose, centralizer,
    cyclotomic_field, decompose, degree, det, diagonal_matrix, dim, direct_sum, divisors,
    domain, dual, eigenspace, eigenspaces, eigenvalues, elem_type, elements, exponent,
    exponents, factor, QQFieldElem, QQPolyRingElem, ZZRingElem, gcd, gen, gens, get_attribute, get_attribute!,
    gmodule, groebner_basis, group_algebra, gset, guess, has_attribute, height_bits, hnf,
    hom, id, ideal, identity_matrix, image, index, inv, involution, irreducible_modules,
    is_abelian, is_central, is_finite, is_invertible, is_isomorphic, is_modular,
    is_rational, is_semisimple, is_simple, is_square, is_subfield, is_subgroup,
    is_independent, is_invertible, iso_oscar_gap,
    jordan_normal_form, kernel, kronecker_product, lcm, leading_coefficient,
    leading_monomial, left_transversal, lex, load, matrix, matrix_algebra, minpoly,
    monomials, multiplication_table, multiplicity, nullspace, nvars, one, orbit,
    orbits, order, parent, permutation_matrix, preimage, primary_decomposition, product, rank, real_solutions, resultant, root_of_unity, roots, rref, save,
    set_attribute!, size, solve, sparse_matrix, splitting_field,
    stabilizer, sub, subst, symbols, symmetric_group, tensor_power, tensor_product, tr,
    trivial_subgroup, unit, zero, zero_matrix, ∘, ⊕, ⊗, AbsSimpleNumField,
    number_of_rows, number_of_columns, is_squarefree, is_commutative,
    gens, center, graph_from_adjacency_matrix, connected_components, weakly_connected_components, Directed, Undirected, morphism, algebra,
    radical, is_zero, minimal_submodules, representation_matrix, QQBarField,
    is_irreducible, polynomial, is_univariate, action, is_equivalent, extension_of_scalars



using InteractiveUtils
using Memoization
using SparseArrays
using Base.Threads


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
export ⋆
export AbstractHomSpace 
export action
export action_matrix
export add_simple! 
export algebra
export algebras
export algebra_structures
export AlgebraObject
export AlgebraMorphism
export ArrowCategory
export ArrowObject
export ArrowMorphism
export associator 
export base_ring 
export base_group
export basis 
export bimodule
export BiModule
export BiModuleCategory
export BiModuleObject
export braiding 
export cat_fr_8122 
export cat_fr_9143 
export Category 
export category
export category_of_right_modules
export category_of_left_modules
export category_of_bimodules
export commutative_algebra_structures
export commutative_algebras
export HomSet 
export HomSpace 
export Morphism 
export Object 
export center 
export center_simples 
export CenterCategory 
export CenterCategory 
export CenterMorphism 
export CenterMorphism 
export CenterObject 
export CenterObject 
export central_objects 
export central_primitive_idempotents
export centralizer
export Cocycle 
export codomain 
export coefficients
export coequilizer
export coev 
export CohSfHomSpace 
export CohSheafMorphism 
export CohSheafObject 
export CohSheaves 
export cokernel 
export compose 
export ConvHomSpace 
export ConvolutionCategory 
export ConvolutionMorphism 
export ConvolutionObject 
export coproduct 
export cyclic_group_3cocycle 
export decompose 
export direct_sum_decomposition
export DeligneProdMorphism 
export DeligneProdObject 
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
export eigenvalues 
export End 
export end_of_free_bimodule
export end_of_free_left_module
export end_of_free_right_module
export end_of_induction
export endomorphism_ring 
export equilizer
export etale_algebra_structures
export etale_algebras
export ev 
export ev_coev 
export exponent
export express_in_basis 
export extension_of_scalars
export factor 
export free_bimodule
export free_left_module
export free_right_module
export Fibonacci 
export Forgetful 
export fpdim 
export Functor 
export FusionCategory 
export fusion_coefficient
export generic_algebra
export getindex 
export graded_vector_spaces
export GradedVectorSpaces 
export GRHomSpace 
export GRepInduction 
export GRepRestriction 
export split_grothendieck_ring 
export GrothendieckRing
export group_algebra
export GroupRepresentationCategory
export GroupRepresentation 
export GroupRepresentationCategory 
export GroupRepresentationCategory 
export GroupRepresentationMorphism 
export GVSHomSpace 
export GVSMorphism 
export GVSObject 
export haagerup_H3 
export half_braiding 
export half_braidings 
export hexagon_axiom
export hom 
export Hom 
export HomFunctor 
export horizontal_direct_sum
export I2 
export I2subcategory 
export id 
export image 
export indecomposable_subobjects
export indecomposables
export induction 
export Induction 
export induction_adjunction
export induction_monad
export induction_restriction
export induction_right_adjunction
export InductionMonad
export int_dim 
export internal_hom
export internal_hom_adjunction
export inv 
export inv_associator 
export inverse_induction_adjunction
export inverse_internal_hom_adjunction
export involution
export is_abelian 
export is_abelian 
export is_additive 
export is_algebra
export is_bimodule
export is_braided
export is_central
export is_epimorphism
export is_equivalent
export is_finite 
export is_fusion 
export is_half_braiding
export is_irreducible
export is_left_module 
export is_linear 
export is_modular 
export is_monoidal 
export is_monomorphism
export is_morita_equivalent
export is_multifusion 
export is_multifusion 
export is_multiring 
export is_multitensor 
export is_relative_braiding
export is_right_module
export is_ring 
export is_ring 
export is_semisimple 
export is_separable
export is_simple 
export is_spherical 
export is_split_semisimple
export is_subobject 
export is_tensor 
export is_zero
export isequivariant 
export isgraded 
export Ising 
export is_invertible
export karoubian_envelope 
export kernel 
export left_action
export left_dual 
export left_inverse 
export LeftModule
export LeftModuleCategory
export LeftModuleObject
export LeftTensorProductFunctor
export left_module
export left_module_category
export left_trace 
export load 
export matrices 
export matrix 
export ModuleCategory
export ModuleMorphism
export ModuleObject
export Monad
export MonadModule
export MonadModuleMorphism
export MonadModules
export Morphism, morphism 
export multiplication
export multiplication_table 
export multiplicity
export NaturalTransformation
export Nat
export normalized_smatrix 
export object 
export object_type 
export one 
export op 
export opposite_category
export opposite_morphism
export opposite_object
export orbit_index 
export orbit_stabilizers 
export pairing 
export parent 
export pentagon_axiom 
export pentagon_equations
export print_multiplication_table 
export print_multiplication_table 
export product 
export product_category 
export product_morphism
export product_object
export pullback
export Pullback 
export PullbackFunctor 
export Pushforward 
export PushforwardFunctor 
export pushout
export pushout_product
export QQBar
export QuantumZZRing
export QuantumZZRingElem
export QZZ
export radical
export rand
export rational_lift 
export regular_representation
export Representation 
export RepresentationCategory 
export representation_category
export restriction 
export Restriction 
export reverse_braiding
export right_action
export right_dual 
export right_inverse
export right_module 
export RightModule
export RightModuleCategory
export RightModuleObject
export RightTensorProductFunctor
export right_trace 
export RingCatMorphism 
export RingSubcategory 
export roots 
export save 
export Semisimplification
export SemisimplifiedObject
export SemisimplifiedMorphism
export semisimplify
export SixJCategory 
export SixJObject 
export separable_algebra_structures
export set_associator! 
export set_braiding! 
export set_canonical_spherical! 
export set_one! 
export set_spherical! 
export set_tensor_product! 
export set_trivial_spherical!
export SetHomSet 
export SetMorphism 
export SetObject 
export Sets 
export ShortExactSequence
export ShortExactSequences
export ShortExactSequenceMorphism
export simple_subobjects 
export simples 
export simples_names 
export six_j_symbols 
export smatrix 
export solve_groebner 
export sort_simples_by_dimension! 
export spherical 
export split
export stalk 
export stalks 
export SubcategoryObject 
export SubcategoryMorphism 
export tambara_yamagami
export tensor_product 
export tensor_power
export tensor_power_category
export TensorPowerCategory
export TensorPowerMorphism
export TensorPowerObject
export tmatrix 
export tr 
export trivial_3_cocylce 
export twist
export twisted_graded_vector_spaces
export TwistedGradedVectorSpaces 
export twisted_group_algebra
export unit
export unitary_cocycle
export sl2_representations
export VectorSpaceMorphism 
export VectorSpaceObject 
export VectorSpaces 
export verlinde_category
export vertical_direct_sum
export VSHomSpace 
export VSObject 
export zero 
export zero_morphism 
export ZPlusRing, ℕRing, ℤ₊Ring
export ZPlusRingElem, ℕRingElem, ℤ₊RingElem


include("CategoryFramework/AbstractTypes.jl")
include("CategoryFramework/AbstractMethods.jl")
include("CategoryFramework/FrameworkChecks.jl")
include("CategoryFramework/ProductCategory.jl")
include("CategoryFramework/Fallbacks.jl")
include("CategoryFramework/LimitsAndDiagrams/Limits.jl")
include("CategoryFramework/OppositeCategory.jl")
include("CategoryFramework/DeligneTensorProduct.jl")
include("CategoryFramework/Semisimplification.jl")
include("CategoryFramework/ArrowCategory.jl")
include("CategoryFramework/ChainComplexes/ChainComplexes.jl")
include("CategoryFramework/ChainComplexes/ShortExactSequences.jl")


include("Utility/QQBar_Polynomials.jl")
include("Utility/SolveGroebner.jl")
include("Utility/QuantumIntegers.jl")
include("Utility/Technicallities.jl")


include("Examples/GradedVectorSpaces/VectorSpaces.jl")
include("Examples/GradedVectorSpaces/Cocycles.jl")
include("Examples/GradedVectorSpaces/GradedVectorSpaces.jl")
include("Examples/GradedVectorSpaces/Unitary-3-cocycle.jl")
include("Examples/Sets/Sets.jl")
include("Examples/GroupRepresentations/Representations.jl")
include("Examples/GroupRepresentations/GroupRepresentations.jl")
include("Examples/ConvolutionCategory/CoherentSheaves.jl")
include("Examples/ConvolutionCategory/ConvolutionCategory.jl")
include("CategoryFramework/NaturalTransformations.jl")


include("TensorCategoryFramework/AbstractTensorMethods.jl")
include("TensorCategoryFramework/FusionCategory.jl")
#include("structures/MultiFusionCategories/FusionCategoryExperimental.jl")
include("TensorCategoryFramework/6j-Solver.jl")
include("TensorCategoryFramework/Skeletization.jl")
include("TensorCategoryFramework/PentagonAxiom.jl")
include("TensorCategoryFramework/HexagonAxion.jl")
include("TensorCategoryFramework/RingSubcategories.jl")
include("TensorCategoryFramework/TensorPowerCategory.jl")
include("TensorCategoryFramework/TensorFunctors.jl")
include("TensorCategoryFramework/Center/Center.jl")
include("TensorCategoryFramework/Center/Induction.jl")
include("TensorCategoryFramework/Center/InductionMonad.jl")
include("TensorCategoryFramework/Center/Centralizer.jl")
include("TensorCategoryFramework/Center/CentralizerInduction.jl")
include("TensorCategoryFramework/Center/CenterChecks.jl")


include("TensorCategoryFramework/InternalModules/InternalAlgebras.jl")
include("TensorCategoryFramework/InternalModules/ModuleCategories.jl")
include("TensorCategoryFramework/InternalModules/ComputationOfAlgebras.jl")
include("TensorCategoryFramework/ModuleCategories/MonadModules.jl")

include("DecategorifiedFramework/multiplication_table.jl")
include("DecategorifiedFramework/ZPlusRings.jl")
include("DecategorifiedFramework/GrothendieckRing.jl")


include("Examples/Verlinde/I2-fusion.jl")
include("Examples/Verlinde/Verlinde.jl")
include("Examples/UqSL2Representations/RepresentationsSL2.jl")
include("Examples/TambaraYamagami/TambaraYamagami.jl")
include("Examples/Fibonacci/FibonacciCategory.jl")
include("Examples/VercleyenSingerland/FR_8211/fr_8211.jl")
include("Examples/VercleyenSingerland/FR_9143/fr_9143.jl")
include("Examples/E6Subfactor/E6subfactor.jl")
include("Examples/Haagerup/HaagerupH3.jl")




#include("Utility/Serialization.jl")


end
