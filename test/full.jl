using TensorCategories
using Oscar, Test

# This file is the authoritative full test suite. Only files reached by an
# active include below contribute coverage; commented-out includes and other
# scripts under test/ remain dormant even if they contain @testsets.

@testset "Test Examples" begin
    include("VectorSpacesTest/VSTest.jl")
    # Re-enabled after verification with Oscar 1.8.1. This was disabled for
    # Oscar 1.8.0 because of GAP MeatAxe issue #6463.
    include("GroupRepresentationTests/GroupRepresentationTests.jl")
    include("SixJCategoryTests/Examples.jl")
    include("UqSl2.jl")
end

@testset "Test Center/Centralizer" begin
    include("CenterTests/InductionTest.jl")
    # include("CenterTests/RepCenterTest.jl")
    include("CenterTests/GradedVectorSpaces.jl")
    include("CentralizerTests/CentralizerVec.jl")
end

@testset "Test generic structures" begin
    include("SixJCategoryTests/RingCatTests.jl")
    include("NaturalityTests.jl")
    include("KrullSchmidtTests.jl")
end

include("InterfaceTests.jl")
include("TensorPowerTests.jl")
include("SerializationTests/OscarSerializationTests.jl")
include("LiteratureTests.jl")

@testset "Test Module Categories" begin
    include("ModuleCategoryTests/ModulesTest.jl")
    include("ModuleCategoryTests/AlgebraTests.jl")
    include("ModuleCategoryTests/Non-semisimpleModules.jl")
end

include("GroupActionsTests/TensorActionTests.jl")
include("GroupActionsTests/EquivariantizationTests.jl")

include("Anyonwiki/QuickTest.jl")
include("Anyonwiki/AnyonwikiTest.jl")

include("CoherentSheaves/CategoryPropertyTests.jl")
#include("CoherentSheaves/ConvolutionCategoryTests.jl")

include("SixJCategoryTests/SymbolConventions.jl")
