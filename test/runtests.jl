using TensorCategories
using Oscar,Test


@testset "Test Examples" begin
    include("VectorSpacesTest/VSTest.jl")
    include("GroupRepresentationTests/GroupRepresentationTests.jl")
    include("SixJCategoryTests/Examples.jl")
    include("UqSl2.jl")
end

@testset "Test Center/Centralizer" begin
    include("CenterTests/InductionTest.jl")
    include("CenterTests/RepCenterTest.jl")
    include("CentralizerTests/CentralizerVec.jl")
end

@testset "Test generic structures" begin
    include("SixJCategoryTests/RingCatTests.jl")
end


@testset "Test Module Categories" begin
    include("ModuleCategoryTests/ModulesTest.jl")
    include("ModuleCategoryTests/AlgebraTests.jl")
    include("ModuleCategoryTests/Non-semisimpleModules.jl")
end

include("Anyonwiki/AnyonwikiTest.jl")

#include("CoherentSheaves/ConvolutionCategoryTests.jl")

