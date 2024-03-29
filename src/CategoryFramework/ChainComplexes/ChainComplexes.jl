#=----------------------------------------------------------
    A framework to construct chain complex categories
    from additive categories 
----------------------------------------------------------=#

abstract type AbstractChainComplex <: Object end
abstract type AbstractChainComplexMorphism <: Morphism end
abstract type AbstractChainComplexCategory <: Category end


struct ChainComplexes <: AbstractChainComplexCategory 
    category::Category    
end

struct ChainComplex <: AbstractChainComplex
end