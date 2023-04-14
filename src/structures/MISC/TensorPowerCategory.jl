mutable struct TensorPowerCategory{T} where T <: CategoryObject <: Category
    generator::T
    simples::Vector{T}

    