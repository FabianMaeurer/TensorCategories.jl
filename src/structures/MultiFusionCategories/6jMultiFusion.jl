struct SixJSymbols
    # function with input i1,i2,i3,i4,j,l
    symbols::Function
end

struct MFCategory{T} <: MultiFusionCategory{T}
    base_ring::Field
    simples::Int64
    six_j::SixJSymbols
end

function MultiFusionCategory(F::Field, multiplication::Array{Int,3})
end
