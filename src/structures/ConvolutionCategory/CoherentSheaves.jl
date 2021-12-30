struct EqCoherentSheaves{T,G} <: MultiFusionCategory{T}

end

struct EqCoherentSheaf{T} <: Object
    parent::EqCoherentSheaves{T}
    stalks::Vector{VectorSpaceObject{T}}
end
