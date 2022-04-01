
abstract type RepresentationCategory <: Category end
#abstract type AlgebraRepresentationCategory{T} <: RepresentationCategory{T} end
#abstract type HopfAlgebraRepresentationCategory <: AlgebraRepresentationCategory end

abstract type Representation <: Object end
# #abstract type GroupRepresentation <: Representation end
# abstract type AlgebraRepresentation <: Representation end
# abstract type HopfAlgebraRepresentation <: Representation end

abstract type RepresentationMorphism <: Morphism end


dim(ρ::Representation) = ρ.dim
base_ring(ρ::Representation) = parent(ρ).base_ring
base_ring(Rep::RepresentationCategory) = Rep.base_ring
