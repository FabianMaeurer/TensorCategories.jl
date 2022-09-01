
abstract type RepresentationCategory <: Category end
#abstract type AlgebraRepresentationCategory{T} <: RepresentationCategory{T} end
#abstract type HopfAlgebraRepresentationCategory <: AlgebraRepresentationCategory end

abstract type Representation <: Object end
# #abstract type GroupRepresentation <: Representation end
# abstract type AlgebraRepresentation <: Representation end
# abstract type HopfAlgebraRepresentation <: Representation end

abstract type RepresentationMorphism <: VectorSpaceMorphism end


intdim(ρ::Representation) = ρ.intdim
dim(ρ::Representation) = base_ring(ρ)(intdim(ρ))
base_ring(ρ::Representation) = parent(ρ).base_ring
base_ring(Rep::RepresentationCategory) = Rep.base_ring
