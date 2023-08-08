
abstract type RepresentationCategory <: Category end

abstract type RepresentationObject <: Object end

abstract type RepresentationMorphism <: VectorSpaceMorphism end


intdim(ρ::RepresentationObject) = ρ.intdim
dim(ρ::RepresentationObject) = base_ring(ρ)(intdim(ρ))
base_ring(ρ::RepresentationObject) = parent(ρ).base_ring
base_ring(Rep::RepresentationCategory) = Rep.base_ring
