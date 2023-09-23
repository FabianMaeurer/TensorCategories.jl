
abstract type RepresentationCategory <: Category end

abstract type RepresentationObject <: Object end

abstract type RepresentationMorphism <: Morphism end


int_dim(ρ::RepresentationObject) = ρ.int_dim
dim(ρ::RepresentationObject) = base_ring(ρ)(int_dim(ρ))
base_ring(ρ::RepresentationObject) = parent(ρ).base_ring
base_ring(Rep::RepresentationCategory) = Rep.base_ring
