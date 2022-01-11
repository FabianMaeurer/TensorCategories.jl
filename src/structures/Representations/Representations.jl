abstract type RepresentationCategory{T} <: TensorCategory{T} end
#abstract type AlgebraRepresentationCategory{T} <: RepresentationCategory{T} end
#abstract type HopfAlgebraRepresentationCategory <: AlgebraRepresentationCategory end

abstract type Representation{T} <: Object end
# #abstract type GroupRepresentation <: Representation end
# abstract type AlgebraRepresentation <: Representation end
# abstract type HopfAlgebraRepresentation <: Representation end

abstract type RepresentationMorphism{T} <: Morphism end


dim(ρ::Representation) = ρ.dim
base_ring(ρ::Representation) = ρ.base_ring
base_ring(Rep::RepresentationCategory) = Rep.base_ring

RepresentationCategory(G::GAPGroup, F::Field) = GroupRepresentationCategory(G,F)
RepresentationCategory(G::GAPGroup) = GroupRepresentationCategory(G)

Representation(G::GAPGroup, args...) = GroupRepresentation(G,args...)
