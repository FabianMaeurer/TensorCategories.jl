
abstract type RepresentationCategory <: Category end

abstract type RepresentationCategoryObject <: CategoryObject end

abstract type RepresentationCategoryMorphism <: VectorSpaceCategoryMorphism end


intdim(ρ::RepresentationCategoryObject) = ρ.intdim
dim(ρ::RepresentationCategoryObject) = base_ring(ρ)(intdim(ρ))
base_ring(ρ::RepresentationCategoryObject) = parent(ρ).base_ring
base_ring(Rep::RepresentationCategory) = Rep.base_ring
