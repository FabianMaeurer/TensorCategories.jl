abstract type RepresentationCategory <: Category end
abstract type GroupRepresentationCategory <: RepresentationCategory end
abstract type AlgebraRepresentationCategory <: GroupRepresentationCategory end
abstract type HopfAlgebraRepresentationCategory <: AlgebraRepresentationCategory end

abstract type Representation end
abstract type GroupRepresentation end
abstract type AlgebraRepresentation end
abstract type HopfAlgebraRepresentation end

abstract type Algebra <: Ring end
abstract type HopfAlgebra <: Ring end

abstract type AlgebraElem{T<:FieldElem} <: RingElem end
abstract type HopfAlgebraElem{T<:FieldElem} <: AlgebraElem{T} end

#----------------------------------------------------------------------
#   Algebras
#----------------------------------------------------------------------
