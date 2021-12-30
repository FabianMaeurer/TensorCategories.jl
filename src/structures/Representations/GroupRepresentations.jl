struct GroupRepresentationCategory{T,G} <: RepresentationCategory{T}
    group::GAPGroup
    base_ring::Field
end

struct GroupRepresentation{T,G} <: Representation{T}
    group::G
    m::GAPGroupHomomorphism{G,MatrixGroup{T,M}} where M <: MatrixElem
    base_ring::Ring
    dim::Int64
end

struct GroupRepresentationMorphism{T,G} <: RepresentationMorphism{T}
    domain::GroupRepresentation{T,G}
    codomain::GroupRepresentation{T,G}
    map::MatElem{T}
end

#-------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------

function GroupRepresentationCategory(G::GAPGroup, F::Field)
    return GroupRepresentationCategory{elem_type(F),typeof(G)}(G,F)
end

function GroupRepresentation(G::GAPGroup, pre_img::Vector, img::Vector)
    F = base_ring(img[1])
    d = size(img[1])[1]
    H = GL(d, F)
    m = hom(G,H, pre_img,H.(img))
    return GroupRepresentation{elem_type(F),typeof(G)}(G,m,F,d)
end



function GroupRepresentation(G::GAPGroup, m::Function)
    F = base_ring(parent(m(G[1])))
    d = size(m(G[1]))[1]
    H = GL(d,F)
    m = hom(G,H,g -> H(m(g)))
    return GroupRepresentation(G,m,F,d)
end


function Morphism(ρ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}, m::MatElem{T}) where {T,G}
    if size(m) != (dim(ρ), dim(τ)) throw(ErrorException("Mismatching dimensions")) end
    if !isequivariant(m,ρ,τ) throw(ErrorException("Map has to be equivariant")) end
    return GroupRepresentationMorphism{T,G}(ρ,τ,m)
end

#-------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------

(ρ::GroupRepresentation)(x) = ρ.m(x)

matrix(f::GroupRepresentationMorphism) = f.map

base_group(Rep::GroupRepresentationCategory) = Rep.group
base_group(ρ::GroupRepresentation) = ρ.group

parent(ρ::GroupRepresentation{T,G}) where {T,G} = GroupRepresentationCategory{T,G}(base_group(ρ), base_ring(ρ))

# function zero(Rep::GroupRepresentationCategory{T,G}) where {T,G}
#     grp = base_group(Rep)
#     F = base_ring(Rep)
#     GroupRepresentation{T,G}(grp,hom(grp,GL(0,F), g -> one(GL(0,F))))
# end

function one(Rep::GroupRepresentationCategory{T,G}) where {T,G}
    grp = base_group(Rep)
    F = base_ring(Rep)
    GroupRepresentation(grp, gens(grp), [one(MatrixSpace(F,1,1)) for _ ∈ gens(grp)])
end

==(ρ::GroupRepresentation, τ::GroupRepresentation) = *([ρ.m(g) == τ.m(g) for g ∈ gens(base_group(ρ))]...)

#-------------------------------------------------------------------------
#   Functionality: Morphisms
#-------------------------------------------------------------------------

function compose(f::GroupRepresentationMorphism{T,G}, g::GroupRepresentationMorphism{T,G}) where {T,G}
    return GroupRepresentationMorphism{T,G}(domain(f),codomain(g), matrix(f)*matrix(g))
end

#-------------------------------------------------------------------------
#   Necessities
#-------------------------------------------------------------------------

function isequivariant(m::MatElem{T}, ρ::GroupRepresentation{T}, τ::GroupRepresentation{T}) where T
    for g ∈ gens(ρ.group)
        if matrix(ρ(g))*m != m*matrix(τ(g))
            return false
        end
    end
    return true
end


#-------------------------------------------------------------------------
#   Tensor Products
#-------------------------------------------------------------------------

function tensor_product(ρ::GroupRepresentation{T}, τ::GroupRepresentation{T}) where T
    @assert ρ.group == τ.group "Mismatching groups"
    G = ρ.group
    generators = gens(G)
    return GroupRepresentation(G, generators, [kronecker_product(matrix(ρ(g)),matrix(τ(g))) for g ∈ generators])
end

#-------------------------------------------------------------------------
#   Direct Sum
#-------------------------------------------------------------------------

function dsum(ρ::GroupRepresentation{T}, τ::GroupRepresentation{T}) where T
    @assert ρ.group == τ.group "Mismatching groups"
    G = ρ.group
    F = base_ring(ρ)

    M1 = MatrixSpace(F,dim(ρ),dim(ρ))
    M2 = MatrixSpace(F,dim(ρ),dim(τ))
    M3 = MatrixSpace(F,dim(τ),dim(ρ))
    M4 = MatrixSpace(F,dim(τ),dim(τ))

    S = Representation(G, gens(G), [[matrix(ρ(g)) zero(M2); zero(M3) matrix(τ(g))] for g ∈ gens(G)])

    incl_ρ = Morphism(ρ, S, [one(M1) zero(M2)])
    incl_τ = Morphism(τ, S, [zero(M3) one(M4)])
    proj_ρ = Morphism(S, ρ, [one(M1); zero(M3)])
    proj_τ = Morphism(S, τ, [zero(M2); one(M4)])

    return S, [incl_ρ, incl_τ], [proj_ρ, proj_τ]
end

#-------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------

function simples(Rep::GroupRepresentationCategory{T,G}) where {T,G}
    grp = base_group(Rep)
    F = base_ring(Rep)

    gap_field = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))
    gap_reps = GAP.Globals.IrreducibleRepresentations(grp.X,gap_field)

    dims = [GAP.Globals.DimensionOfMatrixGroup(GAP.Globals.Range(m)) for m ∈ gap_reps]

    oscar_reps = [GAPGroupHomomorphism(grp, GL(dims[i],F), gap_reps[i]) for i ∈ 1:length(gap_reps)]
    reps = [GroupRepresentation{T,G}(grp,m,F,d) for (m,d) ∈ zip(oscar_reps,dims)]

    return reps
end

#-------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------

struct GRHomSpace{T,G} <: HomSpace{T}
    X::GroupRepresentation{T,G}
    Y::GroupRepresentation{T,G}
    basis::Vector{GroupRepresentationMorphism{T,G}}
    parent::GroupRepresentationCategory{T,G}
end

function Hom(σ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}) where {T,G}
    grp = base_group(σ)
    F = base_ring(σ)
    gap_F = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))

    #Build the modules from σ and τ
    mats_σ = GAP.GapObj([julia_to_gap(σ(g)) for g ∈ gens(grp)])
    mats_τ = GAP.GapObj([julia_to_gap(τ(g)) for g ∈ gens(grp)])
    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
    Mτ = GAP.Globals.GModuleByMats(mats_τ, gap_F)

    # Use GAPs Meat Axe to calculate a basis
    gap_homs = GAP.Globals.MTX.BasisModuleHomomorphisms(Mσ,Mτ)

    dims_m,dims_n = dim(σ), dim(τ)
    mat_homs = [matrix(F,[F(m[i,j]) for i ∈ 1:dims_m, j ∈ 1:dims_n]) for m ∈ gap_homs]

    rep_homs = [Morphism(σ,τ,m) for m ∈ mat_homs]

    return GRHomSpace{T,G}(σ,τ, rep_homs, parent(σ))
end
#-------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------


function show(io::IO, Rep::GroupRepresentationCategory)
    println(io, """Representation Category of $(Rep.group) over $(Rep.base_ring)""")
end

function show(io::IO, ρ::GroupRepresentation)
    println(io,"$(dim(ρ))-dimensional group representation over $(base_ring(ρ)) of $(ρ.group))")
end
