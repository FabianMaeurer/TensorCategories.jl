struct GroupRepresentationCategory{T,G} <: RepresentationCategory{T}
    group::GAPGroup
    base_ring::Field
end

struct GroupRepresentation{T,G} <: Representation{T}
    group::G
    m
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
"""
    RepresentationCategory(G::GAPGroup, F::Field)

Category of finite dimensonal group representations of \\G\\.
"""
function RepresentationCategory(G::GAPGroup, F::Field)
    return GroupRepresentationCategory{elem_type(F),typeof(G)}(G,F)
end

function RepresentationCategory(G::GAPGroup)
    return RepresentationCategory(G, abelian_closure(QQ)[1])
end

"""
    Representation(G::GAPGroup, pre_img::Vector, img::Vector)

Group representation defined by the images of generators of G.
"""
function Representation(G::GAPGroup, pre_img::Vector, img::Vector)
    F = base_ring(img[1])
    d = size(img[1])[1]
    H = GL(d, F)
    m = hom(G,H, pre_img,H.(img))
    return GroupRepresentation{elem_type(F),typeof(G)}(G,m,F,d)
end


"""
    Representation(G::GAPGroup, m::Function)

Group representation defined by m:G -> Mat_n.
"""
function Representation(G::GAPGroup, m::Function)
    F = order(G) == 1 ? base_ring(parent(m(elements(G)[1]))) : base_ring(parent(m(G[1])))
    d = order(G) == 1 ? size(m(elements(G)[1]))[1] : size(m(G[1]))[1]
    H = GL(d,F)
    m = hom(G,H,g -> H(m(g)))
    return GroupRepresentation{elem_type(F),typeof(G)}(G,m,F,d)
end


"""
    Morphism(ρ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}, m::MatElem{T}; check = true) where {T,G}

Morphism between representations defined by ``m``. If check == false equivariancy
will not be checked.
"""
function Morphism(ρ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}, m::MatElem{T}; check = true) where {T,G}
    if size(m) != (dim(ρ), dim(τ)) throw(ErrorException("Mismatching dimensions")) end
    if check
        if !isequivariant(m,ρ,τ) throw(ErrorException("Map has to be equivariant")) end
    end
    return GroupRepresentationMorphism{T,G}(ρ,τ,m)
end

#-------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------
"""
    issemisimple(C::GroupRepresentationCategory)

Return true if C is semisimple else false.
"""
issemisimple(C::GroupRepresentationCategory) = gcd(characteristic(base_ring(C)), order(base_group(C))) == 1

function (ρ::GroupRepresentation)(x)
    if ρ.m == 0
        F = base_ring(ρ)
        return GL(0,F)(zero(MatrixSpace(F,0,0)))
    elseif order(ρ.group) == 1
        return one(codomain(ρ.m))
    else
        return ρ.m(x)
    end
end
matrix(f::GroupRepresentationMorphism) = f.map

base_group(Rep::GroupRepresentationCategory) = Rep.group
base_group(ρ::GroupRepresentation) = ρ.group

"""
    parent(ρ::GroupRepresentation{T,G}) where {T,G}

Return the parent representation category of ρ.
"""
parent(ρ::GroupRepresentation{T,G}) where {T,G} = GroupRepresentationCategory{T,G}(base_group(ρ), base_ring(ρ))

"""
    zero(Rep::GroupRepresentationCategory{T,G}) where {T,G}

Return the zero reprensentation.
"""
function zero(Rep::GroupRepresentationCategory{T,G}) where {T,G}
    grp = base_group(Rep)
    F = base_ring(Rep)
    GroupRepresentation{T,G}(grp,0,F,0)
end

"""
    one(Rep::GroupRepresentationCategory{T,G}) where {T,G}

Return the trivial representation.
"""
function one(Rep::GroupRepresentationCategory{T,G}) where {T,G}
    grp = base_group(Rep)
    F = base_ring(Rep)
    if order(grp) == 1 return Representation(grp,x -> one(MatrixSpace(F,1,1))) end
    Representation(grp, gens(grp), [one(MatrixSpace(F,1,1)) for _ ∈ gens(grp)])
end

"""
    id(ρ::GroupRepresentation{T,G}) where {T,G}

Return the identity on ρ.
"""
function id(ρ::GroupRepresentation{T,G}) where {T,G}
    return GroupRepresentationMorphism{T,G}(ρ,ρ,one(MatrixSpace(base_ring(ρ),dim(ρ),dim(ρ))))
end

function ==(ρ::GroupRepresentation, τ::GroupRepresentation)
    if ρ.m == 0 || τ.m == 0
        return ρ.m == 0 && τ.m == 0
    elseif order(ρ.group) == 1
        if order(τ.group) == 1
            return dim(τ) == dim(ρ)
        end
        return false
    end
    return *([ρ.m(g) == τ.m(g) for g ∈ gens(base_group(ρ))]...)
end

function ==(C::RepresentationCategory, D::RepresentationCategory)
    return C.group == D.group && C.base_ring == D.base_ring
end

function ==(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    return domain(f) == domain(g) && codomain(f) == codomain(g) && f.map == g.map
end

"""
    isisomorphic(σ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}) where {T,G}

Check whether σ and τ are isomorphic. If true return the isomorphism.
"""
function isisomorphic(σ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}) where {T,G}
    @assert parent(σ) == parent(τ) "Mismatching parents"

    if dim(σ) != dim(τ) return false, nothing end
    if dim(σ) == 0 return true, zero_morphism(σ,τ) end

    F = base_ring(σ)
    grp = σ.group

    if order(grp) == 1 return true, Morphism(σ,τ,one(MatrixSpace(F,dim(σ),dim(τ)))) end

    gap_F = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))

    #Build the modules from σ and τ
    mats_σ = GAP.GapObj([julia_to_gap(σ(g)) for g ∈ gens(grp)])
    mats_τ = GAP.GapObj([julia_to_gap(τ(g)) for g ∈ gens(grp)])


    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
    Mτ = GAP.Globals.GModuleByMats(mats_τ, gap_F)


    iso = GAP.Globals.MTX.IsomorphismModules(Mσ,Mτ)

    if iso == GAP.Globals.fail return false,nothing end

    m = matrix(F,[F(iso[i,j]) for i ∈ 1:dim(σ), j ∈ 1:dim(τ)])
    return true, Morphism(σ,τ,m)
end
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
    if dim(ρ)*dim(τ) == 0 return true end
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

"""
    tensor_product(ρ::GroupRepresentation{T}, τ::GroupRepresentation{T}) where T

Return the tensor product of representations.
"""
function tensor_product(ρ::GroupRepresentation{T}, τ::GroupRepresentation{T}) where T
    @assert ρ.group == τ.group "Mismatching groups"

    if ρ.m == 0 || τ.m == 0 return zero(parent(ρ)) end

    G = ρ.group

    if order(G) == 1
        g = elements(G)[1]
        return Representation(G, x -> kronecker_product(matrix(ρ(g)),matrix(τ(g))))
    end
    generators = gens(G)
    return Representation(G, generators, [kronecker_product(matrix(ρ(g)),matrix(τ(g))) for g ∈ generators])
end

"""
    tensor_product(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)

Return the tensor product of morphisms of representations.
"""
function tensor_product(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    dom = domain(f) ⊗ domain(g)
    codom = codomain(f) ⊗ codomain(g)

    m = kronecker_product(matrix(f),matrix(g))
    return Morphism(dom,codom, m)
end

#-------------------------------------------------------------------------
#   Direct Sum
#-------------------------------------------------------------------------

"""
    dsum(ρ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}, morphisms::Bool = false) where {T,G}

Return the direct sum of representations. If morphisms is set true inclusion and
projection morphisms are also returned.
"""
function dsum(ρ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}, morphisms::Bool = false) where {T,G}
    @assert ρ.group == τ.group "Mismatching groups"

    grp = ρ.group
    F = base_ring(ρ)

    if ρ.m == 0
        if !morphisms return τ end
        return τ,[GroupRepresentationMorphism(ρ,τ,zero(MatrixSpace(F,0,dim(τ)))), id(τ)], [GroupRepresentationMorphism(τ,ρ,zero(MatrixSpace(F,dim(τ),0))), id(τ)]
    elseif τ.m == 0
        if !morphisms return ρ end
        return ρ,[id(ρ), GroupRepresentationMorphism(τ,ρ,zero(MatrixSpace(F,0,dim(ρ)))), id(τ)], [id(ρ), GroupRepresentationMorphism(ρ,τ,zero(MatrixSpace(F,dim(ρ),0)))]
    end

    M1 = MatrixSpace(F,dim(ρ),dim(ρ))
    M2 = MatrixSpace(F,dim(ρ),dim(τ))
    M3 = MatrixSpace(F,dim(τ),dim(ρ))
    M4 = MatrixSpace(F,dim(τ),dim(τ))

    generators = order(grp) == 1 ? elements(grp) : gens(grp)

    S = Representation(grp, generators, [[matrix(ρ(g)) zero(M2); zero(M3) matrix(τ(g))] for g ∈ generators])

    if !morphisms return S end

    incl_ρ = Morphism(ρ, S, [one(M1) zero(M2)])
    incl_τ = Morphism(τ, S, [zero(M3) one(M4)])
    proj_ρ = Morphism(S, ρ, [one(M1); zero(M3)])
    proj_τ = Morphism(S, τ, [zero(M2); one(M4)])

    return S, [incl_ρ, incl_τ], [proj_ρ, proj_τ]
end

"""
    dsum(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)

Direct sum of morphisms of representations.
"""
function dsum(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    dom,_,_ = domain(f)⊕domain(g)
    codom,_,_ = codomain(f)⊕codomain(g)
    F = base_ring(domain(f))

    z1 = zero(MatrixSpace(F, dim(domain(f)), dim(codomain(g))))
    z2 = zero(MatrixSpace(F, dim(domain(g)), dim(codomain(f))))

    m = [matrix(f) z1; z2 matrix(g)]

    return Morphism(dom,codom, m)
end

product(X::GroupRepresentation,Y::GroupRepresentation, morphisms = false) = morphisms ? dsum(X,Y, true)[[1,3]] : dsum(X,Y)
coproduct(X::GroupRepresentation,Y::GroupRepresentation, morphisms = false) = morphisms ? dsum(X,Y, true)[[1,2]] : dsum(X,Y)


#-------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------

"""
    simples(Rep::GroupRepresentationCategory{T,G}) where {T,G}

Return a list of the simple objects in Rep.
"""
function simples(Rep::GroupRepresentationCategory{T,G}) where {T,G}
    grp = base_group(Rep)
    F = base_ring(Rep)

    if order(grp) == 1 return [one(Rep)] end

    if characteristic(F) == 0
        mods = irreducible_modules(grp)
        reps = [Representation(grp,gens(grp),[matrix(x) for x ∈ action(m)]) for m ∈ mods]
        return reps
    else

        gap_field = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))
        gap_reps = GAP.Globals.IrreducibleRepresentations(grp.X,gap_field)

        dims = [GAP.Globals.DimensionOfMatrixGroup(GAP.Globals.Range(m)) for m ∈ gap_reps]

        oscar_reps = [GAPGroupHomomorphism(grp, GL(dims[i],F), gap_reps[i]) for i ∈ 1:length(gap_reps)]
        reps = [GroupRepresentation{T,G}(grp,m,F,d) for (m,d) ∈ zip(oscar_reps,dims)]

        return reps
    end
end

"""
    decompose(σ::GroupRepresentation)

Decompose the representation into a direct sum of simple objects. Return a
list of tuples with simple objects and multiplicities.
"""
function decompose(σ::GroupRepresentation)
    F = base_ring(σ)
    if dim(σ) == 0 return [] end
    G = σ.group

    if order(G) == 1 return [(one(parent(σ)),dim(σ))] end

    M = to_gap_module(σ,F)
    ret = []
    facs = GAP.Globals.MTX.CollectedFactors(M)
    d = dim(σ)
    for m ∈ facs
        imgs = [matrix(F,[F(n[i,j]) for i ∈ 1:length(n), j ∈ 1:length(n)]) for n ∈ m[1].generators]
        ret = [ret;(Representation(G,gens(G),imgs),gap_to_julia(m[2]))]
    end
    ret
end

#-------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------

struct GRHomSpace{T,G} <: HomSpace{T}
    X::GroupRepresentation{T,G}
    Y::GroupRepresentation{T,G}
    basis::Vector{GroupRepresentationMorphism{T,G}}
    parent::VectorSpaces{T}
end

"""
    Hom(σ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}) where {T,G}

Return the hom-space of the representations as a vector space.
"""
function Hom(σ::GroupRepresentation{T,G}, τ::GroupRepresentation{T,G}) where {T,G}
    grp = base_group(σ)
    F = base_ring(σ)

    if dim(σ)*dim(τ) == 0 return GRHomSpace{T,G}(σ,τ,GroupRepresentationMorphism{T,G}[],VectorSpaces(F)) end

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

    rep_homs = [Morphism(σ,τ,m,check = false) for m ∈ mat_homs]

    return GRHomSpace{T,G}(σ,τ, rep_homs, VectorSpaces(F))
end

function zero(H::GRHomSpace{T,G}) where {T,G}
    dom = H.X
    codom = H.Y
    m = zero(MatrixSpace(base_ring(dom),dim(dom),dim(codom)))
    return Morphism(dom,codom,m)
end

function zero_morphism(X::GroupRepresentation, Y::GroupRepresentation)
    m = zero(MatrixSpace(base_ring(X),dim(X),dim(Y)))
    return Morphism(X,Y,m)
end

#-------------------------------------------------------------------------
#   Restriction and Induction Functor
#-------------------------------------------------------------------------

function restriction(ρ::GroupRepresentation{T,G}, H::GAPGroup) where {T,G}
    b,f = issubgroup(ρ.group, H)

    if b == false throw(ErrorException("Not a subgroup")) end
    if ρ.m == 0 return zero(RepresentationCategory(H,base_ring(ρ))) end
    h = hom(H,codomain(ρ.m), gens(H), [ρ(f(g)) for g ∈ gens(H)])
    return GroupRepresentation{T,G}(H, h, base_ring(ρ), dim(ρ))
end

function restriction(f::GroupRepresentationMorphism, H::GAPGroup)
    @show domain(f).group, H
    return Morphism(restriction(domain(f),H), restriction(codomain(f),H), matrix(f))
end

function induction(ρ::GroupRepresentation{T,S}, G::GAPGroup) where {T,S}
    H = ρ.group

    if H == G return ρ end

    if !issubgroup(G, H)[1] throw(ErrorException("Not a supergroup")) end

    if ρ.m == 0 return zero(RepresentationCategory(G,base_ring(ρ))) end

    transversal = left_transversal(G,H)

    g = order(G) == 1 ? elements(G) : gens(G)

    ji = [[findfirst(x -> g[k]*t ∈ orbit(gset(H, (y,g) -> y*g, G), x), transversal) for t ∈ transversal] for k ∈ 1:length(g)]
    g_ji = [[transversal[i] for i ∈ m] for m ∈ ji]

    hi = [[inv(g_ji[k][i])*g[k]*transversal[i] for i ∈ 1:length(transversal)] for k ∈ 1:length(g)]

    images = []
    d = dim(ρ)
    n = length(transversal)*d
    for i ∈ 1:length(g)
        m = zero(MatrixSpace(base_ring(ρ), n, n))

        for j ∈ 1:length(transversal)
            m[ (ji[i][j]-1)*d+1:ji[i][j]*d, (j-1)*d+1:j*d] = matrix(ρ(hi[i][j]))
        end
        images = [images; m]
    end
    return Representation(G, g, images)
end

function induction(f::GroupRepresentationMorphism, G::GAPGroup)
    dom = induction(domain(f), G)
    codom = induction(codomain(f), G)
    return Morphism(dom,codom, dsum([Morphism(matrix(f)) for i ∈ 1:Int64(index(G,domain(f).group))]).m)
end
#-------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------


function show(io::IO, Rep::GroupRepresentationCategory)
    print(io, """Representation Category of $(Rep.group) over $(Rep.base_ring)""")
end

function show(io::IO, ρ::GroupRepresentation)
    print(io,"$(dim(ρ))-dimensional group representation over $(base_ring(ρ)) of $(ρ.group))")
end


function to_gap_module(σ::GroupRepresentation,F::Field)
    grp = σ.group
    gap_F = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))
    mats_σ = GAP.GapObj([julia_to_gap(σ(g)) for g ∈ gens(grp)])
    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
end
