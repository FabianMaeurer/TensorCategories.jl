struct GroupRepresentationCategory <: RepresentationCategory
    group::GAPGroup
    base_ring::Field
end

struct GroupRepresentation <: Representation
    parent::GroupRepresentationCategory
    group::GAPGroup
    m
    base_ring::Ring
    intdim::Int
end

struct GroupRepresentationMorphism <: RepresentationMorphism
    domain::GroupRepresentation
    codomain::GroupRepresentation
    map::MatElem
end

istensor(::GroupRepresentationCategory) = true
isfusion(C::GroupRepresentationCategory) = mod(order(C.group),characteristic(base_ring(C))) != 0

const Group_Rep_Cache = Dict{Any,Dict{Any,Any}}()
#-------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------
"""
    RepresentationCategory(G::GAPGroup, F::Field)

Category of finite dimensonal group representations of \\G\\.
"""
function RepresentationCategory(G::GAPGroup, F::Field)
    return GroupRepresentationCategory(G,F)
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
    return GroupRepresentation(RepresentationCategory(G,F),G,m,F,d)
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
    return GroupRepresentation(RepresentationCategory(G,F),G,m,F,d)
end


"""
    Morphism(ρ::GroupRepresentation, τ::GroupRepresentation, m::MatElem; check = true)

Morphism between representations defined by ``m``. If check == false equivariancy
will not be checked.
"""
function Morphism(ρ::GroupRepresentation, τ::GroupRepresentation, m::MatElem; check = true)
    if size(m) != (dim(ρ), dim(τ)) throw(ErrorException("Mismatching dimensions")) end
    if check
        if !isequivariant(m,ρ,τ) throw(ErrorException("Map has to be equivariant")) end
    end
    return GroupRepresentationMorphism(ρ,τ,m)
end

#-------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------
"""
    issemisimple(C::GroupRepresentationCategory)

Return true if C is semisimple else false.
"""
function issemisimple(C::GroupRepresentationCategory) 
    if characteristic(base_ring(C)) == 0
        return true
    end 
    gcd(characteristic(base_ring(C)), order(base_group(C))) == 1
end

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
    parent(ρ::GroupRepresentation)

Return the parent representation category of ρ.
"""
parent(ρ::GroupRepresentation) = ρ.parent

"""
    zero(Rep::GroupRepresentationCategory)

Return the zero reprensentation.
"""
function zero(Rep::GroupRepresentationCategory)
    grp = base_group(Rep)
    F = base_ring(Rep)
    GroupRepresentation(Rep,grp,0,F,0)
end

"""
    one(Rep::GroupRepresentationCategory)

Return the trivial representation.
"""
function one(Rep::GroupRepresentationCategory)
    grp = base_group(Rep)
    F = base_ring(Rep)
    if order(grp) == 1 return Representation(grp,x -> one(MatrixSpace(F,1,1))) end
    Representation(grp, gens(grp), [one(MatrixSpace(F,1,1)) for _ ∈ gens(grp)])
end

"""
    id(ρ::GroupRepresentation)

Return the identity on ρ.
"""
function id(ρ::GroupRepresentation)
    return GroupRepresentationMorphism(ρ,ρ,one(MatrixSpace(base_ring(ρ),intdim(ρ),intdim(ρ))))
end

function ==(ρ::GroupRepresentation, τ::GroupRepresentation)
    if ρ.m == 0 || τ.m == 0
        return ρ.m == 0 && τ.m == 0
    elseif order(ρ.group) == 1
        if order(τ.group) == 1
            return intdim(τ) == intdim(ρ)
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
    isisomorphic(σ::GroupRepresentation, τ::GroupRepresentation)

Check whether σ and τ are isomorphic. If true return the isomorphism.
"""
function isisomorphic(σ::GroupRepresentation, τ::GroupRepresentation)
    @assert parent(σ) == parent(τ) "Mismatching parents"

    if intdim(σ) != intdim(τ) return false, nothing end
    if intdim(σ) == 0 return true, zero_morphism(σ,τ) end

    F = base_ring(σ)
    grp = σ.group

    if order(grp) == 1 return true, Morphism(σ,τ,one(MatrixSpace(F,intdim(σ),intdim(τ)))) end

    gap_F = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))

    #Build the modules from σ and τ
    mats_σ = GAP.GapObj([GAP.julia_to_gap(σ(g)) for g ∈ gens(grp)])
    mats_τ = GAP.GapObj([GAP.julia_to_gap(τ(g)) for g ∈ gens(grp)])


    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
    Mτ = GAP.Globals.GModuleByMats(mats_τ, gap_F)


    iso = GAP.Globals.MTX.IsomorphismModules(Mσ,Mτ)

    if iso == GAP.Globals.fail return false,nothing end

    m = matrix(F,[F(iso[i,j]) for i ∈ 1:intdim(σ), j ∈ 1:intdim(τ)])
    return true, Morphism(σ,τ,m)
end

function dual(ρ::GroupRepresentation)
    G = base_group(ρ)
    F = base_ring(ρ)
    if intdim(ρ) == 0 return ρ end
    generators = order(G) == 1 ? elements(G) : gens(G)
    return Representation(G, generators, [transpose(matrix(ρ(inv(g)))) for g ∈ generators])
end

function ev(ρ::GroupRepresentation)
    dom = dual(ρ) ⊗ ρ
    cod = one(parent(ρ))
    F = base_ring(ρ)
    m = matrix(ev(VectorSpaceObject(F,intdim(ρ))))
    return Morphism(dom,cod,m)
end

function coev(ρ::GroupRepresentation)
    dom = one(parent(ρ))
    cod = ρ ⊗ dual(ρ)
    F = base_ring(ρ)
    m = matrix(coev(VectorSpaceObject(F,intdim(ρ))))
    return Morphism(dom,cod, m)
end
#-------------------------------------------------------------------------
#   Functionality: Morphisms
#-------------------------------------------------------------------------

function compose(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"
    return GroupRepresentationMorphism(domain(f),codomain(g), matrix(f)*matrix(g))
end

associator(σ::GroupRepresentation, τ::GroupRepresentation, ρ::GroupRepresentation) = id(σ⊗τ⊗ρ)

*(x, f::GroupRepresentationMorphism) = Morphism(domain(f),codomain(f),x*f.map)

function +(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    return Morphism(domain(f), codomain(f), f.map + g.map)
end

function (F::Field)(f::GroupRepresentationMorphism)
    D = domain(f)
    C = codomain(f)
    if intdim(D) == intdim(C) == 1
        return F(f.map[1,1])
    else
        throw(ErrorException("Cannot coerce"))
    end
end

#-------------------------------------------------------------------------
#   Functionality: (Co)Kernel
#-------------------------------------------------------------------------

function kernel(f::GroupRepresentationMorphism)
    ρ = domain(f)
    G = base_group(ρ)
    F = base_ring(ρ)

    d,k = kernel(f.map, side = :left)
    k = k[1:d,:]

    if d == 0
        return zero(parent(ρ)), zero_morphism(zero(parent(ρ)), ρ)
    end

    k_inv = transpose(solve_left(transpose(k), one(MatrixSpace(F,d,d))))

    generators = order(G) == 1 ? elements(G) : gens(G)

    images = [k*matrix(ρ(g))*k_inv for g ∈ generators]

    K = Representation(G, generators, images)

    return K, Morphism(K,ρ,k)
end

function cokernel(f::GroupRepresentationMorphism)
    ρ = codomain(f)
    G = base_group(ρ)
    F = base_ring(ρ)
    d,c = kernel(f.map, side = :right)
    c = c[:,1:d]

    if d == 0
        return zero(parent(ρ)), zero_morphism(ρ,zero(parent(ρ)))
    end

    c_inv = solve_left(c, one(MatrixSpace(F,d,d)))

    generators = order(G) == 1 ? elements(G) : gens(G)

    images = [c_inv*matrix(ρ(g))*c for g ∈ generators]
    C = Representation(G, generators, images)
    return C, Morphism(ρ,C,c)
end
#-------------------------------------------------------------------------
#   Necessities
#-------------------------------------------------------------------------

function isequivariant(m::MatElem, ρ::GroupRepresentation, τ::GroupRepresentation)
    if intdim(ρ)*intdim(τ) == 0 return true end
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
    tensor_product(ρ::GroupRepresentation, τ::GroupRepresentation)

Return the tensor product of representations.
"""
function tensor_product(ρ::GroupRepresentation, τ::GroupRepresentation)
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

function braiding(X::GroupRepresentation, Y::GroupRepresentation)
    F = base_ring(X)
    n,m = intdim(X),intdim(Y)
    map = zero(MatrixSpace(F,n*m,n*m))
    for i ∈ 1:n, j ∈ 1:m
        v1 = matrix(F,transpose([k == i ? 1 : 0 for k ∈ 1:n]))
        v2 = matrix(F,transpose([k == j ? 1 : 0 for k ∈ 1:m]))
        map[(j-1)*n + i, :] = kronecker_product(v1,v2)
    end
    return Morphism(X⊗Y, Y⊗X, transpose(map))
end

spherical(X::GroupRepresentation) = id(X)
#-------------------------------------------------------------------------
#   Direct Sum
#-------------------------------------------------------------------------

"""
    dsum(ρ::GroupRepresentation, τ::GroupRepresentation, morphisms::Bool = false)

Return the direct sum of representations. If morphisms is set true inclusion and
projection morphisms are also returned.
"""
function dsum(ρ::GroupRepresentation, τ::GroupRepresentation, morphisms::Bool = false)
    @assert ρ.group == τ.group "Mismatching groups"

    grp = ρ.group
    F = base_ring(ρ)

    if ρ.m == 0
        if !morphisms return τ end
        return τ,[GroupRepresentationMorphism(ρ,τ,zero(MatrixSpace(F,0,intdim(τ)))), id(τ)], [GroupRepresentationMorphism(τ,ρ,zero(MatrixSpace(F,intdim(τ),0))), id(τ)]
    elseif τ.m == 0
        if !morphisms return ρ end
        return ρ,[id(ρ), GroupRepresentationMorphism(τ,ρ,zero(MatrixSpace(F,0,intdim(ρ)))), id(τ)], [id(ρ), GroupRepresentationMorphism(ρ,τ,zero(MatrixSpace(F,intdim(ρ),0)))]
    end

    M1 = MatrixSpace(F,intdim(ρ),intdim(ρ))
    M2 = MatrixSpace(F,intdim(ρ),intdim(τ))
    M3 = MatrixSpace(F,intdim(τ),intdim(ρ))
    M4 = MatrixSpace(F,intdim(τ),intdim(τ))

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

    dom = domain(f)⊕domain(g)
    codom = codomain(f)⊕codomain(g)
    F = base_ring(domain(f))

    z1 = zero(MatrixSpace(F, intdim(domain(f)), intdim(codomain(g))))
    z2 = zero(MatrixSpace(F, intdim(domain(g)), intdim(codomain(f))))

    m = [matrix(f) z1; z2 matrix(g)]

    return Morphism(dom,codom, m)
end

# product(X::GroupRepresentation,Y::GroupRepresentation, morphisms = false) = morphisms ? dsum(X,Y, true)[[1,3]] : dsum(X,Y)
# coproduct(X::GroupRepresentation,Y::GroupRepresentation, morphisms = false) = morphisms ? dsum(X,Y, true)[[1,2]] : dsum(X,Y)


#-------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------

"""
    simples(Rep::GroupRepresentationCategory)

Return a list of the simple objects in Rep.
"""
function simples(Rep::GroupRepresentationCategory)
    if simples ∈ keys(Group_Rep_Cache) 
        if Rep ∈ keys(Group_Rep_Cache[simples])
            return Group_Rep_Cache[simples][Rep]
        end
    else
        Group_Rep_Cache[simples] = Dict{Any,Any}()
    end

    grp = base_group(Rep)
    F = base_ring(Rep)

    if order(grp) == 1 return [one(Rep)] end

    #gap_field = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))
    gap_field = codomain(iso_oscar_gap(F))
    gap_reps = GAP.Globals.IrreducibleRepresentations(grp.X,gap_field)

    intdims = [GAP.Globals.DimensionOfMatrixGroup(GAP.Globals.Range(m)) for m ∈ gap_reps]
 
    oscar_reps = [GAPGroupHomomorphism(grp, GL(intdims[i],F), gap_reps[i]) for i ∈ 1:length(gap_reps)]
    reps = [GroupRepresentation(Rep,grp,m,F,d) for (m,d) ∈ zip(oscar_reps,intdims)]

    
    push!(Group_Rep_Cache[simples],Rep => reps)

    return reps
end

"""
    decompose(σ::GroupRepresentation)

Decompose the representation into a direct sum of simple objects. Return a
list of tuples with simple objects and multiplicities.
"""
function decompose(σ::GroupRepresentation)
    F = base_ring(σ)
    if intdim(σ) == 0 return [] end
    G = σ.group

    if order(G) == 1 return [(one(parent(σ)),intdim(σ))] end

    M = to_gap_module(σ,F)
    ret = []
    facs = GAP.Globals.MTX.CollectedFactors(M)
    d = intdim(σ)
    for m ∈ facs
        imgs = [matrix(F,[F(n[i,j]) for i ∈ 1:length(n), j ∈ 1:length(n)]) for n ∈ m[1].generators]
        ret = [ret;(Representation(G,gens(G),imgs),GAP.gap_to_julia(m[2]))]
    end
    ret
end

function is_simple(σ::GroupRepresentation)
    length(decompose(σ)) == 1
end
#-------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------

struct GRHomSpace<: AbstractHomSpace
    X::GroupRepresentation
    Y::GroupRepresentation
    basis::Vector{GroupRepresentationMorphism}
    parent::VectorSpaces
end

"""
    Hom(σ::GroupRepresentation, τ::GroupRepresentation)

Return the hom-space of the representations as a vector space.
"""
function Hom(σ::GroupRepresentation, τ::GroupRepresentation)
    grp = base_group(σ)
    F = base_ring(σ)

    if intdim(σ)*intdim(τ) == 0 return GRHomSpace(σ,τ,GroupRepresentationMorphism[],VectorSpaces(F)) end

    gap_to_F = iso_oscar_gap(F)
    gap_F = codomain(gap_to_F)
    generators = order(grp) == 1 ? elements(grp) : gens(grp)

    #Build the modules from σ and τ
    mats_σ = GAP.GapObj([GAP.julia_to_gap(σ(g)) for g ∈ generators])
    mats_τ = GAP.GapObj([GAP.julia_to_gap(τ(g)) for g ∈ generators])
    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
    Mτ = GAP.Globals.GModuleByMats(mats_τ, gap_F)

    # Use GAPs Meat Axe to calculate a basis
    gap_homs = GAP.Globals.MTX.BasisModuleHomomorphisms(Mσ,Mτ)

    intdims_m,intdims_n = intdim(σ), intdim(τ)

    mat_homs = [matrix(F,[preimage(gap_to_F, m[i,j]) for i ∈ 1:intdims_m, j ∈ 1:intdims_n]) for m ∈ gap_homs]

    rep_homs = [Morphism(σ,τ,m,check = false) for m ∈ mat_homs]

    return GRHomSpace(σ,τ, rep_homs, VectorSpaces(F))
end

function zero(H::GRHomSpace)
    dom = H.X
    codom = H.Y
    m = zero(MatrixSpace(base_ring(dom),intdim(dom),intdim(codom)))
    return Morphism(dom,codom,m)
end

function zero_morphism(X::GroupRepresentation, Y::GroupRepresentation)
    m = zero(MatrixSpace(base_ring(X),intdim(X),intdim(Y)))
    return Morphism(X,Y,m)
end

#-------------------------------------------------------------------------
#   Restriction and Induction Functor
#-------------------------------------------------------------------------

function restriction(ρ::GroupRepresentation, H::GAPGroup)
    b,f = issubgroup(ρ.group, H)
    RepH = RepresentationCategory(H,base_ring(ρ))
    if b == false throw(ErrorException("Not a subgroup")) end
    if ρ.m == 0 return zero(RepH) end
    h = hom(H,codomain(ρ.m), gens(H), [ρ(f(g)) for g ∈ gens(H)])
    return GroupRepresentation(RepH, H, h, base_ring(ρ), intdim(ρ))
end

function restriction(f::GroupRepresentationMorphism, H::GAPGroup)
    if domain(f).group == H return f end
    return Morphism(restriction(domain(f),H), restriction(codomain(f),H), matrix(f))
end

function induction(ρ::GroupRepresentation, G::GAPGroup)
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
    d = intdim(ρ)
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
    print(io,"$(intdim(ρ))-dimensional group representation over $(base_ring(ρ)) of $(ρ.group))")
end

function show(io::IO, f::GroupRepresentationMorphism)
    println(io, "Group representation Morphism with defining matrix")
    print(io,f.map)
end

#-------------------------------------------------------------------------
#   Utility
#-------------------------------------------------------------------------


function to_gap_module(σ::GroupRepresentation,F::Field)
    grp = σ.group
    gap_F = codomain(iso_oscar_gap(F))
    mats_σ = GAP.GapObj([GAP.julia_to_gap(σ(g)) for g ∈ gens(grp)])
    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
end

function express_in_basis(f::GroupRepresentationMorphism, basis::Vector{GroupRepresentationMorphism})
    o = one(base_group(domain(f)))
    express_in_basis(Morphism(f.map), [Morphism(g.map) for g in basis])
end
