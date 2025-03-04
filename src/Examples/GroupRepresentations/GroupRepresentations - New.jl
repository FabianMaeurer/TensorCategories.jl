struct GroupRepresentationCategory <: RepresentationCategory
    group::Group
    base_ring::Field
end

struct GroupRepresentation <: RepresentationObject
    parent::GroupRepresentationCategory
    group::Group
    M::Union{GModule, Nothing}
    intdim::Int
end

struct GroupRepresentationMorphism <: RepresentationMorphism
    domain::GroupRepresentation
    codomain::GroupRepresentation
    map::MatElem
end

is_tensor(::GroupRepresentationCategory) = true
is_fusion(C::GroupRepresentationCategory) = mod(order(C.group),characteristic(base_ring(C))) != 0

function Base.hash(C::GroupRepresentationCategory, h::UInt)
    hash((C.group, C.base_ring), h)
end

function Base.hash(σ::GroupRepresentation, h::UInt)
    hash((getfield(σ, s) for s ∈ fieldnames(GroupRepresentation)), h)
end
#-------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------
"""
    representation_category(F::Field, G::Group)

Category of finite dimensonal group representations of ``G``.
"""
function representation_category(F::Field, G::Group)
    return GroupRepresentationCategory(G,F)
end

function representation_category(G::Group)
    return representation_category(abelian_closure(QQ)[1],G)
end

"""
    Representation(G::Group, pre_img::Vector, img::Vector)

Group representation defined by the images of generators of G.
"""
function Representation(G::Group, img::Vector)
    F = base_ring(img[1])
    d = size(img[1])[1]
    M = gmodule(G,img)
    return GroupRepresentation(representation_category(F,G),G,M,d)
end


"""
    Representation(G::Group, m::Function)

Group representation defined by m:G -> Mat_n.
"""
function Representation(G::Group, m::Function)
    @assert order(G) > 1
    F = base_ring(parent(m(G[1])))
    d = size(m(G[1]))[1]
    M = gmodule(G, [m(g) for g ∈ gens(G)])
    return GroupRepresentation(representation_category(F,G),G,M,d)
end


"""
    morphism(ρ::GroupRepresentation, τ::GroupRepresentation, m::MatElem; check = true)

Morphism between representations defined by ``m``. If check == false equivariancy
will not be checked.
"""
function morphism(ρ::GroupRepresentation, τ::GroupRepresentation, m::MatElem; check = true)
    if size(m) != (dim(ρ), dim(τ)) throw(ErrorException("Mismatching dimensions")) end
    if check
        if !isequivariant(m,ρ,τ)
            # TODO: Fix that for left_inverse 
            #  "Map has to be equivariant" 
        end
    end
    return GroupRepresentationMorphism(ρ,τ,m)
end

#-------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------
"""
    is_semisimple(C::GroupRepresentationCategory)

Return true if C is semisimple else false.
"""
function is_semisimple(C::GroupRepresentationCategory) 
    if characteristic(base_ring(C)) == 0
        return true
    end 
    gcd(characteristic(base_ring(C)), order(base_group(C))) == 1
end

function get_index(ρ::GroupRepresentation, i)
    F = base_ring(ρ)
    if ρ.m === nothing
        return GL(0,F)(zero(matrix_space(F,0,0)))
    elseif order(ρ.group) == 1
        return zero_matrix(F, int_dim(ρ), int_dim(ρ))
    else
        return ρ.M.ac[i].matrix
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
    GroupRepresentation(Rep,grp,nothing,0)
end

"""
    one(Rep::GroupRepresentationCategory)

Return the trivial representation.
"""
function one(Rep::GroupRepresentationCategory)
    grp = base_group(Rep)
    F = base_ring(Rep)
    if order(grp) == 1 
        return Representation(grp,x -> one(matrix_space(F,1,1))) end
    Representation(grp, gens(grp), [one(matrix_space(F,1,1)) for _ ∈ gens(grp)])
end

"""
    id(ρ::GroupRepresentation)

Return the identity on ρ.
"""
function id(ρ::GroupRepresentation)
    return GroupRepresentationMorphism(ρ,ρ,one(matrix_space(base_ring(ρ),intdim(ρ),intdim(ρ))))
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

#function regular_representation(C::GroupRepresentationCategory)


"""
    is_isomorphic(σ::GroupRepresentation, τ::GroupRepresentation)

Check whether σ and τ are isomorphic. If true return the isomorphism.
"""
#=  =# function is_isomorphic(σ::GroupRepresentation, τ::GroupRepresentation)
    @assert parent(σ) == parent(τ) "Mismatching parents"

    if intdim(σ) != intdim(τ) return false, nothing end
    if intdim(σ) == 0 return true, zero_morphism(σ,τ) end

    F = base_ring(σ)
    grp = σ.group

    if order(grp) == 1 return true, morphism(σ,τ,one(matrix_space(F,intdim(σ),intdim(τ)))) end

    gap_F = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))

    #Build the modules from σ and τ
    mats_σ = GAP.GapObj([GAP.julia_to_gap(σ(g)) for g ∈ gens(grp)])
    mats_τ = GAP.GapObj([GAP.julia_to_gap(τ(g)) for g ∈ gens(grp)])


    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
    Mτ = GAP.Globals.GModuleByMats(mats_τ, gap_F)


    iso = GAP.Globals.MTX.IsomorphismModules(Mσ,Mτ)

    if iso == GAP.Globals.fail return false,nothing end

    m = matrix(F,[F(iso[i,j]) for i ∈ 1:intdim(σ), j ∈ 1:intdim(τ)])
    return true, morphism(σ,τ,m)
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
    return morphism(dom,cod,m)
end

function coev(ρ::GroupRepresentation)
    dom = one(parent(ρ))
    cod = ρ ⊗ dual(ρ)
    F = base_ring(ρ)
    m = matrix(coev(VectorSpaceObject(F,intdim(ρ))))
    return morphism(dom,cod, m)
end
#-------------------------------------------------------------------------
#   Functionality: Morphisms
#-------------------------------------------------------------------------

function compose(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"
    return GroupRepresentationMorphism(domain(f),codomain(g), matrix(f)*matrix(g))
end

inv(f::GroupRepresentationMorphism)= morphism(codomain(f), domain(f), inv(matrix(f)))

associator(σ::GroupRepresentation, τ::GroupRepresentation, ρ::GroupRepresentation) = id(σ⊗τ⊗ρ)

*(x, f::GroupRepresentationMorphism) = morphism(domain(f),codomain(f),x*f.map)

function +(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    return morphism(domain(f), codomain(f), f.map + g.map)
end

#-------------------------------------------------------------------------
#   Functionality: (Co)Kernel
#-------------------------------------------------------------------------

function kernel(f::GroupRepresentationMorphism)
    ρ = domain(f)
    G = base_group(ρ)
    F = base_ring(ρ)

    k = kernel(f.map, side = :left)
    d = number_of_rows(k)

    if d == 0
        return zero(parent(ρ)), zero_morphism(zero(parent(ρ)), ρ)
    end

    k_inv = transpose(solve(transpose(k), one(matrix_space(F,d,d))))

    generators = order(G) == 1 ? elements(G) : gens(G)

    images = [k*matrix(ρ(g))*k_inv for g ∈ generators]

    K = Representation(G, generators, images)

    return K, morphism(K,ρ,k)
end

function cokernel(f::GroupRepresentationMorphism)
    ρ = codomain(f)
    G = base_group(ρ)
    F = base_ring(ρ)
    c = kernel(f.map, side = :right)
    d = number_of_columns(c)

    if d == 0
        return zero(parent(ρ)), zero_morphism(ρ,zero(parent(ρ)))
    end

    c_inv = solve(c, one(matrix_space(F,d,d)))

    generators = order(G) == 1 ? elements(G) : gens(G)

    images = [c_inv*matrix(ρ(g))*c for g ∈ generators]
    C = Representation(G, generators, images)
    return C, morphism(ρ,C,c)
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
    return morphism(dom,codom, m)
end

function braiding(X::GroupRepresentation, Y::GroupRepresentation)
    F = base_ring(X)
    n,m = intdim(X),intdim(Y)
    map = zero(matrix_space(F,n*m,n*m))
    for i ∈ 1:n, j ∈ 1:m
        v1 = matrix(F,transpose([k == i ? 1 : 0 for k ∈ 1:n]))
        v2 = matrix(F,transpose([k == j ? 1 : 0 for k ∈ 1:m]))
        map[(j-1)*n + i, :] = kronecker_product(v1,v2)
    end
    return morphism(X⊗Y, Y⊗X, transpose(map))
end

spherical(X::GroupRepresentation) = id(X)
#-------------------------------------------------------------------------
#   Direct Sum
#-------------------------------------------------------------------------

"""
    direct_sum(ρ::GroupRepresentation, τ::GroupRepresentation, morphisms::Bool = false)

Return the direct sum of representations with the corresponding injections und projections.
"""
function direct_sum(ρ::GroupRepresentation, τ::GroupRepresentation)
    @assert ρ.group == τ.group "Mismatching groups"

    grp = ρ.group
    F = base_ring(ρ)

    if ρ.m == 0
        return τ,[GroupRepresentationMorphism(ρ,τ,zero(matrix_space(F,0,intdim(τ)))), id(τ)], [GroupRepresentationMorphism(τ,ρ,zero(matrix_space(F,intdim(τ),0))), id(τ)]
    elseif τ.m == 0
        return ρ,[id(ρ), GroupRepresentationMorphism(τ,ρ,zero(matrix_space(F,0,intdim(ρ)))), id(τ)], [id(ρ), GroupRepresentationMorphism(ρ,τ,zero(matrix_space(F,intdim(ρ),0)))]
    end

    M1 = matrix_space(F,intdim(ρ),intdim(ρ))
    M2 = matrix_space(F,intdim(ρ),intdim(τ))
    M3 = matrix_space(F,intdim(τ),intdim(ρ))
    M4 = matrix_space(F,intdim(τ),intdim(τ))

    generators = order(grp) == 1 ? elements(grp) : gens(grp)

    S = Representation(grp, generators, [[matrix(ρ(g)) zero(M2); zero(M3) matrix(τ(g))] for g ∈ generators])


    incl_ρ = morphism(ρ, S, [one(M1) zero(M2)])
    incl_τ = morphism(τ, S, [zero(M3) one(M4)])
    proj_ρ = morphism(S, ρ, [one(M1); zero(M3)])
    proj_τ = morphism(S, τ, [zero(M2); one(M4)])

    return S, [incl_ρ, incl_τ], [proj_ρ, proj_τ]
end


"""
    direct_sum(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)

Direct sum of morphisms of representations.
"""
function direct_sum(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)

    dom = domain(f)⊕domain(g)
    codom = codomain(f)⊕codomain(g)
    F = base_ring(domain(f))

    z1 = zero(matrix_space(F, intdim(domain(f)), intdim(codomain(g))))
    z2 = zero(matrix_space(F, intdim(domain(g)), intdim(codomain(f))))

    m = [matrix(f) z1; z2 matrix(g)]

    return morphism(dom,codom, m)
end

# product(X::GroupRepresentation,Y::GroupRepresentation, morphisms = false) = morphisms ? direct_sum(X,Y, true)[[1,3]] : direct_sum(X,Y)
# coproduct(X::GroupRepresentation,Y::GroupRepresentation, morphisms = false) = morphisms ? direct_sum(X,Y, true)[[1,2]] : direct_sum(X,Y)


#-------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------

"""
    simples(Rep::GroupRepresentationCategory)

Return a list of the simples objects in Rep.
"""
#=  =# function simples(Rep::GroupRepresentationCategory)
    grp = base_group(Rep)
    F = base_ring(Rep)

    if order(grp) == 1 return [one(Rep)] end

    #gap_field = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))
    gap_field = codomain(iso_oscar_gap(F))
    gap_reps = GAP.Globals.IrreducibleRepresentations(grp.X,gap_field)

    intdims = [GAP.Globals.DimensionOfMatrixGroup(GAP.Globals.Range(m)) for m ∈ gap_reps]
 
    oscar_reps = [GAPGroupHomomorphism(grp, GL(intdims[i],F), gap_reps[i]) for i ∈ 1:length(gap_reps)]
    reps = [GroupRepresentation(Rep,grp,m,F,d) for (m,d) ∈ zip(oscar_reps,intdims)]

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

    rep_homs = [morphism(σ,τ,m,check = false) for m ∈ mat_homs]

    return GRHomSpace(σ,τ, rep_homs, VectorSpaces(F))
end

function zero(H::GRHomSpace)
    dom = H.X
    codom = H.Y
    m = zero(matrix_space(base_ring(dom),intdim(dom),intdim(codom)))
    return morphism(dom,codom,m)
end

function zero_morphism(X::GroupRepresentation, Y::GroupRepresentation)
    m = zero(matrix_space(base_ring(X),intdim(X),intdim(Y)))
    return morphism(X,Y,m)
end

#-------------------------------------------------------------------------
#   Restriction and Induction Functor
#-------------------------------------------------------------------------

function restriction(ρ::GroupRepresentation, H::Group)
    b,f = issubgroup(ρ.group, H)
    RepH = representation_category(base_ring(ρ),H)
    if b == false throw(ErrorException("Not a sub")) end
    if ρ.m == 0 return zero(RepH) end
    h = hom(H,codomain(ρ.m), gens(H), [ρ(f(g)) for g ∈ gens(H)])
    return GroupRepresentation(RepH, H, h, base_ring(ρ), intdim(ρ))
end

function restriction(f::GroupRepresentationMorphism, H::Group)
    if domain(f).group == H return f end
    return morphism(restriction(domain(f),H), restriction(codomain(f),H), matrix(f))
end

function induction(ρ::GroupRepresentation, G::Group)
    H = ρ.group

    if H == G return ρ end

    if !issubgroup(G, H)[1] throw(ErrorException("Not a supergroup")) end

    if ρ.m == 0 return zero(representation_category(base_ring(ρ),G)) end

    transversal = left_transversal(G,H)

    g = order(G) == 1 ? elements(G) : gens(G)

    ji = [[findfirst(x -> g[k]*t ∈ orbit(gset(H, (y,g) -> y*g, G), x), transversal) for t ∈ transversal] for k ∈ 1:length(g)]
    g_ji = [[transversal[i] for i ∈ m] for m ∈ ji]

    hi = [[inv(g_ji[k][i])*g[k]*transversal[i] for i ∈ 1:length(transversal)] for k ∈ 1:length(g)]

    images = []
    d = intdim(ρ)
    n = length(transversal)*d
    for i ∈ 1:length(g)
        m = zero(matrix_space(base_ring(ρ), n, n))

        for j ∈ 1:length(transversal)
            m[ (ji[i][j]-1)*d+1:ji[i][j]*d, (j-1)*d+1:j*d] = matrix(ρ(hi[i][j]))
        end
        images = [images; m]
    end
    return Representation(G, g, images)
end

function induction(f::GroupRepresentationMorphism, G::Group)
    dom = induction(domain(f), G)
    codom = induction(codomain(f), G)
    return morphism(dom,codom, direct_sum([morphism(matrix(f)) for i ∈ 1:Int64(index(G,domain(f).group))]).m)
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
    express_in_basis(morphism(f.map), [morphism(g.map) for g in basis])
end
