@attributes mutable struct GroupRepresentationCategory <: RepresentationCategory
    group::GAPGroup
    base_ring::Field

    function GroupRepresentationCategory(G::GAPGroup, F::Field) 
        C = new(G,F)
        if characteristic(F) !== 0 && rem(order(G), characteristic(F)) == 0 
            set_attribute!(C, :semisimple, false)
            set_attribute!(C, :tensor, true)
        else
            set_attribute!(C, :fusion, true)
        end
        C
    end
end

struct GroupRepresentation <: RepresentationObject
    parent::GroupRepresentationCategory
    group::GAPGroup
    m
    base_ring::Ring
    int_dim::Int
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

int_dim(ρ::GroupRepresentation) = ρ.int_dim
#-------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------
"""
    representation_category(F::Field, G::GAPGroup)

Category of finite dimensonal group representations of ``G``.
"""
function representation_category(F::Field, G::GAPGroup)
    return GroupRepresentationCategory(G,F)
end

function representation_category(G::GAPGroup)
    return representation_category(abelian_closure(QQ)[1], G)
end

"""
    Representation(G::GAPGroup, pre_img::Vector, img::Vector)

Group representation defined by the images of generators of G.
"""
function Representation(G::GAPGroup, pre_img::Vector, img::Vector; check::Bool = true)
    F = base_ring(img[1])
    Representation(representation_category(F, G),pre_img, img, check = check)
end

function Representation(C::GroupRepresentationCategory, pre_img::Vector, img::Vector; check::Bool = true)
    F = base_ring(img[1])
    d = size(img[1])[1]
    H = GL(d, F)
    G = base_group(C)
    m = hom(G,H, pre_img, H.(img), check = check)
    GroupRepresentation(C,G,m,F,d)
end

"""
    Representation(G::GAPGroup, m::Function)

Group representation defined by m:G -> Mat_n.
"""
function Representation(G::GAPGroup, m::Function)
    F = order(G) == 1 ? base_ring(parent(m(elements(G)[1]))) : base_ring(parent(m(G[1])))
    Representation(representation_category(F,G),m)
end

function Representation(C::GroupRepresentationCategory, m::Function)
    G = base_group(C)
    F = order(G) == 1 ? base_ring(parent(m(elements(G)[1]))) : base_ring(parent(m(G[1])))
    d = order(G) == 1 ? size(m(elements(G)[1]))[1] : size(m(G[1]))[1]
    H = GL(d,F)
    m = hom(G,H,g -> H(m(g)))
    return GroupRepresentation(C,G,m,F,d)
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

function (ρ::GroupRepresentation)(x)
    if ρ.m == 0
        F = base_ring(ρ)
        return GL(0,F)(zero(matrix_space(F,0,0)))
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
    if order(grp) == 1 return Representation(Rep,x -> one(matrix_space(F,1,1))) end
    Representation(Rep, gens(grp), [one(matrix_space(F,1,1)) for _ ∈ gens(grp)], check = false)
end

"""
    id(ρ::GroupRepresentation)

Return the identity on ρ.
"""
function id(ρ::GroupRepresentation)
    return GroupRepresentationMorphism(ρ,ρ,one(matrix_space(base_ring(ρ),int_dim(ρ),int_dim(ρ))))
end

function ==(ρ::GroupRepresentation, τ::GroupRepresentation)
    if ρ.m == 0 || τ.m == 0
        return ρ.m == 0 && τ.m == 0
    elseif order(ρ.group) == 1
        if order(τ.group) == 1
            return int_dim(τ) == int_dim(ρ)
        end
        return false
    end
    if base_group(ρ) != base_group(τ)
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
#= @memoize Dict =# function is_isomorphic(σ::GroupRepresentation, τ::GroupRepresentation)
    @assert parent(σ) == parent(τ) "Mismatching parents"

    if int_dim(σ) != int_dim(τ) return false, nothing end
    if int_dim(σ) == 0 return true, zero_morphism(σ,τ) end

    F = base_ring(σ)
    grp = σ.group

    if order(grp) == 1 return true, morphism(σ,τ,one(matrix_space(F,int_dim(σ),int_dim(τ)))) end

    gap_F = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))

    #Build the modules from σ and τ
    mats_σ = GAP.GapObj([GAP.GapObj(σ(g)) for g ∈ gens(grp)])
    mats_τ = GAP.GapObj([GAP.GapObj(τ(g)) for g ∈ gens(grp)])


    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
    Mτ = GAP.Globals.GModuleByMats(mats_τ, gap_F)


    iso = GAP.Globals.MTX.IsomorphismModules(Mσ,Mτ)

    if iso == GAP.Globals.fail return false,nothing end

    m = matrix(F,[F(iso[i,j]) for i ∈ 1:int_dim(σ), j ∈ 1:int_dim(τ)])
    return true, morphism(σ,τ,m)
end

function dual(ρ::GroupRepresentation)
    G = base_group(ρ)
    F = base_ring(ρ)
    if int_dim(ρ) == 0 return ρ end
    generators = order(G) == 1 ? elements(G) : gens(G)
    return Representation(parent(ρ), generators, [transpose(matrix(ρ(inv(g)))) for g ∈ generators], check = false)
end

function ev(ρ::GroupRepresentation)
    dom = dual(ρ) ⊗ ρ
    cod = one(parent(ρ))
    F = base_ring(ρ)
    m = matrix(ev(VectorSpaceObject(F,int_dim(ρ))))
    return morphism(dom,cod,m, check = false)
end

function coev(ρ::GroupRepresentation)
    dom = one(parent(ρ))
    cod = ρ ⊗ dual(ρ)
    F = base_ring(ρ)
    m = matrix(coev(VectorSpaceObject(F,int_dim(ρ))))
    return morphism(dom,cod, m, check = false)
end
#-------------------------------------------------------------------------
#   Functionality: Morphisms
#-------------------------------------------------------------------------

function compose(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"
    return GroupRepresentationMorphism(domain(f),codomain(g), matrix(f)*matrix(g))
end

inv(f::GroupRepresentationMorphism)= morphism(codomain(f), domain(f), inv(matrix(f)), check = false)

associator(σ::GroupRepresentation, τ::GroupRepresentation, ρ::GroupRepresentation) = id(σ⊗τ⊗ρ)

*(x, f::GroupRepresentationMorphism) = morphism(domain(f),codomain(f),x*f.map, check = false)

function +(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    return morphism(domain(f), codomain(f), f.map + g.map, check = false)
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

    K = Representation(parent(f), generators, images, check = false)

    return K, morphism(K,ρ,k, check = false)
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
    C = Representation(parent(f), generators, images, check = false)
    return C, morphism(ρ,C,c,check = false)
end
#-------------------------------------------------------------------------
#   Necessities
#-------------------------------------------------------------------------

function isequivariant(m::MatElem, ρ::GroupRepresentation, τ::GroupRepresentation)
    if int_dim(ρ)*int_dim(τ) == 0 return true end
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
        return Representation(parent(ρ), x -> kronecker_product(matrix(ρ(g)),matrix(τ(g))))
    end
    generators = gens(G)
    return Representation(parent(ρ), generators, [kronecker_product(matrix(ρ(g)),matrix(τ(g))) for g ∈ generators], check = false)
end

"""
    tensor_product(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)

Return the tensor product of morphisms of representations.
"""
function tensor_product(f::GroupRepresentationMorphism, g::GroupRepresentationMorphism)
    dom = domain(f) ⊗ domain(g)
    codom = codomain(f) ⊗ codomain(g)

    m = kronecker_product(matrix(f),matrix(g))
    return morphism(dom,codom, m, check = false)
end

function braiding(X::GroupRepresentation, Y::GroupRepresentation)
    F = base_ring(X)
    n,m = int_dim(X),int_dim(Y)
    map = zero(matrix_space(F,n*m,n*m))
    for i ∈ 1:n, j ∈ 1:m
        v1 = matrix(F,transpose([k == i ? 1 : 0 for k ∈ 1:n]))
        v2 = matrix(F,transpose([k == j ? 1 : 0 for k ∈ 1:m]))
        map[(j-1)*n + i, :] = kronecker_product(v1,v2)
    end
    return morphism(X⊗Y, Y⊗X, transpose(map),check = false)
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
        return τ,[GroupRepresentationMorphism(ρ,τ,zero(matrix_space(F,0,int_dim(τ)))), id(τ)], [GroupRepresentationMorphism(τ,ρ,zero(matrix_space(F,int_dim(τ),0))), id(τ)]
    elseif τ.m == 0
        return ρ,[id(ρ), GroupRepresentationMorphism(τ,ρ,zero(matrix_space(F,0,int_dim(ρ)))), id(τ)], [id(ρ), GroupRepresentationMorphism(ρ,τ,zero(matrix_space(F,int_dim(ρ),0)))]
    end

    M1 = matrix_space(F,int_dim(ρ),int_dim(ρ))
    M2 = matrix_space(F,int_dim(ρ),int_dim(τ))
    M3 = matrix_space(F,int_dim(τ),int_dim(ρ))
    M4 = matrix_space(F,int_dim(τ),int_dim(τ))

    generators = order(grp) == 1 ? elements(grp) : gens(grp)

    S = Representation(parent(ρ), generators, [[matrix(ρ(g)) zero(M2); zero(M3) matrix(τ(g))] for g ∈ generators], check = false)


    incl_ρ = morphism(ρ, S, [one(M1) zero(M2)], check = false)
    incl_τ = morphism(τ, S, [zero(M3) one(M4)], check = false)
    proj_ρ = morphism(S, ρ, [one(M1); zero(M3)], check = false)
    proj_τ = morphism(S, τ, [zero(M2); one(M4)], check = false)

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

    z1 = zero(matrix_space(F, int_dim(domain(f)), int_dim(codomain(g))))
    z2 = zero(matrix_space(F, int_dim(domain(g)), int_dim(codomain(f))))

    m = [matrix(f) z1; z2 matrix(g)]

    return morphism(dom,codom, m, check = false)
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
function simples(Rep::GroupRepresentationCategory)
    
    return get_attribute!(Rep, :simples) do
        #@show Rep.group
        grp = base_group(Rep)
        F = base_ring(Rep)

        if order(grp) == 1 return [one(Rep)] end

        #gap_field = GAP.Globals.FiniteField(Int(characteristic(F)), degree(F))
        gap_field = codomain(iso_oscar_gap(F))
        gap_reps = if is_finite(F) 
            GAP.Globals.IrreducibleRepresentations(grp.X,gap_field)
        else
            GAP.Globals.IrreducibleRepresentations(grp.X)
        end

        int_dims = [GAP.Globals.DimensionOfMatrixGroup(GAP.Globals.Range(m)) for m ∈ gap_reps]
    
        oscar_reps = [GAPGroupHomomorphism(grp, GL(int_dims[i],F), gap_reps[i]) for i ∈ 1:length(gap_reps)]
        reps = [GroupRepresentation(Rep,grp,m,F,d) for (m,d) ∈ zip(oscar_reps,int_dims)]

        return reps

    end
end


"""
    decompose(σ::GroupRepresentation)

Decompose the representation into a direct sum of simple objects. Return a
list of tuples with simple objects and multiplicities.
"""
#= @memoize Dict =# function decompose(σ::GroupRepresentation)
    F = base_ring(σ)
    if int_dim(σ) == 0 return [] end
    G = σ.group

    if order(G) == 1 return [(one(parent(σ)),int_dim(σ))] end

    M = to_gap_module(σ,F)
    ret = Object[]
    facs = GAP.Globals.MTX.Indecomposition(M)
    d = int_dim(σ)
    for m ∈ facs
        imgs = [matrix(F,[F(n[i,j]) for i ∈ 1:length(n), j ∈ 1:length(n)]) for n ∈ m[2].generators]
        ret = [ret; Representation(parent(σ),gens(G),imgs, check = false)]
    end
    uniques = unique_indecomposables(ret)
    [(s, length(findall(r -> is_isomorphic(s,r)[1], ret))) for s ∈ uniques]
end

# function indecomposable_subobjects(ρ::GroupRepresentation)
#     [x for (x,k) ∈ decompose(ρ)]
# end

function simple_subobjects(σ::GroupRepresentation)
    F = base_ring(σ)
    if int_dim(σ) == 0 return [] end
    G = σ.group

    if order(G) == 1 return [(one(parent(σ)),int_dim(σ))] end

    M = to_gap_module(σ,F)
    ret = []
    facs = GAP.Globals.MTX.CollectedFactors(M)
    d = int_dim(σ)
    for m ∈ facs
        imgs = [matrix(F,[F(n[i,j]) for i ∈ 1:length(n), j ∈ 1:length(n)]) for n ∈ m[1].generators]
        ret = [ret; Representation(parent(σ),gens(G),imgs, check = false)]
    end
    ret
end

function regular_representation(C::GroupRepresentationCategory)
    G = base_group(C)
    H,_ = trivial_subgroup(G)
    RepH = representation_category(base_ring(C),H)
    I = induction(one(RepH), G)
    GroupRepresentation(C, I.group, I.m, I.base_ring, I.int_dim)
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
#= @memoize Dict =# function Hom(σ::GroupRepresentation, τ::GroupRepresentation)
    grp = base_group(σ)
    F = base_ring(σ)

    if int_dim(σ)*int_dim(τ) == 0 return GRHomSpace(σ,τ,GroupRepresentationMorphism[],VectorSpaces(F)) end

    gap_to_F = iso_oscar_gap(F)
    gap_F = codomain(gap_to_F)
    generators = order(grp) == 1 ? elements(grp) : gens(grp)

    #Build the modules from σ and τ
    mats_σ = GAP.GapObj([GAP.GapObj(σ(g)) for g ∈ generators])
    mats_τ = GAP.GapObj([GAP.GapObj(τ(g)) for g ∈ generators])
    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
    Mτ = GAP.Globals.GModuleByMats(mats_τ, gap_F)

    # Use GAPs Meat Axe to calculate a basis
    gap_homs = GAP.Globals.MTX.BasisModuleHomomorphisms(Mσ,Mτ)

    int_dims_m,int_dims_n = int_dim(σ), int_dim(τ)

    mat_homs = [matrix(F,[preimage(gap_to_F, m[i,j]) for i ∈ 1:int_dims_m, j ∈ 1:int_dims_n]) for m ∈ gap_homs]

    rep_homs = [morphism(σ,τ,m,check = false) for m ∈ mat_homs]

    return GRHomSpace(σ,τ, rep_homs, VectorSpaces(F))
end

function zero(H::GRHomSpace)
    dom = H.X
    codom = H.Y
    m = zero(matrix_space(base_ring(dom),int_dim(dom),int_dim(codom)))
    return morphism(dom,codom,m, check = false)
end

function zero_morphism(X::GroupRepresentation, Y::GroupRepresentation)
    m = zero(matrix_space(base_ring(X),int_dim(X),int_dim(Y)))
    return morphism(X,Y,m, check = false)
end

#-------------------------------------------------------------------------
#   Restriction and Induction Functor
#-------------------------------------------------------------------------

function restriction(ρ::GroupRepresentation, C::GroupRepresentationCategory)
    H = base_group(C)
    res = restriction(ρ, H)
    GroupRepresentation(C, res.group, res.m, res.base_ring, res.int_dim)
end

function restriction(ρ::GroupRepresentation, H::GAPGroup)
    b,f = is_subgroup(H, ρ.group)
    RepH = representation_category(base_ring(ρ),H)
    if b == false throw(ErrorException("Not a sub")) end
    if ρ.m == 0 return zero(RepH) end
    h = hom(H,codomain(ρ.m), gens(H), [ρ(f(g)) for g ∈ gens(H)], check = false)
    return GroupRepresentation(RepH, H, h, base_ring(ρ), int_dim(ρ))
end

function restriction(f::GroupRepresentationMorphism, C::GroupRepresentationCategory)
    H = base_group(C)
    if domain(f).group == H return f end
    return morphism(restriction(domain(f),C), restriction(codomain(f),C), matrix(f), check = false)
end

function restriction(f::GroupRepresentationMorphism, H::GAPGroup)
    restriction(f, GroupRepresentationCategory(H, base_ring(f)))
end

function induction(ρ::GroupRepresentation, C::GroupRepresentationCategory)
    G = base_group(C)
    ind = induction(ρ, G)
    GroupRepresentation(C,ind.group, ind.m, ind.base_ring, ind.int_dim)
end

function induction(ρ::GroupRepresentation, G::GAPGroup)
    H = ρ.group

    if H == G return ρ end

    if !is_subgroup(H,G)[1] throw(ErrorException("Not a supergroup")) end

    if ρ.m == 0 return zero(representation_category(base_ring(ρ),G)) end

    transversal = left_transversal(G,H)

    g = order(G) == 1 ? elements(G) : gens(G)

    ji = [[findfirst(x -> g[k]*t ∈ orbit(gset(H, (y,g) -> y*g, G), x), transversal) for t ∈ transversal] for k ∈ 1:length(g)]
    g_ji = [[transversal[i] for i ∈ m] for m ∈ ji]

    hi = [[inv(g_ji[k][i])*g[k]*transversal[i] for i ∈ 1:length(transversal)] for k ∈ 1:length(g)]

    images = []
    d = int_dim(ρ)
    n = length(transversal)*d
    for i ∈ 1:length(g)
        m = zero(matrix_space(base_ring(ρ), n, n))

        for j ∈ 1:length(transversal)
            m[ (ji[i][j]-1)*d+1:ji[i][j]*d, (j-1)*d+1:j*d] = matrix(ρ(hi[i][j]))
        end
        images = [images; m]
    end
    return Representation(G, g, images, check = false)
end

function induction(f::GroupRepresentationMorphism, C::GroupRepresentationCategory)
    G = base_group(C)
    dom = induction(domain(f), C)
    codom = induction(codomain(f), C)
    return morphism(dom,codom, direct_sum([morphism(matrix(f)) for i ∈ 1:Int64(index(G,domain(f).group))]).m, check = false)
end

function induction(f::GroupRepresentationMorphism, G::GAPGroup)
    induction(f, GroupRepresentationCategory(G,base_ring(f)))
end


#-------------------------------------------------------------------------------
#   Restriction and Induction Functors
#-------------------------------------------------------------------------------

struct GRepRestriction <: AbstractFunctor
    domain::GroupRepresentationCategory
    codomain::GroupRepresentationCategory
    obj_map
    mor_map
end

function Restriction(C::GroupRepresentationCategory, D::GroupRepresentationCategory)
    @assert base_ring(C) == base_ring(D) "Not compatible"
    #@assert issubgroup(base_group(D), base_group(C))[1] "Not compatible"
    obj_map = X -> restriction(X, base_group(D))
    mor_map = f -> restriction(f, base_group(D))
    return GRepRestriction(C,D,obj_map,mor_map)
end

struct GRepInduction <: AbstractFunctor
    domain::GroupRepresentationCategory
    codomain::GroupRepresentationCategory
    obj_map
    mor_map
end

function Induction(C::GroupRepresentationCategory, D::GroupRepresentationCategory)
    @assert base_ring(C) == base_ring(D) "Not compatible"
    #@assert issubgroup(base_group(C), base_group(D))[1] "Not compatible"
    obj_map = X -> induction(X, base_group(D))
    mor_map = f -> induction(f, base_group(D))
    return GRepInduction(C,D, obj_map, mor_map)
end

function show(io::IO, F::GRepRestriction)
    print(io,"Restriction functor from $(domain(F)) to $(codomain(F)).")
end

function show(io::IO, F::GRepInduction)
    print(io,"Induction functor from $(domain(F)) to $(codomain(F)).")
end

#-------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------


function show(io::IO, Rep::GroupRepresentationCategory)
    print(io, """Representation Category of $(Rep.group) over $(Rep.base_ring)""")
end

function show(io::IO, ρ::GroupRepresentation)
    print(io,"$(int_dim(ρ))-dimensional group representation over $(base_ring(ρ)) of $(ρ.group))")
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
    mats_σ = GAP.GapObj([GAP.GapObj(σ(g)) for g ∈ gens(grp)])
    Mσ = GAP.Globals.GModuleByMats(mats_σ, gap_F)
end

function express_in_basis(f::GroupRepresentationMorphism, basis::Vector{GroupRepresentationMorphism})
    o = one(base_group(domain(f)))
    express_in_basis(morphism(f.map), [morphism(g.map) for g in basis])
end
