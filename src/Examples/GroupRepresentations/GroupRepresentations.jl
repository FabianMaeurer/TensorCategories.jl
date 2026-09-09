@attributes mutable struct GroupRepresentationCategory <: RepresentationCategory
    group::Group
    base_ring::Field

    function GroupRepresentationCategory(G::Group, F::Field) 
        C = new(G,F)
        set_attribute!(C, :semisimple,
            characteristic(F) == 0 || rem(order(G), characteristic(F)) != 0)
        set_attribute!(C, :tensor, true)
        C
    end
end

struct GroupRepresentation <: RepresentationObject
    parent::GroupRepresentationCategory
    group::Group
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
# For finite G, Rep_k(G) is the module category of the finite-dimensional
# algebra kG, independently of Maschke's semisimplicity condition.
is_finite(C::GroupRepresentationCategory) = is_finite(base_group(C))
# The ordinary flip is equivariant for the diagonal action in every
# characteristic, so Rep(G) is symmetric monoidal.
is_braided(::GroupRepresentationCategory) = true
is_weak_fusion(C::GroupRepresentationCategory) = is_semisimple(C)
is_weak_multifusion(C::GroupRepresentationCategory) = is_weak_fusion(C)
function is_fusion(C::GroupRepresentationCategory)
    get_attribute!(C, :fusion) do
        is_weak_fusion(C) && is_split_semisimple(C)
    end
end

"""
    is_split_semisimple(C::GroupRepresentationCategory)

Return whether `C` is semisimple and every simple representation is
absolutely simple over its coefficient field.
"""
function is_split_semisimple(C::GroupRepresentationCategory)
    get_attribute!(C, :split_semisimple) do
        is_semisimple(C) || return false
        F = base_ring(C)
        G = base_group(C)
        if F isa QQField
            # Over QQ, split semisimplicity is equivalent to every ordinary
            # irreducible character being rational-valued with Schur index one.
            return all(χ -> Oscar.degree_of_character_field(χ) == 1 &&
                            Oscar.schur_index(χ) == 1,
                       collect(Oscar.character_table(G)))
        elseif F isa Oscar.QQAbField || F isa QQBarField
            return true
        end
        all(X -> int_dim(End(X)) == 1, simples(C))
    end
end

"""
    is_finite_representation_type(C::GroupRepresentationCategory)

Return whether `C` has only finitely many isomorphism classes of
indecomposable objects. In modular characteristic this uses Higman's
criterion: a Sylow subgroup for the characteristic is cyclic.
"""
function is_finite_representation_type(C::GroupRepresentationCategory)
    is_semisimple(C) && return true
    p = Int(characteristic(base_ring(C)))
    Oscar.is_cyclic(Oscar.sylow_subgroup(base_group(C), p)[1])
end

"""
    splitting_field(C::GroupRepresentationCategory)

Return a splitting field for the finite-group representations in `C`.
The returned field contains a primitive root of unity whose order is the
prime-to-characteristic part of the exponent of the group. By Brauer's
splitting theorem this is sufficient, but it need not be a minimal splitting
field.
"""
function splitting_field(C::GroupRepresentationCategory)
    F = base_ring(C)
    Oscar.is_exact_type(F) || throw(ArgumentError(
        "splitting_field requires an exact coefficient field"))
    F isa Oscar.QQAbField && return F
    F isa QQBarField && return F
    m = Int(exponent(base_group(C)))
    p = Int(characteristic(F))
    if p != 0
        while m % p == 0
            m = div(m, p)
        end
    end
    m == 1 && return F
    Fx, x = polynomial_ring(F, "x"; cached=false)
    Oscar.splitting_field(x^m - one(Fx))
end

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
    representation_category(F::Field, G::Group)

Category of finite-dimensional group representations of ``G`` over ``F``.
"""
function representation_category(F::Field, G::Group)
    return GroupRepresentationCategory(G,F)
end

function representation_category(G::Group)
    return representation_category(abelian_closure(QQ)[1], G)
end

"""
    rep(F::Field, G::Group)
    rep(G::Group)

Short aliases for `representation_category(F,G)` and
`representation_category(G)`, respectively. The one-argument form uses the
abelian closure of the rational field.
"""
rep(F::Field, G::Group) = representation_category(F,G)
rep(G::Group) = representation_category(G)

function extension_of_scalars(C::GroupRepresentationCategory,L::Field;
                              embedding=_scalar_extension_embedding(base_ring(C),L))
    representation_category(L,base_group(C))
end

function extension_of_scalars(X::GroupRepresentation,L::Field,
                              D::GroupRepresentationCategory=
                                  extension_of_scalars(parent(X),L);
                              embedding=_scalar_extension_embedding(base_ring(X),L))
    base_ring(D) == L && base_group(D) == base_group(X) ||
        throw(ArgumentError("incompatible target representation category"))
    int_dim(X) == 0 && return zero(D)
    G = base_group(X)
    generators = order(G) == 1 ? elements(G) : gens(G)
    matrices = [matrix(L,int_dim(X),int_dim(X),
                       embedding.(collect(matrix(X(g))))) for g in generators]
    Representation(D,generators,matrices)
end

function extension_of_scalars(f::GroupRepresentationMorphism,L::Field,
                              D::GroupRepresentationCategory=
                                  extension_of_scalars(parent(f),L);
                              embedding=_scalar_extension_embedding(base_ring(f),L))
    morphism(extension_of_scalars(domain(f),L,D;embedding),
             extension_of_scalars(codomain(f),L,D;embedding),
             matrix(L,size(matrix(f))...,embedding.(collect(matrix(f)))))
end

"""
    Representation(G::Group, pre_img::Vector, img::Vector)

Group representation defined by the images of generators of G. The images
must define a representation; pass `check=true` to verify the group relations.
"""
function Representation(G::Group, pre_img::Vector, img::Vector; check::Bool = false)
    F = base_ring(img[1])
    Representation(representation_category(F, G),pre_img, img, check = check)
end

function Representation(C::GroupRepresentationCategory, pre_img::Vector, img::Vector; check::Bool = false)
    isempty(img) && throw(ArgumentError(
        "supply generator images, including an identity image for the trivial group"))
    length(pre_img) == length(img) || throw(ArgumentError(
        "mismatching generator and image lists"))
    F = base_ring(C)
    all(m -> base_ring(m) == F, img) || throw(ArgumentError(
        "generator matrices have the wrong base field"))
    d = size(img[1], 1)
    all(m -> size(m) == (d, d), img) || throw(ArgumentError(
        "generator matrices must be square and have equal sizes"))
    # A finite-group representation over an infinite field lands in the
    # finitely generated matrix image. Constructing GL(d,F) instead fails,
    # for example, because GL(1,QQ) has no known finite generating set.
    H = is_finite(F) ? GL(d, F) : Oscar.matrix_group(F, d, img; check = check)
    G = base_group(C)
    m = hom(G,H, pre_img, H.(img), check = check)
    GroupRepresentation(C,G,m,F,d)
end

"""
    Representation(G::Group, m::Function)

Group representation defined by m:G -> Mat_n.
"""
function Representation(G::Group, m::Function)
    F = order(G) == 1 ? base_ring(parent(m(elements(G)[1]))) : base_ring(parent(m(G[1])))
    Representation(representation_category(F,G),m)
end

function Representation(C::GroupRepresentationCategory, m::Function)
    G = base_group(C)
    generators = order(G) == 1 ? elements(G) : gens(G)
    images = [m(g) for g in generators]
    Representation(C, generators,
        [a isa MatElem ? a : matrix(a) for a in images])
end
"""
    morphism(ρ::GroupRepresentation, τ::GroupRepresentation, m::MatElem; check = false)

Morphism between representations defined by ``m``. Equivariance is assumed by
default; pass `check=true` to verify it.
"""
function morphism(ρ::GroupRepresentation, τ::GroupRepresentation, m::MatElem; check = false)
    parent(ρ) == parent(τ) || throw(ArgumentError("Mismatching parents"))
    base_ring(m) == base_ring(ρ) || throw(ArgumentError("Mismatching base fields"))
    size(m) == (int_dim(ρ), int_dim(τ)) || throw(ErrorException("Mismatching dimensions"))
    check && !isequivariant(m, ρ, τ) && throw(ArgumentError("Map must be equivariant"))
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

tr(ρ::GroupRepresentationMorphism) = morphism(one(parent(ρ)), one(parent(ρ)), matrix(base_ring(ρ),1,1,[tr(matrix(ρ))]))

is_zero(f::GroupRepresentationMorphism) = iszero(matrix(f))

is_monomorphism(f::GroupRepresentationMorphism) = rank(matrix(f)) == int_dim(domain(f))
is_epimorphism(f::GroupRepresentationMorphism) = rank(matrix(f)) == int_dim(codomain(f))

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
    parent(ρ) == parent(τ) || return false
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
#=  =# function is_isomorphic(σ::GroupRepresentation, τ::GroupRepresentation)
    parent(σ) == parent(τ) || return false, nothing

    if int_dim(σ) != int_dim(τ) return false, nothing end
    if int_dim(σ) == 0 return true, zero_morphism(σ,τ) end

    F = base_ring(σ)
    grp = σ.group

    if order(grp) == 1 return true, morphism(σ,τ,one(matrix_space(F,int_dim(σ),int_dim(τ)))) end

    if !is_finite(F)
        characteristic(F) == 0 || throw(ArgumentError(
            "isomorphism search over infinite positive-characteristic fields is unsupported"))
        B = basis(Hom(σ,τ))
        isempty(B) && return false,nothing
        # A basis map is often already invertible, avoiding a symbolic
        # determinant in the common case.
        for f in B
            !iszero(det(matrix(f))) && return true,f
        end
        R,x = polynomial_ring(F,length(B))
        M = sum((x[i]*matrix(B[i]) for i in eachindex(B));
                init=zero_matrix(R,int_dim(σ),int_dim(σ)))
        determinant = det(M)
        iszero(determinant) && return false,nothing
        # A nonzero polynomial of degree at most n in each variable cannot
        # vanish on {0,...,n}^r in characteristic zero. Specialize greedily
        # while retaining a nonzero polynomial.
        values = elem_type(F)[]
        for i in eachindex(B)
            found = false
            for a in 0:int_dim(σ)
                candidate = [R.(values);R(a);x[i+1:end]]
                if !iszero(Oscar.evaluate(determinant,candidate))
                    push!(values,F(a))
                    found = true
                    break
                end
            end
            found || error("nonzero determinant specialization failed")
        end
        f = sum((values[i]*B[i] for i in eachindex(B));
                init=zero_morphism(σ,τ))
        return true,f
    end

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
    #@assert codomain(f) == domain(g) "Morphisms not compatible"
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
        return ρ,[id(ρ), GroupRepresentationMorphism(τ,ρ,zero(matrix_space(F,0,int_dim(ρ))))], [id(ρ), GroupRepresentationMorphism(ρ,τ,zero(matrix_space(F,int_dim(ρ),0)))]
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
    simples(Rep::GroupRepresentationCategory; backend=:auto)

Return representatives of simple objects. `backend=:gap` uses GAP's ordinary
or modular representation routines. `backend=:hecke` obtains the simple
factors of the regular module with Hecke's matrix-module algorithms. The
default `:auto` uses GAP and reports when a nonsplit characteristic-zero field
requires a rational/Schur-index backend.
"""
function simples(Rep::GroupRepresentationCategory; backend::Symbol=:auto)
    backend in (:auto, :gap, :hecke) || throw(ArgumentError(
        "backend must be :auto, :gap, or :hecke"))
    key = backend === :auto ? :simples : Symbol(:simples_, backend)
    get_attribute!(Rep, key) do
        grp, F = base_group(Rep), base_ring(Rep)
        order(grp) == 1 && return [one(Rep)]
        backend === :hecke && return first.(_composition_factors_hecke(
            regular_representation(Rep)))

        gap_field = is_finite(F) ? codomain(iso_oscar_gap(F)) : nothing
        gap_reps = is_finite(F) ?
            GAP.Globals.IrreducibleRepresentations(grp.X, gap_field) :
            GAP.Globals.IrreducibleRepresentations(grp.X)
        try
            generators = gens(grp)
            [Representation(Rep, generators,
                [matrix(F, GAP.Globals.Image(m, g.X)) for g in generators];
                check=false) for m in gap_reps]
        catch
            backend === :gap && rethrow()
            throw(ArgumentError(
                "GAP's absolutely irreducible representations are not realizable over this coefficient field; try backend=:hecke for a field-aware matrix-module computation"))
        end
    end
end

function _hecke_module(σ::GroupRepresentation)
    G = base_group(σ)
    generators = order(G) == 1 ? elements(G) : gens(G)
    if !is_finite(base_ring(σ)) && length(generators) > 1
        throw(ArgumentError(
            "Hecke's matrix-module algorithms are not currently reliable for multiple generators over an infinite field"))
    end
    Oscar.Hecke.Amodule([matrix(σ(g)) for g in generators])
end

function _representation_from_hecke(C::GroupRepresentationCategory, M)
    Representation(C, gens(base_group(C)), Oscar.Hecke.action_of_gens(M);
                   check=false)
end

function _composition_factors_hecke(σ::GroupRepresentation)
    C = parent(σ)
    [(_representation_from_hecke(C,M), Int(n)) for (M,n) in
        Oscar.Hecke.composition_factors_with_multiplicity(_hecke_module(σ))]
end

"""
    indecomposables(C::GroupRepresentationCategory; backend=:auto)

Enumerate all indecomposable representations when supported. In a semisimple
category these are the simple representations. The finite-representation-type
predicate is available in modular characteristic, but the installed backends
do not provide a general enumeration there.
"""
function indecomposables(C::GroupRepresentationCategory; backend::Symbol=:auto)
    is_semisimple(C) && return simples(C; backend)
    is_finite_representation_type(C) || throw(ArgumentError(
        "the representation category has infinite representation type"))
    throw(ArgumentError(
        "enumeration of all indecomposable modular representations is not implemented"))
end


"""
    decompose(σ::GroupRepresentation; backend=:auto)

Decompose the representation into a direct sum of indecomposable objects.
Return a list of tuples with indecomposable objects and multiplicities.
"""
#=  =# function decompose(σ::GroupRepresentation; backend::Symbol=:auto)
    backend in (:auto, :gap, :hecke) || throw(ArgumentError(
        "backend must be :auto, :gap, or :hecke"))
    F = base_ring(σ)
    if int_dim(σ) == 0 return [] end
    G = σ.group

    if order(G) == 1 return [(one(parent(σ)),int_dim(σ))] end

    if backend === :hecke
        is_semisimple(parent(σ)) || throw(ArgumentError(
            "Hecke computes composition factors, not indecomposable summands in a nonsemisimple category"))
        return _composition_factors_hecke(σ)
    elseif !is_finite(F)
        backend === :gap && throw(ArgumentError(
            "GAP's MeatAxe backend requires a finite coefficient field"))
        is_semisimple(parent(σ)) || throw(ArgumentError(
            "Krull-Schmidt decomposition over this field is unsupported"))
        return _semisimple_representation_factors(σ)
    end

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

"""
    composition_factors(σ::GroupRepresentation; backend=:auto)

Return simple composition factors with their Jordan--Hölder multiplicities.
These are subquotients, not necessarily subobjects or direct summands.
"""
function composition_factors(σ::GroupRepresentation; backend::Symbol=:auto)
    backend in (:auto, :gap, :hecke) || throw(ArgumentError(
        "backend must be :auto, :gap, or :hecke"))
    F = base_ring(σ)
    if int_dim(σ) == 0 return Tuple{GroupRepresentation,Int}[] end
    G = σ.group

    if order(G) == 1 return [(one(parent(σ)),int_dim(σ))] end

    backend === :hecke && return _composition_factors_hecke(σ)
    if !is_finite(F)
        backend === :gap && throw(ArgumentError(
            "GAP's MeatAxe backend requires a finite coefficient field"))
        is_semisimple(parent(σ)) || throw(ArgumentError(
            "composition factors over this field are unsupported"))
        return _semisimple_representation_factors(σ)
    end

    M = to_gap_module(σ,F)
    ret = Tuple{GroupRepresentation,Int}[]
    # GAP Reference Manual 69.7-11: CollectedFactors records frequencies.
    facs = GAP.Globals.MTX.CollectedFactors(M)
    for m ∈ facs
        imgs = [matrix(F,[F(n[i,j]) for i ∈ 1:length(n), j ∈ 1:length(n)]) for n ∈ m[1].generators]
        push!(ret,(Representation(parent(σ),gens(G),imgs,check=false),Int(m[2])))
    end
    ret
end

function _semisimple_representation_factors(σ::GroupRepresentation)
    [(S, div(int_dim(Hom(S,σ)), int_dim(End(S))))
     for S in simples(parent(σ)) if int_dim(Hom(S,σ)) != 0]
end

"""
    simple_subobjects(σ::GroupRepresentation)

Return the simple isomorphism types in the socle, without multiplicities.
"""
function simple_subobjects(σ::GroupRepresentation; backend::Symbol=:auto)
    [S for (S,_) in composition_factors(σ; backend) if int_dim(Hom(S,σ)) != 0]
end

function is_simple(σ::GroupRepresentation; backend::Symbol=:auto)
    backend in (:auto, :gap, :hecke) || throw(ArgumentError(
        "backend must be :auto, :gap, or :hecke"))
    int_dim(σ) == 0 && return false
    int_dim(σ) == 1 && return true
    order(base_group(σ)) == 1 && return false
    if backend === :hecke
        return Oscar.Hecke.meataxe(_hecke_module(σ))[1]
    end
    backend === :gap && !is_finite(base_ring(σ)) && throw(ArgumentError(
        "GAP's MeatAxe backend requires a finite coefficient field"))
    if !(base_ring(σ) isa FinField)
        characteristic(base_ring(σ)) == 0 || throw(ArgumentError(
            "irreducibility over this field is unsupported"))
        # Maschke and Schur: in characteristic zero a division endomorphism
        # algebra certifies simplicity. Over QQ the commutative case can be
        # decided as a field; noncommutative division algebras need a separate
        # backend and must not be confused with matrix algebras.
        E = End(σ)
        int_dim(E) == 1 && return true
        A = endomorphism_ring_by_basis(σ,basis(E))
        is_commutative(A) && return is_simple(A)
        throw(ArgumentError(
            "irreducibility with noncommutative endomorphism algebra requires a division-algebra backend"))
    end
    # GAP Reference Manual 69.5-1: test irreducibility itself; one factor
    # type does not imply that a module has composition length one.
    Bool(GAP.Globals.MTX.IsIrreducible(to_gap_module(σ,base_ring(σ))))
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
#=  =# function Hom(σ::GroupRepresentation, τ::GroupRepresentation)
    parent(σ) == parent(τ) || throw(ArgumentError(
        "representations belong to different categories"))
    grp = base_group(σ)
    F = base_ring(σ)

    if int_dim(σ)*int_dim(τ) == 0 return GRHomSpace(σ,τ,GroupRepresentationMorphism[],VectorSpaces(F)) end

    if !is_finite(F)
        # Row-vector convention: rho(g)M=M tau(g). These homogeneous exact
        # linear equations compute Hom without a finite-field MeatAxe.
        n,m = int_dim(σ),int_dim(τ)
        generators = order(grp) == 1 ? elements(grp) : gens(grp)
        E = [matrix(F,n,m,[Int(i == a && j == b)
                           for i in 1:n,j in 1:m])
             for a in 1:n for b in 1:m]
        A = zero_matrix(F,length(generators)*n*m,n*m)
        for (k,g) in enumerate(generators),(l,e) in enumerate(E)
            c = collect(matrix(σ(g))*e-e*matrix(τ(g)))[:]
            for i in eachindex(c)
                A[(k-1)*n*m+i,l] = c[i]
            end
        end
        d,N = nullspace(A)
        maps = GroupRepresentationMorphism[morphism(σ,τ,
            sum((N[i,j]*E[i] for i in eachindex(E));
                init=zero_matrix(F,n,m))) for j in 1:d]
        return GRHomSpace(σ,τ,maps,VectorSpaces(F))
    end

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

function restriction(ρ::GroupRepresentation, H::Group)
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

function restriction(f::GroupRepresentationMorphism, H::Group)
    restriction(f, GroupRepresentationCategory(H, base_ring(f)))
end

function induction(ρ::GroupRepresentation, C::GroupRepresentationCategory)
    G = base_group(C)
    ind = induction(ρ, G)
    GroupRepresentation(C,ind.group, ind.m, ind.base_ring, ind.int_dim)
end

function induction(ρ::GroupRepresentation, G::Group)
    H = ρ.group

    if H == G return ρ end

    if !is_subgroup(H,G)[1] throw(ErrorException("Not a supergroup")) end

    if ρ.m == 0 return zero(representation_category(base_ring(ρ),G)) end

    transversal = left_transversal(G,H)

    g = order(G) == 1 ? elements(G) : gens(G)

    ji = [[findfirst(x -> g[k]*t ∈ orbit(gset(H, (y,g) -> y*g, G), x), transversal) for t ∈ transversal] for k ∈ 1:length(g)]
    g_ji = [[transversal[i] for i ∈ m] for m ∈ ji]

    hi = [[inv(g_ji[k][i])*g[k]*transversal[i] for i ∈ 1:length(transversal)] for k ∈ 1:length(g)]

    images = MatElem[]
    d = int_dim(ρ)
    n = length(transversal)*d
    for i ∈ 1:length(g)
        m = zero(matrix_space(base_ring(ρ), n, n))

        for j ∈ 1:length(transversal)
            m[ (ji[i][j]-1)*d+1:ji[i][j]*d, (j-1)*d+1:j*d] = matrix(ρ(hi[i][j]))
        end
        push!(images, m)
    end
    return Representation(G, g, images, check = false)
end

function induction(f::GroupRepresentationMorphism, C::GroupRepresentationCategory)
    G = base_group(C)
    dom = induction(domain(f), C)
    codom = induction(codomain(f), C)
    return morphism(dom,codom, direct_sum([morphism(matrix(f)) for i ∈ 1:Int64(index(G,domain(f).group))]).m, check = false)
end

function induction(f::GroupRepresentationMorphism, G::Group)
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
