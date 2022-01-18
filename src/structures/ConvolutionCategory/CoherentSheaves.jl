struct CohSheaves{T,G} <: MultiTensorCategory{T}
    group::G
    base_ring::Field
    GSet::GSet
    orbit_reps
    orbit_stabilizers
end

struct CohSheaf{T,G} <: Object
    parent::CohSheaves{T}
    stalks::Vector{GroupRepresentation{T,G}}
end

struct CohSheafMorphism{T,G} <: Morphism
    domain::CohSheaf{T,G}
    codomain::CohSheaf{T,G}
    m::Vector{GroupRepresentationMorphism{T,G}}
end

#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------
"""
    CohSheaves(X::GSet,F::Field)

The category of ``G``-equivariant coherent sheafes on ``X``.
"""
function CohSheaves(X::GSet, F::Field)
    G = X.group
    orbit_reps = [O.seeds[1] for O ∈ orbits(X)]
    orbit_stabilizers = [stabilizer(G,x,X.action_function)[1] for x ∈ orbit_reps]
    return CohSheaves{elem_type(F), typeof(G)}(G, F, X, orbit_reps, orbit_stabilizers)
end

"""
    CohSheaves(X, F::Field)

The category of coherent sheafes on ``X``.
"""
function CohSheaves(X,F::Field)
    G = symmetric_group(1)
    return CohSheaves(gset(G,X), F)
end

#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------
"""
    issemisimple(C::CohSheaves)

Return whether ``C``is semisimple.
"""
issemisimple(C::CohSheaves) = gcd(order(C.group), characteristic(base_ring(C))) == 1

"""
    stalks(X::CohSheaf)

Return the stalks of ``X``.
"""
stalks(X::CohSheaf) = X.stalks

function orbit_index(X::CohSheaf, y)
    i = findfirst(x -> y ∈ x, orbits(parent(X).GSet))
end

function orbit_index(C::CohSheaves, y)
    i = findfirst(x -> y ∈ x, orbits(C.GSet))
end

function stalk(X::CohSheaf,y)
    return stalks(X)[orbit_index(X,y)]
end

"""
    zero(C::CohSheaves)

Return the zero sheaf on the ``G``-set.
"""
zero(C::CohSheaves) = CohSheaf(C,[zero(RepresentationCategory(H,base_ring(C))) for H ∈ C.orbit_stabilizers])

zero_morphism(X::CohSheaf, Y::CohSheaf) = CohSheafMorphism(X,Y,[zero(Hom(x,y)) for (x,y) ∈ zip(stalks(X),stalks(Y))])

function ==(X::CohSheaf, Y::CohSheaf)
    if parent(X) != parent(Y) return false end
    for (s,r) ∈ zip(stalks(X),stalks(Y))
        if s != r return false end
    end
    return true
end

"""
    isisomorphic(X::CohSheaf{T,G}, Y::CohSheaf{T,G}) where {T,G}

Check whether ``X``and ``Y`` are isomorphic and the isomorphism if possible.
"""
function isisomorphic(X::CohSheaf{T,G}, Y::CohSheaf{T,G}) where {T,G}
    m = GroupRepresentationMorphism{T,G}[]
    for (s,r) ∈ zip(stalks(X),stalks(Y))
        b, iso = isisomorphic(s,r)
        if !b return false, nothing end
        m = [m; iso]
    end
    return true, CohSheafMorphism{T,G}(X,Y,m)
end

==(f::CohSheafMorphism, g::CohSheafMorphism) = f.m == f.m

id(X::CohSheaf{T,G}) where {T,G} = CohSheafMorphism{T,G}(X,X,[id(s) for s ∈ stalks(X)])

associator(X::CohSheaf, Y::CohSheaf, Z::CohSheaf) = id(X⊗Y⊗Z)

#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    dsum(X::CohSheaf{T,G}, Y::CohSheaf{T,G}, morphisms::Bool = false) where {T,G}

Return the direct sum of sheaves. Return also the inclusion and projection if
morphisms = true.
"""
function dsum(X::CohSheaf{T,G}, Y::CohSheaf{T,G}, morphisms::Bool = false) where {T,G}
    sums = [dsum(x,y,true) for (x,y) ∈ zip(stalks(X), stalks(Y))]
    Z = CohSheaf{T,G}(parent(X), [s[1] for s ∈ sums])

    if !morphisms return Z end

    ix = [CohSheafMorphism{T,G}(x,Z,[s[2][i] for s ∈ sums]) for (x,i) ∈ zip([X,Y],1:2)]
    px = [CohSheafMorphism{T,G}(Z,x,[s[3][i] for s ∈ sums]) for (x,i) ∈ zip([X,Y],1:2)]
    return Z,ix,px
end

"""
    dsum(f::CohSheafMorphism{T,G}, g::CohSheafMorphism{T,G}) where {T,G}

Return the direct sum of morphisms of sheaves.
"""
function dsum(f::CohSheafMorphism{T,G}, g::CohSheafMorphism{T,G}) where {T,G}
    dom = dsum(domain(f), domain(g))[1]
    codom = dsum(codomain(f), codomain(g))[1]
    mors = [dsum(m,n) for (m,n) ∈ zip(f.m,g.m)]
    return CohSheafMorphism{T,G}(dom,codom, mors)
end

product(X::CohSheaf,Y::CohSheaf,projections = false) = projections ? dsum(X,Y,projections)[[1,3]] : dsum(X,Y)
coproduct(X::CohSheaf,Y::CohSheaf,projections = false) = projections ? dsum(X,Y,projections)[[1,2]] : dsum(X,Y)

#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

"""
    tensor_product(X::CohSheaf{T,G}, Y::CohSheaf{T,G}) where {T,G}

Return the tensor product of equivariant coherent sheaves.
"""
function tensor_product(X::CohSheaf{T,G}, Y::CohSheaf{T,G}) where {T,G}
    @assert parent(X) == parent(Y) "Mismatching parents"
    return CohSheaf{T,G}(parent(X), [x⊗y for (x,y) ∈ zip(stalks(X), stalks(Y))])
end

"""
    tensor_product(f::CohSheafMorphism{T,G}, g::CohSheafMorphism{T,G}) where {T,G}

Return the tensor product of morphisms of equivariant coherent sheaves.
"""
function tensor_product(f::CohSheafMorphism{T,G}, g::CohSheafMorphism{T,G}) where {T,G}
    dom = tensor_product(domain(f), domain(g))
    codom = tensor_product(codomain(f), codomain(g))
    mors = [tensor_product(m,n) for (m,n) ∈ zip(f.m,g.m)]
    return CohSheafMorphism{T,G}(dom,codom,mors)
end

"""
    one(C::CohSheaves{T,G}) where {T,G}

Return the one object in ``C``.
"""
function one(C::CohSheaves{T,G}) where {T,G}
    return CohSheaf{T,G}(C,[one(RepresentationCategory(H,base_ring(C))) for H ∈ C.orbit_stabilizers])
end

#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::CohSheafMorphism{T,G}, g::CohSheafMorphism{T,G}) where {T,G}
    dom = domain(f)
    codom = codomain(f)
    mors = [compose(m,n) for (m,n) ∈ zip(f.m,g.m)]
    return CohSheafMorphism{T,G}(dom,codom,mors)
end

#-----------------------------------------------------------------
#   Simple Objects
#-----------------------------------------------------------------
"""
    simples(C::CohSheaves{T,G}) where {T,G}

Return the simple objects of ``C``.
"""
function simples(C::CohSheaves{T,G}) where {T,G}

    simple_objects = CohSheaf[]

    #zero modules
    zero_mods = [zero(RepresentationCategory(H,C.base_ring)) for H ∈ C.orbit_stabilizers]

    for k ∈ 1:length(C.orbit_stabilizers)

        #Get simple objects from the corresponding representation categories
        RepH = RepresentationCategory(C.orbit_stabilizers[k], C.base_ring)

        RepH_simples = simples(RepH)
        for i ∈ 1:length(RepH_simples)
            Hsimple_sheaves = [CohSheaf(C,[k == j ? RepH_simples[i] : zero_mods[j] for j ∈ 1:length(C.orbit_stabilizers)])]
            simple_objects = [simple_objects; Hsimple_sheaves]
        end
    end
    return simple_objects
end

"""
    decompose(X::CohSheaf)

Decompose ``X`` into a direct sum of simple objects with multiplicity.
"""
function decompose(X::CohSheaf)
    ret = []
    C = parent(X)
    zero_mods = [zero(RepresentationCategory(H,C.base_ring)) for H ∈ C.orbit_stabilizers]

    for k ∈ 1:length(C.orbit_reps)
        X_H_facs = decompose(stalks(X)[k])

        ret = [ret; [(CohSheaf(C,[k == j ? Y : zero_mods[j] for j ∈ 1:length(C.orbit_reps)]),d) for (Y,d) ∈ X_H_facs]]
    end
    return ret
end
#-----------------------------------------------------------------
#   Hom Spaces
#-----------------------------------------------------------------

struct CohSfHomSpace{T,G} <: HomSpace{T}
    X::CohSheaf{T,G}
    Y::CohSheaf{T,G}
    basis::Vector{CohSheafMorphism{T,G}}
    parent::VectorSpaces{T}
end

"""
    Hom(X::CohSheaf{T,G}, Y::CohSheaf{T,G}) where {T,G}

Return Hom(``X,Y``) as a vector space.
"""
function Hom(X::CohSheaf{T,G}, Y::CohSheaf{T,G}) where {T,G}
    @assert parent(X) == parent(Y) "Missmatching parents"

    b = CohSheafMorphism{T,G}[]
    for i ∈ 1:length(stalks(X))
        Hxy = Hom(stalks(X)[i],stalks(Y)[i])

        for ρ ∈ basis(Hxy)
            reps = [zero(Hxy) for j ∈ 1:length(stalks(X))]
            reps[i] = ρ
            b = [b; CohSheafMorphism(X,Y,reps)]
        end
    end
    return CohSfHomSpace{T,G}(X,Y,b,VectorSpaces(base_ring(X)))
end

zero(H::CohSfHomSpace) = zero_morphism(H.X,H.Y)


#-----------------------------------------------------------------
#   Pretty Printing
#-----------------------------------------------------------------

function show(io::IO, C::CohSheaves)
    print(io, "Category of equivariant coherent sheaves on $(C.GSet.seeds) over $(C.base_ring).")
end

function show(io::IO, X::CohSheaf)
    print(io, "Equivariant choherent sheaf on $(X.parent.GSet.seeds) over $(base_ring(X)).")
end

function show(io::IO, X::CohSheafMorphism)
    print(io, "Morphism of equivariant choherent sheaves on $(domain(X).parent.GSet.seeds) over $(base_ring(X)).")
end


#-----------------------------------------------------------------
#   Functors
#-----------------------------------------------------------------

struct PullbackFunctor{T} <: Functor
    domain::T
    codomain::T
    obj_map
    mor_map
end

pullb_obj_map(CY,CX,X,f) = CohSheaf(CX, [restriction(stalk(X,f(x)), H) for (x,H) ∈ zip(CX.orbit_reps, CX.orbit_stabilizers)])

function pullb_mor_map(CY,CX,m,f)
    dom = pullb_obj_map(CY,CX,domain(m),f)
    codom = pullb_obj_map(CY,CX,codomain(m),f)
    maps = [restriction(m.m[orbit_index(CY,f(x))], H) for (x,H) ∈ zip(CX.orbit_reps, CX.orbit_stabilizers)]
    CohSheafMorphism(dom, codom, maps)
end

function Pullback(CY::CohSheaves{T,G}, CX::CohSheaves{T,G}, f::Function) where {T,G}
    @assert isequivariant(CX.GSet, CY.GSet, f) "Map not equivariant"

    obj_map = X -> pullb_obj_map(CY,CX,X,f)
    mor_map = m -> pullb_mor_map(CY,CX,m,f)

    return PullbackFunctor{CohSheaves{T,G}}(CY, CX, obj_map, mor_map)
end

struct PushforwardFunctor{T} <: Functor
    domain::T
    codomain::T
    obj_map
    mor_map
end

function pushf_obj_map(CX,CY,X,f)
    stlks = [zero(RepresentationCategory(H,base_ring(CY))) for H ∈ CY.orbit_stabilizers]
    for i ∈ 1:length(CY.orbit_reps)
        y = CY.orbit_reps[i]
        Gy = CY.orbit_stabilizers[i]

        if length([x for x ∈ CX.GSet.seeds if f(x) == y]) == 0 continue end
        fiber = gset(Gy, CX.GSet.action_function, [x for x ∈ CX.GSet.seeds if f(x) == y])

        orbit_reps = [O.seeds[1] for O ∈ orbits(fiber)]

        for j ∈ 1:length(orbit_reps)
            stlks[i] = dsum(stlks[i], induction(stalk(X,orbit_reps[j]), CY.orbit_stabilizers[i]))
        end
    end
    return CohSheaf(CY, stlks)
end

function pushf_mor_map(CX,CY,m,f)
    mor = [zero_morphism(s,r) for (s,r) ∈ zip(stalks(domain(m)),stalks(codomain(m)))]
    for i ∈ 1:length(CY.orbit_reps)
        y = CY.orbit_reps[i]
        Gy = CY.orbit_stabilizers[i]

        if length([x for x ∈ CX.GSet.seeds if f(x) == y]) == 0 continue end
        fiber = gset(Gy, CX.GSet.action_function, [x for x ∈ CX.GSet.seeds if f(x) == y])

        orbit_reps = [O.seeds[1] for O ∈ orbits(fiber)]

        mor[i] = dsum([induction(mor[i],Gy); [induction(m.m[orbit_index(CX,y)], Gy) for y ∈ orbit_reps]]...)
    end
    return CohSheafMorphism(pushf_obj_map(CX,CY,domain(m),f), pushf_obj_map(CX,CY,codomain(m),f), mor)
end

function Pushforward(CX::CohSheaves{T,G}, CY::CohSheaves{T,G}, f::Function) where {T,G}
    @assert isequivariant(CX.GSet, CY.GSet, f) "Map not equivariant"

    return PushforwardFunctor(CX,CY,X -> pushf_obj_map(CX,CY,X,f),m -> pushf_mor_map(CX,CY,m,f))
end

#dummy
function isequivariant(X::GSet, Y::GSet, f::Function)
    true
end
