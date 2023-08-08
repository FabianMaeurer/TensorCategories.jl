struct CohSheaves <: Category
    group::GAPGroup
    base_ring::Field
    GSet::GSet
    orbit_reps
    orbit_stabilizers
end

struct CohSheafObject <: Object
    parent::CohSheaves
    stalks::Vector{GroupRepresentation}
end

struct CohSheafMorphism <: Morphism
    domain::CohSheafObject
    codomain::CohSheafObject
    m::Vector{GroupRepresentationMorphism}
end

is_multitensor(::CohSheaves) = true
is_multifusion(C::CohSheaves) = mod(order(C.group),characteristic(base_ring(C))) != 0

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
    return CohSheaves(G, F, X, orbit_reps, orbit_stabilizers)
end

"""
    CohSheaves(X, F::Field)

The category of coherent sheafes on ``X``.
"""
function CohSheaves(X,F::Field)
    G = symmetric_group(1)
    return CohSheaves(gset(G,X), F)
end

Morphism(X::CohSheafObject, Y::CohSheafObject, m::Vector) = CohSheafMorphism(X,Y,m)


#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------
"""
    is_semisimple(C::CohSheaves)

Return whether ``C``is semisimple.
"""
is_semisimple(C::CohSheaves) = gcd(order(C.group), characteristic(base_ring(C))) == 1

"""
    stalks(X::CohSheafObject)

Return the stalks of ``X``.
"""
stalks(X::CohSheafObject) = X.stalks

orbit_stabilizers(Coh::CohSheaves) = Coh.orbit_stabilizers

function orbit_index(X::CohSheafObject, y)
    i = findfirst(x -> y ∈ x, orbits(parent(X).GSet))
end

function orbit_index(C::CohSheaves, y)
    i = findfirst(x -> y ∈ x, orbits(C.GSet))
end

function stalk(X::CohSheafObject,y)
    return stalks(X)[orbit_index(X,y)]
end

"""
    zero(C::CohSheaves)

Return the zero sheaf on the ``G``-set.
"""
zero(C::CohSheaves) = CohSheafObject(C,[zero(RepresentationCategory(H,base_ring(C))) for H ∈ C.orbit_stabilizers])

"""
    zero_morphism(X::CohSheafObject, Y::CohSheafObject)

Return the zero morphism ```0:X → Y```.
"""
zero_morphism(X::CohSheafObject, Y::CohSheafObject) = CohSheafMorphism(X,Y,[zero(Hom(x,y)) for (x,y) ∈ zip(stalks(X),stalks(Y))])

function ==(X::CohSheafObject, Y::CohSheafObject)
    if parent(X) != parent(Y) return false end
    for (s,r) ∈ zip(stalks(X),stalks(Y))
        if s != r return false end
    end
    return true
end

"""
    is_isomorphic(X::CohSheafObject, Y::CohSheafObject)

Check whether ``X``and ``Y`` are isomorphic and the isomorphism if possible.
"""
function is_isomorphic(X::CohSheafObject, Y::CohSheafObject)
    m = GroupRepresentationMorphism[]
    for (s,r) ∈ zip(stalks(X),stalks(Y))
        b, iso = is_isomorphic(s,r)
        if !b return false, nothing end
        m = [m; iso]
    end
    return true, CohSheafMorphism(X,Y,m)
end

==(f::CohSheafMorphism, g::CohSheafMorphism) = f.m == g.m

"""
    id(X::CohSheafObject)

Return the identity on ```X```.
"""
id(X::CohSheafObject) = CohSheafMorphism(X,X,[id(s) for s ∈ stalks(X)])

"""
    associator(X::CohSheafObject, Y::CohSheafObject, Z::CohSheafObject)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
associator(X::CohSheafObject, Y::CohSheafObject, Z::CohSheafObject) = id(X⊗Y⊗Z)

"""
    dual(X::CohSheafObject)

Return the dual object of ```X```.
"""
dual(X::CohSheafObject) = CohSheafObject(parent(X),[dual(s) for s ∈ stalks(X)])

"""
    ev(X::CohSheafObject)

Return the evaluation morphism ```X∗⊗X → 1```.
"""
function ev(X::CohSheafObject)
    dom = dual(X)⊗X
    cod = one(parent(X))
    return CohSheafMorphism(dom,cod, [ev(s) for s ∈ stalks(X)])
end

"""
    coev(X::CohSheafObject)

Return the coevaluation morphism ```1 → X⊗X∗```.
"""
function coev(X::CohSheafObject)
    dom = one(parent(X))
    cod = X ⊗ dual(X)
    return CohSheafMorphism(dom,cod, [coev(s) for s ∈ stalks(X)])
end

"""
    spherical(X::CohSheafObject)

Return the spherical structure isomorphism ```X → X∗∗```.
"""
spherical(X::CohSheafObject) = Morphism(X,X,[spherical(s) for s ∈ stalks(X)])

"""
    braiding(X::cohSheaf, Y::CohSheafObject)

Return the braiding isomoephism ```X⊗Y → Y⊗X```.
"""
braiding(X::CohSheafObject, Y::CohSheafObject) = Morphism(X⊗Y, Y⊗X, [braiding(x,y) for (x,y) ∈ zip(stalks(X),stalks(Y))])

dim(X::CohSheafObject) = sum(dim.(stalks(X)))

function (F::Field)(f::CohSheafMorphism) 
    @assert is_simple(domain(f)) && domain(f) == codomain(f)
    return sum([F(m) for m ∈ f.m if m != zero_morphism(domain(m),codomain(m))])
end

function is_simple(X::CohSheafObject)
    non_zero = [s for s ∈ stalks(X) if s != zero(parent(s))]
    if length(non_zero) != 1
        return false
    end
    return is_simple(non_zero[1])
end

#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    direct_sum(X::CohSheafObject, Y::CohSheafObject)

Return the direct sum of sheaves. Return also the inclusion and projection.
"""
function direct_sum(X::CohSheafObject, Y::CohSheafObject)
    sums = [direct_sum(x,y) for (x,y) ∈ zip(stalks(X), stalks(Y))]
    Z = CohSheafObject(parent(X), [s[1] for s ∈ sums])

    ix = [CohSheafMorphism(x,Z,[s[2][i] for s ∈ sums]) for (x,i) ∈ zip([X,Y],1:2)]
    px = [CohSheafMorphism(Z,x,[s[3][i] for s ∈ sums]) for (x,i) ∈ zip([X,Y],1:2)]
    return Z,ix,px
end

"""
    direct_sum(f::CohSheafMorphism, g::CohSheafMorphism)

Return the direct sum of morphisms of sheaves.
"""
function direct_sum(f::CohSheafMorphism, g::CohSheafMorphism)
    dom = direct_sum(domain(f), domain(g))[1]
    codom = direct_sum(codomain(f), codomain(g))[2]
    mors = [direct_sum(m,n) for (m,n) ∈ zip(f.m,g.m)]
    return CohSheafMorphism(dom,codom, mors)
end


#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

"""
    kernel(f::CohSheafMorphism)

Return a tuple ```(K,k)``` where ```K``` is the kernel object and ```k``` is the inclusion.
"""
function kernel(f::CohSheafMorphism)
    kernels = [kernel(g) for g ∈ f.m]
    K = CohSheafObject(parent(domain(f)), [k for (k,_) ∈ kernels])
    return K, Morphism(K, domain(f), [m for (_,m) ∈ kernels])
end

"""
    cokernel(f::CohSheafMorphism)

Return a tuple ```(C,c)``` where ```C``` is the kernel object and ```c``` is the projection.
"""
function cokernel(f::CohSheafMorphism)
    cokernels = [cokernel(g) for g ∈ f.m]
    C = CohSheafObject(parent(domain(f)), [c for (c,_) ∈ cokernels])
    return C, Morphism(codomain(f), C, [m for (_,m) ∈ cokernels])
end

#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

"""
    tensor_product(X::CohSheafObject, Y::CohSheafObject)

Return the tensor product of equivariant coherent sheaves.
"""
function tensor_product(X::CohSheafObject, Y::CohSheafObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    return CohSheafObject(parent(X), [x⊗y for (x,y) ∈ zip(stalks(X), stalks(Y))])
end

"""
    tensor_product(f::CohSheafMorphism, g::CohSheafMorphism)

Return the tensor product of morphisms of equivariant coherent sheaves.
"""
function tensor_product(f::CohSheafMorphism, g::CohSheafMorphism)
    dom = tensor_product(domain(f), domain(g))
    codom = tensor_product(codomain(f), codomain(g))

    mors = [tensor_product(m,n) for (m,n) ∈ zip(f.m,g.m)]
    return CohSheafMorphism(dom,codom,mors)
end

"""
    one(C::CohSheaves)

Return the one object in ``C``.
"""
function one(C::CohSheaves)
    return CohSheafObject(C,[one(RepresentationCategory(H,base_ring(C))) for H ∈ C.orbit_stabilizers])
end

#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

"""
    compose(f::CohSheafMorphism, g::CohSheafMorphism)

Return the composition ```g∘f```.
"""
function compose(f::CohSheafMorphism, g::CohSheafMorphism)
    dom = domain(f)
    codom = codomain(f)
    mors = [compose(m,n) for (m,n) ∈ zip(f.m,g.m)]
    return CohSheafMorphism(dom,codom,mors)
end

function +(f::CohSheafMorphism, g::CohSheafMorphism)
    #@assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    return Morphism(domain(f), codomain(f), [fm + gm for (fm,gm) ∈ zip(f.m,g.m)])
end

function *(x,f::CohSheafMorphism)
    Morphism(domain(f),codomain(f),x .* f.m)
end

matrices(f::CohSheafMorphism) = matrix.(f.m)

"""
    inv(f::CohSheafMorphism)

Retrn the inverse morphism of ```f```.
"""
function inv(f::CohSheafMorphism)
    return Morphism(codomain(f), domain(f), [inv(g) for g in f.m])
end


function matrix(f::CohSheafMorphism)
    diagonal_matrix(matrices(f))
end

function left_inverse(f::CohSheafMorphism)
    return Morphism(codomain(f),domain(f), [left_inverse(g) for g ∈ f.m])
end
#-----------------------------------------------------------------
#   Simple Objects
#-----------------------------------------------------------------
"""
    simples(C::CohSheaves)

Return the simple objects of ``C``.
"""
function simples(C::CohSheaves)

    simple_objects = CohSheafObject[]

    #zero modules
    zero_mods = [zero(RepresentationCategory(H,C.base_ring)) for H ∈ C.orbit_stabilizers]

    for k ∈ 1:length(C.orbit_stabilizers)

        #Get simple objects from the corresponding representation categories
        RepH = RepresentationCategory(C.orbit_stabilizers[k], C.base_ring)

        RepH_simples = simples(RepH)
        for i ∈ 1:length(RepH_simples)
            Hsimple_sheaves = [CohSheafObject(C,[k == j ? RepH_simples[i] : zero_mods[j] for j ∈ 1:length(C.orbit_stabilizers)])]
            simple_objects = [simple_objects; Hsimple_sheaves]
        end
    end
    return simple_objects
end

"""
    decompose(X::CohSheafObject)

Decompose ``X`` into a direct sum of simple objects with multiplicity.
"""
function decompose(X::CohSheafObject)
    ret = []
    C = parent(X)
    zero_mods = [zero(RepresentationCategory(H,C.base_ring)) for H ∈ C.orbit_stabilizers]

    for k ∈ 1:length(C.orbit_reps)
        X_H_facs = decompose(stalks(X)[k])

        ret = [ret; [(CohSheafObject(C,[k == j ? Y : zero_mods[j] for j ∈ 1:length(C.orbit_reps)]),d) for (Y,d) ∈ X_H_facs]]
    end
    return ret
end
#-----------------------------------------------------------------
#   Hom Spaces
#-----------------------------------------------------------------

struct CohSfCategoryHomSpace <: AbstractCategoryHomSpace
    X::CohSheafObject
    Y::CohSheafObject
    basis::Vector{CohSheafMorphism}
    parent::VectorSpaces
end

"""
    Hom(X::CohSheafObject, Y::CohSheafObject)

Return Hom(``X,Y``) as a vector space.
"""
function Hom(X::CohSheafObject, Y::CohSheafObject)
    @assert parent(X) == parent(Y) "Missmatching parents"

    b = CohSheafMorphism[]
    H = [Hom(stalks(X)[i],stalks(Y)[i]) for i ∈ 1:length(stalks(X))]
    for i ∈ 1:length(stalks(X))
        for ρ ∈ basis(H[i])
            reps = [zero(H[j]) for j ∈ 1:length(stalks(X))]
            reps[i] = ρ
            b = [b; CohSheafMorphism(X,Y,reps)]
        end
    end
    return CohSfCategoryHomSpace(X,Y,b,VectorSpaces(base_ring(X)))
end

zero(H::CohSfCategoryHomSpace) = zero_morphism(H.X,H.Y)


#-----------------------------------------------------------------
#   Pretty Printing
#-----------------------------------------------------------------

function show(io::IO, C::CohSheaves)
    print(io, "Category of equivariant coherent sheaves on $(C.GSet.seeds) over $(C.base_ring)")
end

function show(io::IO, X::CohSheafObject)
    print(io, "Equivariant choherent sheaf on $(X.parent.GSet.seeds) over $(base_ring(X))")
end

function show(io::IO, X::CohSheafMorphism)
    print(io, "Morphism of equivariant choherent sheaves on $(domain(X).parent.GSet.seeds) over $(base_ring(X))")
end


#-----------------------------------------------------------------
#   Functors
#-----------------------------------------------------------------

struct PullbackFunctor <: AbstractFunctor
    domain::Category
    codomain::Category
    obj_map
    mor_map
end

pullb_obj_map(CY,CX,X,f) = CohSheafObject(CX, [restriction(stalk(X,f(x)), H) for (x,H) ∈ zip(CX.orbit_reps, CX.orbit_stabilizers)])

function pullb_mor_map(CY,CX,m,f)
    dom = pullb_obj_map(CY,CX,domain(m),f)
    codom = pullb_obj_map(CY,CX,codomain(m),f)
    maps = [restriction(m.m[orbit_index(CY,f(x))], H) for (x,H) ∈ zip(CX.orbit_reps, CX.orbit_stabilizers)]
    CohSheafMorphism(dom, codom, maps)
end

"""
    Pullback(C::CohSheaves, D::CohSheaves, f::Function)

Return the pullback functor ```C → D``` defined by the ```G```-set map ```f::X → Y```.
"""
function Pullback(CY::CohSheaves, CX::CohSheaves, f::Function)
    @assert isequivariant(CX.GSet, CY.GSet, f) "Map not equivariant"

    obj_map = X -> pullb_obj_map(CY,CX,X,f)
    mor_map = m -> pullb_mor_map(CY,CX,m,f)

    return PullbackFunctor(CY, CX, obj_map, mor_map)
end

struct PushforwardFunctor <: AbstractFunctor
    domain::Category
    codomain::Category
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
            stlks[i] = direct_sum(stlks[i], induction(stalk(X,orbit_reps[j]), CY.orbit_stabilizers[i]))[1]
        end
    end
    return CohSheafObject(CY, stlks)
end

function pushf_mor_map(CX,CY,m,f)
    mor = GroupRepresentationMorphism[]

    for i ∈ 1:length(CY.orbit_reps)
        y = CY.orbit_reps[i]
        Gy = CY.orbit_stabilizers[i]

        if length([x for x ∈ CX.GSet.seeds if f(x) == y]) == 0 continue end
        fiber = gset(Gy, CX.GSet.action_function, [x for x ∈ CX.GSet.seeds if f(x) == y])

        orbit_reps = [O.seeds[1] for O ∈ orbits(fiber)]

        mor = [mor; direct_sum([induction(m.m[orbit_index(CX,y)], Gy) for y ∈ orbit_reps]...)]
    end
    return CohSheafMorphism(pushf_obj_map(CX,CY,domain(m),f), pushf_obj_map(CX,CY,codomain(m),f), mor)
end

"""
    Pushforward(C::CohSheaves, D::CohSheaves, f::Function)

Return the push forward functor ```C → D``` defined by the ```G```-set map ```f::X → Y```.
"""
function Pushforward(CX::CohSheaves, CY::CohSheaves, f::Function)
    @assert isequivariant(CX.GSet, CY.GSet, f) "Map not equivariant"

    return PushforwardFunctor(CX,CY,X -> pushf_obj_map(CX,CY,X,f),m -> pushf_mor_map(CX,CY,m,f))
end

#dummy
function isequivariant(X::GSet, Y::GSet, f::Function)
    true
end

function show(io::IO,F::PushforwardFunctor)
    print(io, "Pushforward functor from $(domain(F)) to $(codomain(F))")
end

function show(io::IO,F::PullbackFunctor)
    print(io, "Pullback functor from $(domain(F)) to $(codomain(F))")
end
