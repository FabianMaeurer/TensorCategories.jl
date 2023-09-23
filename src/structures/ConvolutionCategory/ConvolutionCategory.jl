struct ConvolutionCategory <: Category
    group::GAPGroup
    base_ring::Field
    GSet::GSet
    squaredGSet::GSet
    cubedGSet::GSet
    squaredCoh::CohSheaves
    cubedCoh::CohSheaves
    projectors::Vector
end

struct ConvolutionObject <: Object
    sheaf::CohSheafObject
    parent::ConvolutionCategory
end

struct ConvolutionMorphism <: Morphism
    domain::ConvolutionObject
    codomain::ConvolutionObject
    m::CohSheafMorphism
end
is_tensor(::ConvolutionCategory) = true
is_fusion(C::ConvolutionCategory) = mod(order(C.group),characteristic(base_ring(C))) != 0
"""
    ConvolutionCategory(X::GSet, K::Field)

Return the category of equivariant coherent sheaves on ``X`` with convolution product.
"""
function ConvolutionCategory(X::GSet, K::Field)
    G = X.group
    sqX = gset(G,(x,g) -> Tuple(X.action_function(xi,g) for xi ∈ x), [(x,y) for x ∈ X.seeds, y ∈ X.seeds][:])
    cuX = gset(G,(x,g) -> Tuple(X.action_function(xi,g) for xi ∈ x), [(x,y,z) for x ∈ X.seeds, y ∈ X.seeds, z ∈ X.seeds][:])
    sqCoh = CohSheaves(sqX,K)
    cuCoh = CohSheaves(cuX,K)

    p12 = x -> (x[1],x[2])
    p13 = x -> (x[1],x[3])
    p23 = x -> (x[2],x[3])

    P12 = Pullback(sqCoh, cuCoh, p12)
    P13 = Pushforward(cuCoh, sqCoh, p13)
    P23 = Pullback(sqCoh, cuCoh, p23)
    return ConvolutionCategory(G,K,X,sqX,cuX,sqCoh,cuCoh, [P12, P13, P23])
end

"""
    ConvolutionCategory(X, K::Field)

Return the category of coherent sheaves on ``X`` with convolution product.
"""
function ConvolutionCategory(X, K::Field)
    G = symmetric_group(1)
    return ConvolutionCategory(gset(G,X), K)
end

function ConvolutionCategory(X, G::GAPGroup, K::Field)
    ConvolutionCategory(gset(G,X), K)
end

Morphism(D::ConvolutionObject, C::ConvolutionObject, m:: CohSheafMorphism) = ConvolutionMorphism(D,C,m)

function Base.hash(C::ConvolutionCategory, h::UInt)
    hash((getfield(C, s) for s ∈ fieldnames(ConvolutionCategory)), h)
end

function Base.hash(X::ConvolutionObject, h::UInt)
    hash((getfield(σ, s) for s ∈ fieldnames(ConvolutionObject)), h)
end

#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------

"""
    is_semisimple(C::ConvolutionCategory)

Check whether ``C`` semisimple.
"""
is_semisimple(C::ConvolutionCategory) = gcd(order(C.group), characteristic(base_ring(C))) == 1

"""
    orbit_stabilizers(C::ConvolutionCategory)

Return the stabilizers of representatives of the orbits.
"""
orbit_stabilizers(C::ConvolutionCategory) = C.squaredCoh.orbit_stabilizers
orbit_index(X::ConvolutionObject, y) = orbit_index(X.sheaf, y)
orbit_index(X::ConvolutionCategory, y) = orbit_index(X.squaredCoh, y)

"""
    stalks(X::ConvolutionObject)

Return the stalks of the sheaf ``X``.
"""
stalks(X::ConvolutionObject) = stalks(X.sheaf)
stalk(X::ConvolutionObject, x) = stalk(X.sheaf,x)

==(X::ConvolutionObject,Y::ConvolutionObject) = X.sheaf == Y.sheaf

==(f::ConvolutionMorphism, g::ConvolutionMorphism) = f.m == f.m

"""
    is_isomorphic(X::ConvolutionObject, Y::ConvolutionObject)

Check whether ``X``and ``Y``are isomorphic. Return the isomorphism if true.
"""
function is_isomorphic(X::ConvolutionObject, Y::ConvolutionObject)
    b, iso = is_isomorphic(X.sheaf, Y.sheaf)
    if !b return false, nothing end
    return true, ConvolutionMorphism(X,Y,iso)
end

id(X::ConvolutionObject) = ConvolutionMorphism(X,X,id(X.sheaf))

function associator(X::ConvolutionObject, Y::ConvolutionObject, Z::ConvolutionObject)
    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)
    if dom == zero(parent(X)) return zero_morphism(dom,cod) end
    S = simples(parent(X))
    #@show codomain(decompose_morphism(cod,S)[1]) == codomain(decompose_morphism(dom,S)[1])
    return inv(direct_sum_decomposition(cod,S)[2]) ∘ direct_sum_decomposition(dom,S)[2]
end
#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    direct_sum(X::ConvolutionObject, Y::ConvolutionObject, morphisms::Bool = false)

documentation
"""
function direct_sum(X::ConvolutionObject, Y::ConvolutionObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Z,ix,px = direct_sum(X.sheaf,Y.sheaf)
    Z = ConvolutionObject(Z,parent(X))


    ix = [ConvolutionMorphism(x,Z,i) for (x,i) ∈ zip([X,Y],ix)]
    px = [ConvolutionMorphism(Z,x,p) for (x,p) ∈ zip([X,Y],px)]
    return Z, ix, px
end

"""
    direct_sum(f::ConvolutionMorphism, g::ConvolutionMorphism)

Return the direct sum of morphisms of coherent sheaves (with convolution product).
"""
function direct_sum(f::ConvolutionMorphism, g::ConvolutionMorphism)
    dom = direct_sum(domain(f), domain(g))
    codom = direct_sum(codomain(f), codomain(g))
    m = direct_sum(f.m,g.m)
    return ConvolutionMorphism(dom,codom,m)
end


"""
    zero(C::ConvolutionCategory)

Return the zero object in Conv(``X``).
"""
zero(C::ConvolutionCategory) = ConvolutionObject(zero(C.squaredCoh),C)


#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

function kernel(f::ConvolutionMorphism)
    K,k = kernel(f.m)
    Conv_K = ConvolutionObject(K,parent(f))
    return ConvolutionObject(K,parent(domain(f))), Morphism(Conv_K, domain(f), k)
end

function cokernel(f::ConvolutionMorphism)
    C,c = cokernel(f.m)
    Conv_C = ConvolutionObject(C,parent(f))
    return ConvolutionObject(C, parent(domain(f))), Morphism(codomain(f), Conv_C, c)
end

#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

"""
    tensor_product(X::ConvolutionObject, Y::ConvolutionObject)

Return the convolution product of coherent sheaves.
"""
function tensor_product(X::ConvolutionObject, Y::ConvolutionObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    p12,p13,p23 = parent(X).projectors

    return ConvolutionObject(p13(p12(X.sheaf)⊗p23(Y.sheaf)),parent(X))
end

"""
    tensor_product(f::ConvolutionMorphism, g::ConvolutionMorphism)

Return the convolution product of morphisms of coherent sheaves.
"""
function tensor_product(f::ConvolutionMorphism, g::ConvolutionMorphism)
    dom = domain(f)⊗domain(g)
    codom = codomain(f)⊗codomain(g)

    p12,p13,p23 = parent(domain(f)).projectors
    return ConvolutionMorphism(dom,codom, p13(p12(f.m)⊗p23(g.m)))
end

"""
    one(C::ConvolutionCategory)

Return the one object in Conv(``X``).
"""
function one(C::ConvolutionCategory)
    F = base_ring(C)

    stlks = [zero(RepresentationCategory(H,F)) for H ∈ orbit_stabilizers(C)]
    diag = [(x,x) for x ∈ C.GSet.seeds]

    for i ∈ [orbit_index(C,d) for d ∈ diag]
        stlks[i] = one(RepresentationCategory(orbit_stabilizers(C)[i], F))
    end
    return ConvolutionObject(CohSheafObject(C.squaredCoh, stlks), C)
end

function dual(X::ConvolutionObject)
    orbit_reps = parent(X).squaredCoh.orbit_reps
    GSet = parent(X).squaredCoh.GSet
    perm = [findfirst(e -> e ∈ orbit(GSet, (y,x)), orbit_reps) for (x,y) ∈ orbit_reps]
    reps = [dual(ρ) for ρ ∈ stalks(X)][perm]
    return ConvolutionObject(CohSheafObject(parent(X.sheaf), reps), parent(X))
end

spherical(X::ConvolutionObject) = id(X)
#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::ConvolutionMorphism,g::ConvolutionMorphism)
    return ConvolutionMorphism(domain(f),codomain(g),compose(f.m,g.m))
end

function zero_morphism(X::ConvolutionObject, Y::ConvolutionObject)
    return ConvolutionMorphism(X,Y,zero_morphism(X.sheaf,Y.sheaf))
end

function +(f::ConvolutionMorphism, g::ConvolutionMorphism)
    Morphism(domain(f),codomain(f), f.m + g.m)
end

function *(x, f::ConvolutionMorphism)
    Morphism(domain(f), codomain(f), x * f.m)
end

function matrices(f::ConvolutionMorphism)
    matrices(f.m)
end

matrix(f::ConvolutionMorphism) = matrix(f.m)

function inv(f::ConvolutionMorphism)
    return Morphism(codomain(f), domain(f), inv(f.m))
end

left_inverse(f::ConvolutionMorphism) = Morphism(codomain(f),domain(f), left_inverse(f.m))
#-----------------------------------------------------------------
#   Simple Objects
#-----------------------------------------------------------------

"""
    simples(C::ConvolutionCategory)

Return a list of simple objects in Conv(``X``).
"""
#= @memoize Dict =# function simples(C::ConvolutionCategory)
    return [ConvolutionObject(sh,C) for sh ∈ simples(C.squaredCoh)]
end

"""
    decompose(X::ConvolutionObject)

Decompose ``X`` into a direct sum of simple objects with multiplicities.
"""
function decompose(X::ConvolutionObject)
    facs = decompose(X.sheaf)
    return [(ConvolutionObject(sh,parent(X)),d) for (sh,d) ∈ facs]
end

function is_simple(X::ConvolutionObject)
    non_zero_stalks = [s for s ∈ stalks(X) if s != zero(parent(s))]
    if length(non_zero_stalks) != 1
        return false
    end
    return is_simple(non_zero_stalks[1])
end 

#-----------------------------------------------------------------
#   Hom Space
#-----------------------------------------------------------------

struct ConvHomSpace <: AbstractHomSpace
    X::ConvolutionObject
    Y::ConvolutionObject
    basis::Vector{ConvolutionMorphism}
    parent::VectorSpaces
end

"""
    Hom(X::ConvolutionObject, Y::ConvolutionObject)

Return Hom(``X,Y``) as a vector space.
"""
function Hom(X::ConvolutionObject, Y::ConvolutionObject)
    @assert parent(X) == parent(Y) "Missmatching parents"
    b = basis(Hom(X.sheaf,Y.sheaf))
    conv_b = [ConvolutionMorphism(X,Y,m) for m ∈ b]

    return ConvHomSpace(X,Y,conv_b, VectorSpaces(base_ring(X)))
end

zero(H::ConvHomSpace) = zero_morphism(H.X,H.Y)


#-----------------------------------------------------------------
#   pretty printing
#-----------------------------------------------------------------

function show(io::IO, C::ConvolutionCategory)
    print(io, """Convolution category over G-set with $(length(C.GSet)) elements.""")
end

function show(io::IO, X::ConvolutionObject)
    print(io, "Object in $(parent(X))")
end

function show(io::IO, f::ConvolutionMorphism)
    print(io, "Morphism in $(parent(domain(f)))")
end
