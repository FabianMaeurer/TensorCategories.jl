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

struct ConvolutionCategoryObject <: CategoryObject
    sheaf::CohSheaf
    parent::ConvolutionCategory
end

struct ConvolutionCategoryMorphism <: CategoryMorphism
    domain::ConvolutionCategoryObject
    codomain::ConvolutionCategoryObject
    m::CohSheafCategoryMorphism
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

Morphism(D::ConvolutionCategoryObject, C::ConvolutionCategoryObject, m:: CohSheafCategoryMorphism) = ConvolutionCategoryMorphism(D,C,m)

function Base.hash(C::ConvolutionCategory, h::UInt)
    hash((getfield(C, s) for s ∈ fieldnames(ConvolutionCategory)), h)
end

function Base.hash(X::ConvolutionCategoryObject, h::UInt)
    hash((getfield(σ, s) for s ∈ fieldnames(ConvolutionCategoryObject)), h)
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
orbit_index(X::ConvolutionCategoryObject, y) = orbit_index(X.sheaf, y)
orbit_index(X::ConvolutionCategory, y) = orbit_index(X.squaredCoh, y)

"""
    stalks(X::ConvolutionCategoryObject)

Return the stalks of the sheaf ``X``.
"""
stalks(X::ConvolutionCategoryObject) = stalks(X.sheaf)
stalk(X::ConvolutionCategoryObject, x) = stalk(X.sheaf,x)

==(X::ConvolutionCategoryObject,Y::ConvolutionCategoryObject) = X.sheaf == Y.sheaf

==(f::ConvolutionCategoryMorphism, g::ConvolutionCategoryMorphism) = f.m == f.m

"""
    is_isomorphic(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject)

Check whether ``X``and ``Y``are isomorphic. Return the isomorphism if true.
"""
function is_isomorphic(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject)
    b, iso = is_isomorphic(X.sheaf, Y.sheaf)
    if !b return false, nothing end
    return true, ConvolutionCategoryMorphism(X,Y,iso)
end

id(X::ConvolutionCategoryObject) = ConvolutionCategoryMorphism(X,X,id(X.sheaf))

function associator(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject, Z::ConvolutionCategoryObject)
    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)
    S = simples(parent(X))
    #@show codomain(decompose_morphism(cod,S)[1]) == codomain(decompose_morphism(dom,S)[1])
    return inv(decompose_morphism(cod,S)[1])∘decompose_morphism(dom,S)[1]
end
#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    direct_sum(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject, morphisms::Bool = false)

documentation
"""
function direct_sum(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject, morphisms::Bool = false)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Z,ix,px = direct_sum(X.sheaf,Y.sheaf,true)
    Z = ConvolutionCategoryObject(Z,parent(X))

    if !morphisms return Z end

    ix = [ConvolutionCategoryMorphism(x,Z,i) for (x,i) ∈ zip([X,Y],ix)]
    px = [ConvolutionCategoryMorphism(Z,x,p) for (x,p) ∈ zip([X,Y],px)]
    return Z, ix, px
end

"""
    direct_sum(f::ConvolutionCategoryMorphism, g::ConvolutionCategoryMorphism)

Return the direct sum of morphisms of coherent sheaves (with convolution product).
"""
function direct_sum(f::ConvolutionCategoryMorphism, g::ConvolutionCategoryMorphism)
    dom = direct_sum(domain(f), domain(g))
    codom = direct_sum(codomain(f), codomain(g))
    m = direct_sum(f.m,g.m)
    return ConvolutionCategoryMorphism(dom,codom,m)
end

product(X::ConvolutionCategoryObject,Y::ConvolutionCategoryObject,projections::Bool = false) = projections ? direct_sum(X,Y,projections)[[1,3]] : direct_sum(X,Y)
coproduct(X::ConvolutionCategoryObject,Y::ConvolutionCategoryObject,projections::Bool = false) = projections ? direct_sum(X,Y,projections)[[1,2]] : direct_sum(X,Y)

"""
    zero(C::ConvolutionCategory)

Return the zero object in Conv(``X``).
"""
zero(C::ConvolutionCategory) = ConvolutionCategoryObject(zero(C.squaredCoh),C)


#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

function kernel(f::ConvolutionCategoryMorphism)
    K,k = kernel(f.m)
    Conv_K = ConvolutionCategoryObject(K,parent(f))
    return ConvolutionCategoryObject(K,parent(domain(f))), Morphism(Conv_K, domain(f), k)
end

function cokernel(f::ConvolutionCategoryMorphism)
    C,c = cokernel(f.m)
    Conv_C = ConvolutionCategoryObject(C,parent(f))
    return ConvolutionCategoryObject(C, parent(domain(f))), Morphism(codomain(f), Conv_C, c)
end

#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

"""
    tensor_product(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject)

Return the convolution product of coherent sheaves.
"""
function tensor_product(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    p12,p13,p23 = parent(X).projectors

    return ConvolutionCategoryObject(p13(p12(X.sheaf)⊗p23(Y.sheaf)),parent(X))
end

"""
    tensor_product(f::ConvolutionCategoryMorphism, g::ConvolutionCategoryMorphism)

Return the convolution product of morphisms of coherent sheaves.
"""
function tensor_product(f::ConvolutionCategoryMorphism, g::ConvolutionCategoryMorphism)
    dom = domain(f)⊗domain(g)
    codom = codomain(f)⊗codomain(g)

    p12,p13,p23 = parent(domain(f)).projectors
    return ConvolutionCategoryMorphism(dom,codom, p13(p12(f.m)⊗p23(g.m)))
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
    return ConvolutionCategoryObject(CohSheaf(C.squaredCoh, stlks), C)
end

function dual(X::ConvolutionCategoryObject)
    orbit_reps = parent(X).squaredCoh.orbit_reps
    GSet = parent(X).squaredCoh.GSet
    perm = [findfirst(e -> e ∈ orbit(GSet, (y,x)), orbit_reps) for (x,y) ∈ orbit_reps]
    reps = [dual(ρ) for ρ ∈ stalks(X)][perm]
    return ConvolutionCategoryObject(CohSheaf(parent(X.sheaf), reps), parent(X))
end

spherical(X::ConvolutionCategoryObject) = id(X)
#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::ConvolutionCategoryMorphism,g::ConvolutionCategoryMorphism)
    return ConvolutionCategoryMorphism(domain(f),codomain(g),compose(f.m,g.m))
end

function zero_morphism(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject)
    return ConvolutionCategoryMorphism(X,Y,zero_morphism(X.sheaf,Y.sheaf))
end

function +(f::ConvolutionCategoryMorphism, g::ConvolutionCategoryMorphism)
    Morphism(domain(f),codomain(f), f.m + g.m)
end

function *(x, f::ConvolutionCategoryMorphism)
    Morphism(domain(f), codomain(f), x * f.m)
end

function matrices(f::ConvolutionCategoryMorphism)
    matrices(f.m)
end

matrix(f::ConvolutionCategoryMorphism) = matrix(f.m)

function inv(f::ConvolutionCategoryMorphism)
    return Morphism(codomain(f), domain(f), inv(f.m))
end

left_inverse(f::ConvolutionCategoryMorphism) = Morphism(codomain(f),domain(f), left_inverse(f.m))
#-----------------------------------------------------------------
#   Simple CategoryObjects
#-----------------------------------------------------------------

"""
    simples(C::ConvolutionCategory)

Return a list of simple objects in Conv(``X``).
"""
@memoize Dict function simples(C::ConvolutionCategory)
    return [ConvolutionCategoryObject(sh,C) for sh ∈ simples(C.squaredCoh)]
end

"""
    decompose(X::ConvolutionCategoryObject)

Decompose ``X`` into a direct sum of simple objects with multiplicities.
"""
function decompose(X::ConvolutionCategoryObject)
    facs = decompose(X.sheaf)
    return [(ConvolutionCategoryObject(sh,parent(X)),d) for (sh,d) ∈ facs]
end

function is_simple(X::ConvolutionCategoryObject)
    non_zero_stalks = [s for s ∈ stalks(X) if s != zero(parent(s))]
    if length(non_zero_stalks) != 1
        return false
    end
    return is_simple(non_zero_stalks[1])
end 

#-----------------------------------------------------------------
#   Hom Space
#-----------------------------------------------------------------

struct ConvCategoryHomSpace <: AbstractCategoryHomSpace
    X::ConvolutionCategoryObject
    Y::ConvolutionCategoryObject
    basis::Vector{ConvolutionCategoryMorphism}
    parent::VectorSpaces
end

"""
    Hom(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject)

Return Hom(``X,Y``) as a vector space.
"""
function Hom(X::ConvolutionCategoryObject, Y::ConvolutionCategoryObject)
    @assert parent(X) == parent(Y) "Missmatching parents"
    b = basis(Hom(X.sheaf,Y.sheaf))
    conv_b = [ConvolutionCategoryMorphism(X,Y,m) for m ∈ b]

    return ConvCategoryHomSpace(X,Y,conv_b, VectorSpaces(base_ring(X)))
end

zero(H::ConvCategoryHomSpace) = zero_morphism(H.X,H.Y)


#-----------------------------------------------------------------
#   pretty printing
#-----------------------------------------------------------------

function show(io::IO, C::ConvolutionCategory)
    print(io, """Convolution category over G-set with $(length(C.GSet)) elements.""")
end

function show(io::IO, X::ConvolutionCategoryObject)
    print(io, "CategoryObject in $(parent(X))")
end

function show(io::IO, f::ConvolutionCategoryMorphism)
    print(io, "Morphism in $(parent(domain(f)))")
end
