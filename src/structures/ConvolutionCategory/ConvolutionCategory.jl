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
    sheaf::CohSheaf
    parent::ConvolutionCategory
end

struct ConvolutionMorphism <: Morphism
    domain::ConvolutionObject
    codomain::ConvolutionObject
    m::CohSheafMorphism
end

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

Morphism(D::ConvolutionObject, C::ConvolutionObject, m:: CohSheafMorphism) = ConvolutionMorphism(D,C,m)
#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------

"""
    issemisimple(C::ConvolutionCategory)

Check whether ``C`` semisimple.
"""
issemisimple(C::ConvolutionCategory) = gcd(order(C.group), characteristic(base_ring(C))) == 1

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

==(X::ConvolutionObject,Y::ConvolutionObject) = X.sheaf == Y.sheaf

==(f::ConvolutionMorphism, g::ConvolutionMorphism) = f.m == f.m

"""
    isisomorphic(X::ConvolutionObject, Y::ConvolutionObject)

Check whether ``X``and ``Y``are isomorphic. Return the isomorphism if true.
"""
function isisomorphic(X::ConvolutionObject, Y::ConvolutionObject)
    b, iso = isisomorphic(X.sheaf, Y.sheaf)
    if !b return false, nothing end
    return true, ConvolutionMorphism(X,Y,iso)
end

id(X::ConvolutionObject) = ConvolutionMorphism(X,X,id(X.sheaf))

function associator(X::ConvolutionObject, Y::ConvolutionObject, Z::ConvolutionObject)
    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)
    _,f = isisomorphic(dom,cod)
    return f
end
#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    dsum(X::ConvolutionObject, Y::ConvolutionObject, morphisms::Bool = false)

documentation
"""
function dsum(X::ConvolutionObject, Y::ConvolutionObject, morphisms::Bool = false)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Z,ix,px = dsum(X.sheaf,Y.sheaf,true)
    Z = ConvolutionObject(Z,parent(X))

    if !morphisms return Z end

    ix = [ConvolutionMorphism(x,Z,i) for (x,i) ∈ zip([X,Y],ix)]
    px = [ConvolutionMorphism(Z,x,p) for (x,p) ∈ zip([X,Y],px)]
    return Z, ix, px
end

"""
    dsum(f::ConvolutionMorphism, g::ConvolutionMorphism)

Return the direct sum of morphisms of coherent sheaves (with convolution product).
"""
function dsum(f::ConvolutionMorphism, g::ConvolutionMorphism)
    dom = dsum(domain(f), domain(g))
    codom = dsum(codomain(f), codomain(g))
    m = dsum(f.m,g.m)
    return ConvolutionMorphism(dom,codom,m)
end

product(X::ConvolutionObject,Y::ConvolutionObject,projections = false) = projections ? dsum(X,Y,projections)[[1,3]] : dsum(X,Y)
coproduct(X::ConvolutionObject,Y::ConvolutionObject,projections = false) = projections ? dsum(X,Y,projections)[[1,2]] : dsum(X,Y)

"""
    zero(C::ConvolutionCategory)

Return the zero object in Conv(``X``).
"""
zero(C::ConvolutionCategory) = ConvolutionObject(zero(C.squaredCoh),C)

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
    return ConvolutionObject(CohSheaf(C.squaredCoh, stlks), C)
end

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
#-----------------------------------------------------------------
#   Simple Objects
#-----------------------------------------------------------------

"""
    simples(C::ConvolutionCategory)

Return a list of simple objects in Conv(``X``).
"""
function simples(C::ConvolutionCategory)
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

#-----------------------------------------------------------------
#   Hom Space
#-----------------------------------------------------------------

struct ConvHomSpace <: HomSpace
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
