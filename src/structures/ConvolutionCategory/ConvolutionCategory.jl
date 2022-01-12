struct ConvolutionCategory{T,G} <: MultiTensorCategory{T}
    group::G
    base_ring::Field
    GSet::GSet
    squaredGSet::GSet
    cubedGSet::GSet
    squaredCoh::CohSheaves{T,G}
    cubedCoh::CohSheaves{T,G}
    projectors::Vector
end

struct ConvolutionObject{T,G} <: Object
    sheaf::CohSheaf{T,G}
    parent::ConvolutionCategory{T,G}
end

struct ConvolutionMorphism{T,G} <: Morphism
    domain::ConvolutionObject{T,G}
    codomain::ConvolutionObject{T,G}
    m::CohSheafMorphism{T,G}
end

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
    return ConvolutionCategory{elem_type(K),typeof(G)}(G,K,X,sqX,cuX,sqCoh,cuCoh, [P12, P13, P23])
end

function ConvolutionCategory(X, K::Field)
    G = symmetric_group(1)
    return ConvolutionCategory(gset(G,X), K)
end
#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------
issemisimple(C::ConvolutionCategory) = gcd(order(C.group), characteristic(base_ring(C))) == 1

orbit_stabilizers(C::ConvolutionCategory) = C.squaredCoh.orbit_stabilizers
orbit_index(X::ConvolutionObject, y) = orbit_index(X.sheaf, y)
orbit_index(X::ConvolutionCategory, y) = orbit_index(X.squaredCoh, y)
stalks(X::ConvolutionObject) = stalks(X.sheaf)

==(X::ConvolutionObject,Y::ConvolutionObject) = X.sheaf == Y.sheaf

function isisomorphic(X::ConvolutionObject{T,G}, Y::ConvolutionObject{T,G}) where {T,G}
    b, iso = isisomorphic(X.sheaf, Y.sheaf)
    if !b return false, nothing end
    return true, ConvolutionMorphism{T,G}(X,Y,iso)
end

#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

function dsum(X::ConvolutionObject{T,G}, Y::ConvolutionObject{T,G}, morphisms::Bool = false) where {T,G}
    @assert parent(X) == parent(Y) "Mismatching parents"
    Z,ix,px = dsum(X.sheaf,Y.sheaf,true)
    Z = ConvolutionObject{T,G}(Z,parent(X))

    if !morphisms return Z end

    ix = [ConvolutionMorphism{T,G}(x,Z,i) for (x,i) ∈ zip([X,Y],ix)]
    px = [ConvolutionMorphism{T,G}(Z,x,p) for (x,p) ∈ zip([X,Y],px)]
    return Z, ix, px
end

function dsum(f::ConvolutionMorphism{T,G}, g::ConvolutionMorphism{T,G}) where {T,G}
    dom = dsum(domain(f), domain(g))[1]
    codom = dsum(codomain(f), codomain(g))[1]
    m = dsum(f.m,g.m)
    return ConvolutionMorphism{T,G}(dom,codom,m)
end

product(X::ConvolutionObject,Y::ConvolutionObject,projections = false) = projections ? dsum(X,Y,projections)[[1,3]] : dsum(X,Y)
coproduct(X::ConvolutionObject,Y::ConvolutionObject,projections = false) = projections ? dsum(X,Y,projections)[[1,2]] : dsum(X,Y)


zero(C::ConvolutionCategory{T,G}) where {T,G} = ConvolutionObject(zero(C.squaredCoh),C)

#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

function tensor_product(X::ConvolutionObject{T,G}, Y::ConvolutionObject{T,G}) where {T,G}
    @assert parent(X) == parent(Y) "Mismatching parents"
    p12,p13,p23 = parent(X).projectors

    return ConvolutionObject{T,G}(p13(p12(X.sheaf)⊗p23(Y.sheaf)),parent(X))
end

function tensor_product(f::ConvolutionMorphism{T,G}, g::ConvolutionMorphism{T,G}) where {T,G}
    dom = domain(f)⊗domain(g)
    codom = codomain(f)⊗codomain(g)

    p12,p13,p23 = parent(domain(f)).projectors
    return ConvolutionMorphism{T,G}(dom,codom, p13(p12(f.m)⊗p23(g.m)))
end

function one(C::ConvolutionCategory{T,G}) where {T,G}
    F = base_ring(C)

    stlks = [zero(RepresentationCategory(H,F)) for H ∈ orbit_stabilizers(C)]
    diag = [(x,x) for x ∈ C.GSet.seeds]

    for i ∈ [orbit_index(C,d) for d ∈ diag]
        stlks[i] = one(RepresentationCategory(orbit_stabilizers(C)[i], F))
    end
    return ConvolutionObject{T,G}(CohSheaf{T,G}(C.squaredCoh, stlks), C)
end

#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::ConvolutionMorphism{T,G},g::ConvolutionMorphism{T,G}) where {T,G}
    return ConvolutionMorphism{T,G}(domain(f),codomain(g),compose(f.m,g.m))
end

function zero_morphism(X::ConvolutionObject{T,G}, Y::ConvolutionObject{T,G}) where {T,G}
    return ConvolutionMorphism{T,G}(X,Y,zero_morphism(X.sheaf,Y.sheaf))
end

#-----------------------------------------------------------------
#   Simple Objects
#-----------------------------------------------------------------

function simples(C::ConvolutionCategory{T,G}) where {T,G}
    return [ConvolutionObject{T,G}(sh,C) for sh ∈ simples(C.squaredCoh)]
end

function decompose(X::ConvolutionObject)
    facs = decompose(X.sheaf)
    return [(ConvolutionObject(sh,parent(X)),d) for (sh,d) ∈ facs]
end

#-----------------------------------------------------------------
#   Hom Space
#-----------------------------------------------------------------

struct ConvHomSpace{T,G} <: HomSpace{T}
    X::ConvolutionObject{T,G}
    Y::ConvolutionObject{T,G}
    basis::Vector{ConvolutionMorphism{T,G}}
    parent::VectorSpaces{T}
end

function Hom(X::ConvolutionObject{T,G}, Y::ConvolutionObject{T,G}) where {T,G}
    @assert parent(X) == parent(Y) "Missmatching parents"
    b = basis(Hom(X.sheaf,Y.sheaf))
    conv_b = [ConvolutionMorphism{T,G}(X,Y,m) for m ∈ b]

    return ConvHomSpace{T,G}(X,Y,conv_b, VectorSpaces(base_ring(X)))
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
