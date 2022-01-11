
"""
    VectorSpaces{T}(K::S) where T <: FieldElem

The category of finite dimensional vector spaces over K.
"""
struct VectorSpaces{T<:FieldElem} <: TensorCategory{T}
    base_ring::Field
end


struct VSObject{T<:FieldElem} <: VectorSpaceObject{T}
    basis::Vector
    parent::VectorSpaces{T}
end

struct VSMorphism{T<:FieldElem} <: VectorSpaceMorphism{T}
    m::MatElem{T}
    domain::VectorSpaceObject{T}
    codomain::VectorSpaceObject{T}
end

features(::VectorSpaces) = [:semisimple,:abelian,:linear,:monoidal, :additive]

#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------

function VectorSpaces(K::T) where T <: Field
    S = elem_type(K)
    return VectorSpaces{S}(K)
end

# function (Vec::VectorSpaces{T})(V::FreeModule{T}) where T <: FieldElem
#     return VectorSpaceObject{T,FreeModule{T}}(V,Vec)
# end
#

"""
    VectorSpaceObject(Vec::VectorSpaces, n::Int64)
    VectorSpaceObject(K::Field, n::Int)
    VectorSpaceObject(Vec::VectorSpaces, basis::Vector)
    VectorSpaceObject(K::Field, basis::Vector)

The n-dimensional vector space with basis v1,..,vn (or other specified basis)
"""
function VectorSpaceObject(Vec::VectorSpaces{T}, n::Int) where T
    basis = ["v$i" for i ∈ 1:n]
    return VSObject{T}(basis,Vec)
end


function VectorSpaceObject(K::Field,n::Int)
    Vec = VectorSpaces(K)
    return VectorSpaceObject(Vec,n)
end

function VectorSpaceObject(Vec::VectorSpaces{T}, basis::Vector) where T
    return VSObject{T}(basis, Vec)
end

function VectorSpaceObject(K::F, basis::Vector{S},
        bstring::Vector{String} = String[]) where {F <: Field,S}

    Vec = VectorSpaces(K)
    return VectorSpaceObject(Vec,basis,bstring)
end

"""
    Morphism(X::VectorSpaceObject, Y::VectorSpaceObject, m::MatElem)

Return a morphism in the category of vector spaces defined by m.
"""
function Morphism(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, m::MatElem{T}) where {T}
    if parent(X) != parent(Y)
        throw(ErrorException("Missmatching parents."))
    elseif size(m) != (dim(X),dim(Y))
        throw(ErrorException("Mismatching dimensions"))
    else
        return VSMorphism{T}(m,X,Y)
    end
end

function Morphism(m::MatElem)
    l,n = size(m)
    F = base_ring(m)
    dom = VectorSpaceObject(F,l)
    codom = VectorSpaceObject(F,n)
    return Morphism(dom,codom,m)
end
#
# function VectorSpaceMorphism(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, m::U) where {T,U <: MatrixElem}
#     if parent(X) == parent(Y)
#         f = ModuleHomomorphism(X.V,Y.V,m)
#         return VectorSpaceMorphism{T}(f,X,Y)
#     else
#         throw(ErrorException("Missmatching parents."))
#     end
# end
#
# #-----------------------------------------------------------------
# #   Pretty Printing
# #-----------------------------------------------------------------
function Base.show(io::IO, C::VectorSpaces)
    print(io,"Category of finite dimensional VectorSpaces over $(C.base_ring)")
end

function Base.show(io::IO, V::VectorSpaceObject)
    print(io, "Vector space of dimension $(dim(V)) over $(base_ring(V)).")
end

function Base.show(io::IO, m::VectorSpaceMorphism{T}) where {T}
    print(io, """
Vector space morphism with
Domain:$(domain(m))
Codomain:$(codomain(m))""")
end
#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------

base_ring(V::VectorSpaceObject) = parent(V).base_ring
base_ring(Vec::VectorSpaces) = Vec.base_ring

dim(V::VectorSpaceObject) = length(V.basis)

basis(V::VectorSpaceObject) = V.basis

simples(Vec::VectorSpaces) = [VectorSpaceObject(base_ring(Vec),1)]

decompose(V::VSObject) = [(one(parent(V)),dim(V))]

one(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec),1)

zero(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec), 0)

==(V::VectorSpaces{T},W::VectorSpaces{T}) where T = V.base_ring == W.base_ring

function ==(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where T
    a = X.basis == Y.basis
    return a
end

function isisomorphic(V::VSObject, W::VSObject)
    if parent(V) != parent(W) return false, nothing end
    if dim(V) != dim(W) return false, nothing end

    return true, Morphism(V,W,one(MatrixSpace(base_ring(V),dim(V),dim(V))))
end

#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    dsum(X::VectorSpaceObject{T,S}...) where {T,S <: FreeModule}

Direct sum space of X... together with the embedding morphisms.
"""
function dsum(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, morphisms = false) where {T}
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end
    F = base_ring(X)
    b = [(1,x) for x in basis(X)] ∪ [(2,y) for y in basis(Y)]

    V = VectorSpaceObject(parent(X),b)

    if !morphisms return V end

    ix = Morphism(X,V, matrix(F,[i == j ? 1 : 0 for i ∈ 1:dim(X), j ∈ 1:dim(V)]))
    iy = Morphism(Y,V, matrix(F,[i == j - dim(X) for i ∈ 1:dim(Y), j ∈ 1:dim(V)]))

    px = Morphism(V,X, matrix(F,[i == j ? 1 : 0 for i ∈ 1:dim(V), j ∈ 1:dim(X)]))
    py = Morphism(V,Y, matrix(F,[i == j ? 1 : 0 for i ∈ 1:dim(V), j ∈ 1:dim(Y)]))

    return V,[ix,iy], [px,py]
end

product(X::VectorSpaceObject, Y::VectorSpaceObject, projections = false) = dsum(X,Y, projections)[[1,3]]
coproduct(X::VectorSpaceObject, Y::VectorSpaceObject, injections = false) = dsum(X,Y, injections)[[1,2]]

function dsum(f::VectorSpaceMorphism{T},g::VectorSpaceMorphism{T}) where T
    F = base_ring(domain(f))
    mf,nf = size(f.m)
    mg,ng = size(g.m)
    z1 = zero(MatrixSpace(F,mf,ng))
    z2 = zero(MatrixSpace(F,mg,nf))
    m = vcat(hcat(f.m,z1), hcat(z2,g.m))
    return VSMorphism{T}(m,dsum(domain(f),domain(g))[1],dsum(codomain(f),codomain(g))[1])
end


#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

function tensor_product(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T,S1,S2}
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end
    b = [[(x,y) for x ∈ basis(X), y ∈ basis(Y)]...]
    return VectorSpaceObject(parent(X),b)
end
#
⊗(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T} = tensor_product(X,Y)
#
function tensor_product(f::VectorSpaceMorphism, g::VectorSpaceMorphism) where {T}
    D = tensor_product(domain(f),domain(g))
    C = tensor_product(codomain(f),codomain(g))
    m = kronecker_product(g.m, f.m)
    return VectorSpaceMorphism(D,C,m)
end
#

#
# function compose(m1::VectorSpaceMorphism, m2::VectorSpaceMorphism)
#     return VectorSpaceMorphism(domain(m1), codomain(m2), compose(m1.m,m2.m))
# end
#
# ∘(m1::VectorSpaceMorphism,m2::VectorSpaceMorphism) = compose(m2,m1)
#


#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::VectorSpaceMorphism{T}...) where T
    if [domain(f[i]) == codomain(f[i-1]) for i ∈ 2:length(f)] != trues(length(f)-1)
        throw(ErrorException("Morphisms not compatible"))
    end
    return VSMorphism{T}(*([g.m for g ∈ f]...),domain(f[1]),codomain(f[end]))
end


function ==(f::VectorSpaceMorphism, g::VectorSpaceMorphism)
    a = domain(f) == domain(g)
    b = codomain(f) == codomain(g)
    c = f.m == g.m
    return a && b && c
end

function id(X::VectorSpaceObject{T}) where T
    n = dim(X)
    m = matrix(base_ring(X), [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return VectorSpaceMorphism(X,X,m)
end

inv(f::VectorSpaceMorphism)= VectorSpaceMorphism(domain(f), codomain(f), inv(f.m))

*(λ,f::VSMorphism)  = VectorSpaceMorphism(domain(f),codomain(f),parent(domain(f)).base_ring(λ)*f.m)

isinvertible(f::VSMorphism) = rank(f.m) == dim(domain(f)) == dimension(codomain(f))
#---------------------------------------------------------------------------
#   Associators
#---------------------------------------------------------------------------
#
function associator(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, Z::VectorSpaceObject{T}) where T
    if !(parent(X) == parent(Y) == parent(Z))
        throw(ErrorException("Mismatching parents"))
    end
    n = *(dim.([X,Y,Z])...)
    F = base_ring(X)
    m = matrix(F, [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return VectorSpaceMorphism((X⊗Y)⊗Z, X⊗(Y⊗Z), m)
end


#----------------------------------------------------------------------------
#   Hom Spaces
#----------------------------------------------------------------------------

struct VSHomSpace{T} <: HomSpace{T}
    X::VectorSpaceObject{T}
    Y::VectorSpaceObject{T}
    basis::Vector{VectorSpaceMorphism{T}}
    parent::VectorSpaces{T}
end

function Hom(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where T
    n1,n2 = (dim(X),dim(Y))
    mats = [matrix(base_ring(X), [i==k && j == l ? 1 : 0 for i ∈ 1:n1, j ∈ 1:n2]) for k ∈ 1:n1, l ∈ 1:n2]
    basis = [[VectorSpaceMorphism(X,Y,m) for m ∈ mats]...]
    return VSHomSpace{T}(X,Y,basis,VectorSpaces(base_ring(X)))
end

basis(V::VSHomSpace) = V.basis

zero(V::VSHomSpace) = VectorSpaceMorphism(V.X,V.Y,matrix(base_ring(V.X), [0 for i ∈ 1:dim(V.X), j ∈ 1:dim(V.Y)]))
