
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

function VectorSpaceObject(K::Field, basis::Vector)
    Vec = VectorSpaces(K)
    return VSObject{elem_type(K)}(basis, Vec)
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

"""
    Morphism(m::MatElem)

Vector space morphisms defined by m.
"""
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

"""
    dim(V::VectorSpaceObject) = length(V.basis)

Return the vector space dimension of ``V``.
"""
dim(V::VectorSpaceObject) = length(V.basis)

basis(V::VectorSpaceObject) = V.basis

simples(Vec::VectorSpaces) = [VectorSpaceObject(base_ring(Vec),1)]

decompose(V::VSObject) = [(one(parent(V)),dim(V))]

matrix(f::VSMorphism) = f.m
"""
    one(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec),1)

Return the one-dimensional vector space.
"""
one(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec),1)

"""
    zero(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec), 0)

Return the zero-dimensional vector space.
"""
zero(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec), 0)

==(V::VectorSpaces{T},W::VectorSpaces{T}) where T = V.base_ring == W.base_ring

function ==(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where T
    a = X.basis == Y.basis
    return a
end

"""
    isisomorphic(V::VSObject, W::VSObject)

Check whether ``V`` and ``W``are isomorphic. Return the isomorphisms if existent.
"""
function isisomorphic(V::VectorSpaceObject, W::VectorSpaceObject)
    if parent(V) != parent(W) return false, nothing end
    if dim(V) != dim(W) return false, nothing end

    return true, Morphism(V,W,one(MatrixSpace(base_ring(V),dim(V),dim(V))))
end


dual(V::VectorSpaceObject) = Hom(V,one(parent(V)))

function ev(V::VectorSpaceObject)
    dom = dual(V)⊗V
    cod = one(parent(V))
    m = [matrix(f)[i] for f∈ basis(dual(V)), i ∈ 1:dim(V)]
    Morphism(dom,cod, matrix(base_ring(V), reshape(m,dim(dom),1)))
end

function coev(V::VectorSpaceObject)
    dom = one(parent(V))
    cod = V ⊗ dual(V)
    m = [Int(i==j) for i ∈ 1:dim(V), j ∈ 1:dim(V)][:]
    Morphism(dom,cod, transpose(matrix(base_ring(V), reshape(m,dim(cod),1))))
end

#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    dsum(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, morphisms = false) where {T}

Direct sum of vector spaces together with the embedding morphisms if morphisms = true.
"""
function dsum(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, morphisms::Bool = false) where {T}
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end

    if dim(X) == 0 return Y end
    if dim(Y) == 0 return X end
    
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

product(X::VectorSpaceObject, Y::VectorSpaceObject, projections::Bool = false) = projections ? dsum(X,Y, projections)[[1,3]] : dsum(X,Y)
coproduct(X::VectorSpaceObject, Y::VectorSpaceObject, injections::Bool = false) = injections ? dsum(X,Y, injections)[[1,2]] : dsum(X,Y)

"""
    dsum(f::VectorSpaceMorphism{T},g::VectorSpaceMorphism{T}) where T

Return the direct sum of morphisms of vector spaces.
"""
function dsum(f::VectorSpaceMorphism{T},g::VectorSpaceMorphism{T}) where T
    F = base_ring(domain(f))
    mf,nf = size(f.m)
    mg,ng = size(g.m)
    z1 = zero(MatrixSpace(F,mf,ng))
    z2 = zero(MatrixSpace(F,mg,nf))
    m = vcat(hcat(f.m,z1), hcat(z2,g.m))
    return VSMorphism{T}(m,dsum(domain(f),domain(g)),dsum(codomain(f),codomain(g)))
end


#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

"""
    tensor_product(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T,S1,S2}

Return the tensor product of vector spaces.
"""
function tensor_product(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T,S1,S2}
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end
    b = [[(x,y) for x ∈ basis(X), y ∈ basis(Y)]...]
    return VectorSpaceObject(parent(X),b)
end

"""
    tensor_product(f::VectorSpaceMorphism, g::VectorSpaceMorphism) where {T}

Return the tensor product of vector space morphisms.
"""
function tensor_product(f::VectorSpaceMorphism, g::VectorSpaceMorphism) where {T}
    D = tensor_product(domain(f),domain(g))
    C = tensor_product(codomain(f),codomain(g))
    m = kronecker_product(g.m, f.m)
    return Morphism(D,C,m)
end
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

function +(f::VSMorphism, g::VSMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g)
    return Morphism(domain(f),codomain(f), f.m + g.m)
end

"""
    id(X::VectorSpaceObject{T}) where T

Return the identity on the vector space ``X``.
"""
function id(X::VectorSpaceObject{T}) where T
    n = dim(X)
    m = matrix(base_ring(X), [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return Morphism(X,X,m)
end

inv(f::VectorSpaceMorphism)= Morphism(domain(f), codomain(f), inv(f.m))

*(λ,f::VSMorphism)  = Morphism(domain(f),codomain(f),parent(domain(f)).base_ring(λ)*f.m)

isinvertible(f::VSMorphism) = rank(f.m) == dim(domain(f)) == dimension(codomain(f))
#---------------------------------------------------------------------------
#   Associators
#---------------------------------------------------------------------------
#

"""
    associator(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, Z::VectorSpaceObject{T}) where T

Return the associator isomorphism a::(X⊗Y)⊗Z -> X⊗(Y⊗Z).
"""
function associator(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, Z::VectorSpaceObject{T}) where T
    if !(parent(X) == parent(Y) == parent(Z))
        throw(ErrorException("Mismatching parents"))
    end
    n = *(dim.([X,Y,Z])...)
    F = base_ring(X)
    m = matrix(F, [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return Morphism((X⊗Y)⊗Z, X⊗(Y⊗Z), m)
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

"""
    Hom(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where T

Return the Hom(``X,Y```) as a vector space.
"""
function Hom(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where T
    n1,n2 = (dim(X),dim(Y))
    mats = [matrix(base_ring(X), [i==k && j == l ? 1 : 0 for i ∈ 1:n1, j ∈ 1:n2]) for k ∈ 1:n1, l ∈ 1:n2]
    basis = [[Morphism(X,Y,m) for m ∈ mats]...]
    return VSHomSpace{T}(X,Y,basis,VectorSpaces(base_ring(X)))
end

basis(V::VSHomSpace) = V.basis

zero(V::VSHomSpace) = Morphism(V.X,V.Y,matrix(base_ring(V.X), [0 for i ∈ 1:dim(V.X), j ∈ 1:dim(V.Y)]))

zero_morphism(V::VSObject,W::VSObject) = Morphism(V,W,matrix(base_ring(V), [0 for i ∈ 1:dim(V), j ∈ 1:dim(W)]))
