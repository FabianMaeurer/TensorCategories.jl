
"""
    VectorSpaces{T}(K::S) where T <: FieldElem

The category of finite dimensional vector spaces over K.
"""
struct VectorSpaces <: Category
    base_ring::Field
end


struct VSObject<: VectorSpaceObject
    basis::Vector
    parent::VectorSpaces
end

struct VSMorphism <: VectorSpaceMorphism
    m::MatElem
    domain::VectorSpaceObject
    codomain::VectorSpaceObject
end

isfusion(::VectorSpaces) = true

#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------

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
function VectorSpaceObject(Vec::VectorSpaces, n::Int)
    basis = ["v$i" for i ∈ 1:n]
    return VSObject(basis,Vec)
end


function VectorSpaceObject(K::Field,n::Int)
    Vec = VectorSpaces(K)
    return VectorSpaceObject(Vec,n)
end

function VectorSpaceObject(Vec::VectorSpaces, basis::Vector)
    return VSObject(basis, Vec)
end

function VectorSpaceObject(K::Field, basis::Vector)
    Vec = VectorSpaces(K)
    return VSObject(basis, Vec)
end
"""
    Morphism(X::VectorSpaceObject, Y::VectorSpaceObject, m::MatElem)

Return a morphism in the category of vector spaces defined by m.
"""
function Morphism(X::VectorSpaceObject, Y::VectorSpaceObject, m::MatElem)
    if parent(X) != parent(Y)
        throw(ErrorException("Missmatching parents."))
    elseif size(m) != (dim(X),dim(Y))
        throw(ErrorException("Mismatching dimensions"))
    else
        return VSMorphism(m,X,Y)
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

function Base.show(io::IO, m::VectorSpaceMorphism)
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
dim(V::VectorSpaceObject) = length(basis(V))

basis(V::VectorSpaceObject) = V.basis

simples(Vec::VectorSpaces) = [VectorSpaceObject(base_ring(Vec),1)]

decompose(V::VSObject) = [(one(parent(V)),dim(V))]

matrix(f::VectorSpaceMorphism) = f.m
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

==(V::VectorSpaces,W::VectorSpaces) = V.base_ring == W.base_ring

function ==(X::VectorSpaceObject, Y::VectorSpaceObject) where T
    basis(X) == basis(Y) && base_ring(X) == base_ring(Y)
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
    m = [matrix(f)[i] for f ∈ basis(dual(V)), i ∈ 1:dim(V)]
    Morphism(dom,cod, matrix(base_ring(V), reshape(m,dim(dom),1)))
end

function coev(V::VectorSpaceObject)
    dom = one(parent(V))
    cod = V ⊗ dual(V)
    m = [Int(i==j) for i ∈ 1:dim(V), j ∈ 1:dim(V)][:]
    Morphism(dom,cod, transpose(matrix(base_ring(V), reshape(m,dim(cod),1))))
end

spherical(V::VectorSpaceObject) = Morphism(V,dual(dual(V)), id(V).m)

#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    dsum(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, morphisms = false) where {T}

Direct sum of vector spaces together with the embedding morphisms if morphisms = true.
"""
function dsum(X::VectorSpaceObject, Y::VectorSpaceObject, morphisms::Bool = false)
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end

    if dim(X) == 0 return morphisms ? (Y,[zero_morphism(X,Y), id(Y)], [zero_morphism(Y,X), id(Y)]) : Y end
    if dim(Y) == 0 return morphisms ? (X,[id(X), zero_morphism(Y,X)], [id(X), zero_morphism(X,Y), ]) : X end

    F = base_ring(X)
    b = [(1,x) for x in basis(X)] ∪ [(2,y) for y in basis(Y)]

    V = VectorSpaceObject(parent(X),b)

    if !morphisms return V end

    ix = Morphism(X,V, matrix(F,[i == j ? 1 : 0 for i ∈ 1:dim(X), j ∈ 1:dim(V)]))
    iy = Morphism(Y,V, matrix(F,[i == j - dim(X) for i ∈ 1:dim(Y), j ∈ 1:dim(V)]))

    px = Morphism(V,X, transpose(matrix(ix)))
    py = Morphism(V,Y, transpose(matrix(iy)))

    return V,[ix,iy], [px,py]
end

product(X::VectorSpaceObject, Y::VectorSpaceObject, projections::Bool = false) = projections ? dsum(X,Y, projections)[[1,3]] : dsum(X,Y)
coproduct(X::VectorSpaceObject, Y::VectorSpaceObject, injections::Bool = false) = injections ? dsum(X,Y, injections)[[1,2]] : dsum(X,Y)

"""
    dsum(f::VectorSpaceMorphism{T},g::VectorSpaceMorphism{T}) where T

Return the direct sum of morphisms of vector spaces.
"""
function dsum(f::VectorSpaceMorphism,g::VectorSpaceMorphism)
    F = base_ring(domain(f))
    mf,nf = size(f.m)
    mg,ng = size(g.m)
    z1 = zero(MatrixSpace(F,mf,ng))
    z2 = zero(MatrixSpace(F,mg,nf))
    m = vcat(hcat(f.m,z1), hcat(z2,g.m))
    return VSMorphism(m,dsum(domain(f),domain(g)),dsum(codomain(f),codomain(g)))
end

#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

function kernel(f::VSMorphism)
    F = base_ring(domain(f))
    d,k = kernel(f.m, F, side = :left)
    k = k[1:d,:]
    K = VectorSpaceObject(parent(domain(f)), d)
    return K, Morphism(K,domain(f),k)
end

function cokernel(f::VSMorphism)
    F = base_ring(domain(f))
    d,k = kernel(f.m, F)
    k = k[:,1:d]
    K = VectorSpaceObject(parent(domain(f)), d)
    return K, Morphism(codomain(f), K, k)
end
#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

"""
    tensor_product(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T,S1,S2}

Return the tensor product of vector spaces.
"""
function tensor_product(X::VectorSpaceObject, Y::VectorSpaceObject)
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end
    b = [[(x,y) for y ∈ basis(Y), x ∈ basis(X)]...]
    return VectorSpaceObject(parent(X),b)
end

"""
    tensor_product(f::VectorSpaceMorphism, g::VectorSpaceMorphism)

Return the tensor product of vector space morphisms.
"""
function tensor_product(f::VectorSpaceMorphism, g::VectorSpaceMorphism)
    D = tensor_product(domain(f),domain(g))
    C = tensor_product(codomain(f),codomain(g))
    m = kronecker_product(f.m, g.m)
    return Morphism(D,C,m)
end
#


#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::VectorSpaceMorphism...)
    if [isisomorphic(domain(f[i]), codomain(f[i-1]))[1] for i ∈ 2:length(f)] != trues(length(f)-1)
        throw(ErrorException("Morphisms not compatible"))
    end
    return VSMorphism(*([g.m for g ∈ f]...),domain(f[1]),codomain(f[end]))
end


function ==(f::VectorSpaceMorphism, g::VectorSpaceMorphism)
    a = domain(f) == domain(g)
    b = codomain(f) == codomain(g)
    c = f.m == g.m
    return a && b && c
end

function +(f::VectorSpaceMorphism, g::VectorSpaceMorphism)
    @assert isisomorphic(domain(f),domain(g))[1] && isisomorphic(codomain(f),codomain(g))[1]
    return Morphism(domain(f),codomain(f), f.m + g.m)
end

"""
    id(X::VectorSpaceObject{T}) where T

Return the identity on the vector space ``X``.
"""
function id(X::VectorSpaceObject)
    n = dim(X)
    m = matrix(base_ring(X), [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return Morphism(X,X,m)
end

inv(f::VectorSpaceMorphism)= Morphism(codomain(f), domain(f), inv(matrix(f)))

*(λ,f::VectorSpaceMorphism)  = Morphism(domain(f),codomain(f),parent(domain(f)).base_ring(λ)*f.m)

isinvertible(f::VectorSpaceMorphism) = rank(f.m) == dim(domain(f)) == dimension(codomain(f))

function left_inverse(f::VectorSpaceMorphism)
    k = matrix(f)
    d = rank(k)
    F = base_ring(f)
    k_inv = transpose(solve_left(transpose(k), one(MatrixSpace(F,d,d))))
    return Morphism(codomain(f), domain(f), k_inv)
end

function right_inverse(f::VectorSpaceMorphism)
    k = matrix(f)
    d = rank(k)
    F = base_ring(f)
    c_inv = solve_left(k, one(MatrixSpace(F,d,d)))
    return Morphism(codomain(f),domain(f), c_inv)
end
#---------------------------------------------------------------------------
#   Associators
#---------------------------------------------------------------------------
#

"""
    associator(X::VectorSpaceObject, Y::VectorSpaceObject, Z::VectorSpaceObject)

Return the associator isomorphism a::(X⊗Y)⊗Z -> X⊗(Y⊗Z).
"""
function associator(X::VectorSpaceObject, Y::VectorSpaceObject, Z::VectorSpaceObject)
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

struct VSHomSpace <: HomSpace
    X::VectorSpaceObject
    Y::VectorSpaceObject
    basis::Vector{VectorSpaceMorphism}
    parent::VectorSpaces
end

"""
    Hom(X::VectorSpaceObject, Y::VectorSpaceObject)

Return the Hom(``X,Y```) as a vector space.
"""
function Hom(X::VectorSpaceObject, Y::VectorSpaceObject)
    n1,n2 = (dim(X),dim(Y))
    mats = [matrix(base_ring(X), [i==k && j == l ? 1 : 0 for i ∈ 1:n1, j ∈ 1:n2]) for k ∈ 1:n1, l ∈ 1:n2]
    basis = [[Morphism(X,Y,m) for m ∈ mats]...]
    return VSHomSpace(X,Y,basis,VectorSpaces(base_ring(X)))
end

basis(V::VSHomSpace) = V.basis

zero(V::VSHomSpace) = Morphism(V.X,V.Y,matrix(base_ring(V.X), [0 for i ∈ 1:dim(V.X), j ∈ 1:dim(V.Y)]))

zero_morphism(V::VectorSpaceObject,W::VectorSpaceObject) = Morphism(V,W, zero(MatrixSpace(base_ring(V), dim(V),dim(W))))

function express_in_basis(f::VectorSpaceMorphism, B::Vector{<:VectorSpaceMorphism})
    F = base_ring(f)
    B_mat = matrix(F,hcat([[x for x ∈ b.m][:] for b ∈ B]...))
    f_mat = matrix(F, 1, *(size(f.m)...), [x for x ∈ f.m][:])

    return [x for x ∈ solve_left(transpose(B_mat),f_mat)][:]
end

(F::Field)(f::VectorSpaceMorphism) = F(matrix(f)[1,1])
