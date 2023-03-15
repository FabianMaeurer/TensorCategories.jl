
"""
    VectorSpaces{T}(K::S) where T <: FieldElem

The category of finite dimensional vector spaces over K.
"""
struct VectorSpaces <: Category
    base_ring::Field
end


struct VSCategoryObject<: VectorSpaceCategoryObject
    basis::Vector
    parent::VectorSpaces
end

struct VSCategoryMorphism <: VectorSpaceCategoryMorphism
    m::MatElem
    domain::VectorSpaceCategoryObject
    codomain::VectorSpaceCategoryObject
end

is_fusion(::VectorSpaces) = true

#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------

# function (Vec::VectorSpaces{T})(V::FreeModule{T}) where T <: FieldElem
#     return VectorSpaceCategoryObject{T,FreeModule{T}}(V,Vec)
# end
#

"""
    VectorSpaceCategoryObject(Vec::VectorSpaces, n::Int64)
    VectorSpaceCategoryObject(K::Field, n::Int)
    VectorSpaceCategoryObject(Vec::VectorSpaces, basis::Vector)
    VectorSpaceCategoryObject(K::Field, basis::Vector)

The n-dimensional vector space with basis v1,..,vn (or other specified basis)
"""
function VectorSpaceCategoryObject(Vec::VectorSpaces, n::Int)
    basis = ["v$i" for i ∈ 1:n]
    return VSCategoryObject(basis,Vec)
end


function VectorSpaceCategoryObject(K::Field,n::Int)
    Vec = VectorSpaces(K)
    return VectorSpaceCategoryObject(Vec,n)
end

function VectorSpaceCategoryObject(Vec::VectorSpaces, basis::Vector)
    return VSCategoryObject(basis, Vec)
end

function VectorSpaceCategoryObject(K::Field, basis::Vector)
    Vec = VectorSpaces(K)
    return VSCategoryObject(basis, Vec)
end
"""
    Morphism(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject, m::MatElem)

Return a morphism in the category of vector spaces defined by m.
"""
function Morphism(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject, m::MatElem)
    if parent(X) != parent(Y)
        throw(ErrorException("Missmatching parents."))
    elseif size(m) != (dim(X),dim(Y))
        throw(ErrorException("Mismatching dimensions"))
    else
        return VSCategoryMorphism(m,X,Y)
    end
end

"""
    Morphism(m::MatElem)

Vector space morphisms defined by m.
"""
function Morphism(m::MatElem)
    l,n = size(m)
    F = base_ring(m)
    dom = VectorSpaceCategoryObject(F,l)
    codom = VectorSpaceCategoryObject(F,n)
    return Morphism(dom,codom,m)
end
#
# function VectorSpaceCategoryMorphism(X::VectorSpaceCategoryObject{T}, Y::VectorSpaceCategoryObject{T}, m::U) where {T,U <: MatrixElem}
#     if parent(X) == parent(Y)
#         f = ModuleHomomorphism(X.V,Y.V,m)
#         return VectorSpaceCategoryMorphism{T}(f,X,Y)
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

function Base.show(io::IO, V::VectorSpaceCategoryObject)
    print(io, "Vector space of dimension $(int_dim(V)) over $(base_ring(V)).")
end

function Base.show(io::IO, m::VectorSpaceCategoryMorphism)
    print(io, """
Vector space morphism with
Domain:$(domain(m))
Codomain:$(codomain(m))""")
end
#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------

base_ring(V::VectorSpaceCategoryObject) = parent(V).base_ring
base_ring(Vec::VectorSpaces) = Vec.base_ring

"""
    dim(V::VectorSpaceCategoryObject) = length(V.basis)

Return the vector space dimension of ``V``.
"""
dim(V::VectorSpaceCategoryObject) = base_ring(V)(length(basis(V)))

basis(V::VectorSpaceCategoryObject) = V.basis

simples(Vec::VectorSpaces) = [VectorSpaceCategoryObject(base_ring(Vec),1)]

decompose(V::VSCategoryObject) = [(one(parent(V)),dim(V))]

matrix(f::VectorSpaceCategoryMorphism) = f.m
"""
    one(Vec::VectorSpaces) = VectorSpaceCategoryObject(base_ring(Vec),1)

Return the one-dimensional vector space.
"""
one(Vec::VectorSpaces) = VectorSpaceCategoryObject(base_ring(Vec),1)

"""
    zero(Vec::VectorSpaces) = VectorSpaceCategoryObject(base_ring(Vec), 0)

Return the zero-dimensional vector space.
"""
zero(Vec::VectorSpaces) = VectorSpaceCategoryObject(base_ring(Vec), 0)

==(V::VectorSpaces,W::VectorSpaces) = V.base_ring == W.base_ring

function ==(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject) 
    basis(X) == basis(Y) && base_ring(X) == base_ring(Y)
end

"""
    is_isomorphic(V::VSCategoryObject, W::VSCategoryObject)

Check whether ``V`` and ``W``are isomorphic. Return the isomorphisms if existent.
"""
function is_isomorphic(V::VectorSpaceCategoryObject, W::VectorSpaceCategoryObject)
    if parent(V) != parent(W) return false, nothing end
    if dim(V) != dim(W) return false, nothing end

    return true, Morphism(V,W,one(MatrixSpace(base_ring(V),int_dim(V),int_dim(V))))
end


dual(V::VectorSpaceCategoryObject) = Hom(V,one(parent(V)))

function ev(V::VectorSpaceCategoryObject)
    dom = dual(V)⊗V
    cod = one(parent(V))
    m = [matrix(f)[i] for f ∈ basis(dual(V)), i ∈ 1:int_dim(V)]
    Morphism(dom,cod, matrix(base_ring(V), reshape(m,int_dim(dom),1)))
end

function coev(V::VectorSpaceCategoryObject)
    dom = one(parent(V))
    cod = V ⊗ dual(V)
    m = [Int(i==j) for i ∈ 1:int_dim(V), j ∈ 1:int_dim(V)][:]
    Morphism(dom,cod, transpose(matrix(base_ring(V), reshape(m,int_dim(cod),1))))
end

spherical(V::VectorSpaceCategoryObject) = Morphism(V,dual(dual(V)), id(V).m)

int_dim(V::VectorSpaceCategoryObject) = length(basis(V))
#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    direct_sum(X::VectorSpaceCategoryObject{T}, Y::VectorSpaceCategoryObject{T}, morphisms = false) where {T}

Direct sum of vector spaces together with the embedding morphisms if morphisms = true.
"""
function direct_sum(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject, morphisms::Bool = false)
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end

    if dim(X) == 0 return morphisms ? (Y,[zero_morphism(X,Y), id(Y)], [zero_morphism(Y,X), id(Y)]) : Y end
    if dim(Y) == 0 return morphisms ? (X,[id(X), zero_morphism(Y,X)], [id(X), zero_morphism(X,Y), ]) : X end

    F = base_ring(X)
    b = [(1,x) for x in basis(X)] ∪ [(2,y) for y in basis(Y)]

    V = VectorSpaceCategoryObject(parent(X),b)

    if !morphisms return V end

    ix = Morphism(X,V, matrix(F,[i == j ? 1 : 0 for i ∈ 1:int_dim(X), j ∈ 1:int_dim(V)]))
    iy = Morphism(Y,V, matrix(F,[i == j - int_dim(X) for i ∈ 1:int_dim(Y), j ∈ 1:int_dim(V)]))

    px = Morphism(V,X, transpose(matrix(ix)))
    py = Morphism(V,Y, transpose(matrix(iy)))

    return V,[ix,iy], [px,py]
end

direct_sum(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject) = direct_sum(X,Y,true)

product(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject, projections::Bool = false) = projections ? direct_sum(X,Y, projections)[[1,3]] : direct_sum(X,Y)
coproduct(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject, injections::Bool = false) = injections ? direct_sum(X,Y, injections)[[1,2]] : direct_sum(X,Y)

"""
    direct_sum(f::VectorSpaceCategoryMorphism{T},g::VectorSpaceCategoryMorphism{T}) where T

Return the direct sum of morphisms of vector spaces.
"""
function direct_sum(f::VectorSpaceCategoryMorphism,g::VectorSpaceCategoryMorphism)
    F = base_ring(domain(f))
    mf,nf = size(f.m)
    mg,ng = size(g.m)
    z1 = zero(MatrixSpace(F,mf,ng))
    z2 = zero(MatrixSpace(F,mg,nf))
    m = vcat(hcat(f.m,z1), hcat(z2,g.m))
    return VSCategoryMorphism(m,direct_sum(domain(f),domain(g)),direct_sum(codomain(f),codomain(g)))
end

#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

function kernel(f::VSCategoryMorphism)
    F = base_ring(domain(f))
    d,k = kernel(f.m, F, side = :left)
    k = k[1:d,:]
    K = VectorSpaceCategoryObject(parent(domain(f)), d)
    return K, Morphism(K,domain(f),k)
end

function cokernel(f::VSCategoryMorphism)
    F = base_ring(domain(f))
    d,k = kernel(f.m, F)
    k = k[:,1:d]
    K = VectorSpaceCategoryObject(parent(domain(f)), d)
    return K, Morphism(codomain(f), K, k)
end
#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

"""
    tensor_product(X::VectorSpaceCategoryObject{T}, Y::VectorSpaceCategoryObject{T}) where {T,S1,S2}

Return the tensor product of vector spaces.
"""
function tensor_product(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject)
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end
    b = [[(x,y) for y ∈ basis(Y), x ∈ basis(X)]...]
    return VectorSpaceCategoryObject(parent(X),b)
end

"""
    tensor_product(f::VectorSpaceCategoryMorphism, g::VectorSpaceCategoryMorphism)

Return the tensor product of vector space morphisms.
"""
function tensor_product(f::VectorSpaceCategoryMorphism, g::VectorSpaceCategoryMorphism)
    D = tensor_product(domain(f),domain(g))
    C = tensor_product(codomain(f),codomain(g))
    m = kronecker_product(f.m, g.m)
    return Morphism(D,C,m)
end
#


#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::VectorSpaceCategoryMorphism...)
    if [is_isomorphic(domain(f[i]), codomain(f[i-1]))[1] for i ∈ 2:length(f)] != trues(length(f)-1)
        throw(ErrorException("Morphisms not compatible"))
    end
    return Morphism(domain(f[1]),codomain(f[end]),*([g.m for g ∈ f]...))
end


function ==(f::VectorSpaceCategoryMorphism, g::VectorSpaceCategoryMorphism)
    a = domain(f) == domain(g)
    b = codomain(f) == codomain(g)
    c = f.m == g.m
    return a && b && c
end

function +(f::VectorSpaceCategoryMorphism, g::VectorSpaceCategoryMorphism)
    @assert is_isomorphic(domain(f),domain(g))[1] && is_isomorphic(codomain(f),codomain(g))[1]
    return Morphism(domain(f),codomain(f), f.m + g.m)
end

"""
    id(X::VectorSpaceCategoryObject{T}) where T

Return the identity on the vector space ``X``.
"""
function id(X::VectorSpaceCategoryObject)
    n = int_dim(X)
    m = matrix(base_ring(X), [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return Morphism(X,X,m)
end

inv(f::VectorSpaceCategoryMorphism)= Morphism(codomain(f), domain(f), inv(matrix(f)))

*(λ,f::VectorSpaceCategoryMorphism)  = Morphism(domain(f),codomain(f),parent(domain(f)).base_ring(λ)*f.m)

isinvertible(f::VectorSpaceCategoryMorphism) = rank(f.m) == dim(domain(f)) == dimension(codomain(f))

function left_inverse(f::VectorSpaceCategoryMorphism)
    k = matrix(f)
    d = rank(k)
    F = base_ring(f)
    k_inv = transpose(solve_left(transpose(k), one(MatrixSpace(F,d,d))))
    return Morphism(codomain(f), domain(f), k_inv)
end

function right_inverse(f::VectorSpaceCategoryMorphism)
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
    associator(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject, Z::VectorSpaceCategoryObject)

Return the associator isomorphism a::(X⊗Y)⊗Z -> X⊗(Y⊗Z).
"""
function associator(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject, Z::VectorSpaceCategoryObject)
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

struct VSCategoryHomSpace <: AbstractCategoryHomSpace
    X::VectorSpaceCategoryObject
    Y::VectorSpaceCategoryObject
    basis::Vector{VectorSpaceCategoryMorphism}
    parent::VectorSpaces
end

"""
    Hom(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject)

Return the Hom(``X,Y```) as a vector space.
"""
function Hom(X::VectorSpaceCategoryObject, Y::VectorSpaceCategoryObject)
    n1,n2 = (length(basis(X)), length(basis(Y)))
    mats = [matrix(base_ring(X), [i==k && j == l ? 1 : 0 for i ∈ 1:n1, j ∈ 1:n2]) for k ∈ 1:n1, l ∈ 1:n2]
    base = [[Morphism(X,Y,m) for m ∈ mats]...]
    return VSCategoryHomSpace(X,Y,base,VectorSpaces(base_ring(X)))
end

basis(V::VSCategoryHomSpace) = V.basis

zero(V::VSCategoryHomSpace) = Morphism(V.X,V.Y,matrix(base_ring(V.X), [0 for i ∈ 1:dim(V.X), j ∈ 1:dim(V.Y)]))

zero_morphism(V::VectorSpaceCategoryObject,W::VectorSpaceCategoryObject) = Morphism(V,W, zero(MatrixSpace(base_ring(V), int_dim(V), int_dim(W))))

function express_in_basis(f::VectorSpaceCategoryMorphism, B::Vector{<:VectorSpaceCategoryMorphism})
    F = base_ring(f)
    B_mat = matrix(F,hcat([[x for x ∈ b.m][:] for b ∈ B]...))
    f_mat = matrix(F, 1, *(size(f.m)...), [x for x ∈ f.m][:])

    return [x for x ∈ solve_left(transpose(B_mat),f_mat)][:]
end

(F::Field)(f::VectorSpaceCategoryMorphism) = F(matrix(f)[1,1])
