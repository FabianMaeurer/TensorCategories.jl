
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

is_fusion(::VectorSpaces) = true

#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------

# function (Vec::VectorSpaces{T})(V::FreeModule{T}) where T <: FieldElem
#     return VectorSpaceObject{T,FreeModule{T}}(V,Vec)
# end
#

vector_spaces(K::Field = QQ) = VectorSpaces(K)

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
    morphism(X::VectorSpaceObject, Y::VectorSpaceObject, m::MatElem)

Return a morphism in the category of vector spaces defined by m.
"""
function morphism(X::VectorSpaceObject, Y::VectorSpaceObject, m::MatElem)
    if parent(X) != parent(Y)
        throw(ErrorException("Missmatching parents."))
    elseif size(m) != (int_dim(X),int_dim(Y))
        throw(ErrorException("Mismatching dimensions"))
    else
        return VSMorphism(m,X,Y)
    end
end

"""
    morphism(m::MatElem)

Vector space morphisms defined by m.
"""
function morphism(m::MatElem)
    l,n = size(m)
    F = base_ring(m)
    dom = VectorSpaceObject(F,l)
    codom = VectorSpaceObject(F,n)
    return morphism(dom,codom,m)
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
    print(io, "Vector space of dimension $(int_dim(V)) over $(base_ring(V)).")
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
dim(V::VectorSpaceObject) = base_ring(V)(length(basis(V)))

basis(V::VectorSpaceObject) = V.basis

simples(Vec::VectorSpaces) = [VectorSpaceObject(base_ring(Vec),1)]

decompose(V::VSObject) = [(one(parent(V)),int_dim(V))]

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

function ==(X::VectorSpaceObject, Y::VectorSpaceObject) 
    int_dim(X) == int_dim(Y) && base_ring(X) == base_ring(Y)
end

"""
    is_isomorphic(V::VSObject, W::VSObject)

Check whether ``V`` and ``W``are isomorphic. Return the isomorphisms if existent.
"""
function is_isomorphic(V::VectorSpaceObject, W::VectorSpaceObject)
    if parent(V) != parent(W) return false, nothing end
    if int_dim(V) != int_dim(W) return false, nothing end

    return true, morphism(V,W,one(matrix_space(base_ring(V),int_dim(V),int_dim(V))))
end


dual(V::VectorSpaceObject) = Hom(V,one(parent(V)))

function ev(V::VectorSpaceObject)
    dom = dual(V)⊗V
    cod = one(parent(V))
    m = [matrix(f)[i] for f ∈ basis(dual(V)), i ∈ 1:int_dim(V)]
    morphism(dom,cod, matrix(base_ring(V), reshape(m,int_dim(dom),1)))
end

function coev(V::VectorSpaceObject)
    dom = one(parent(V))
    cod = V ⊗ dual(V)
    m = [Int(i==j) for i ∈ 1:int_dim(V), j ∈ 1:int_dim(V)][:]
    morphism(dom,cod, transpose(matrix(base_ring(V), reshape(m,int_dim(cod),1))))
end

spherical(V::VectorSpaceObject) = morphism(V,dual(dual(V)), id(V).m)

int_dim(V::VectorSpaceObject) = length(basis(V))
#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    direct_sum(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T}

Direct sum of vector spaces together with the embedding morphisms.
"""
function direct_sum(X::VectorSpaceObject, Y::VectorSpaceObject,)
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end

    if int_dim(X) == 0 return (Y,[zero_morphism(X,Y), id(Y)], [zero_morphism(Y,X), id(Y)])  end
    if int_dim(Y) == 0 return (X,[id(X), zero_morphism(Y,X)], [id(X), zero_morphism(X,Y), ]) end

    F = base_ring(X)
    b = [(1,x) for x in basis(X)] ∪ [(2,y) for y in basis(Y)]

    V = VectorSpaceObject(parent(X),b)

    ix = morphism(X,V, matrix(F,[i == j ? 1 : 0 for i ∈ 1:int_dim(X), j ∈ 1:int_dim(V)]))
    iy = morphism(Y,V, matrix(F,[i == j - int_dim(X) for i ∈ 1:int_dim(Y), j ∈ 1:int_dim(V)]))

    px = morphism(V,X, transpose(matrix(ix)))
    py = morphism(V,Y, transpose(matrix(iy)))

    return V,[ix,iy], [px,py]
end



"""
    direct_sum(f::VectorSpaceMorphism{T},g::VectorSpaceMorphism{T}) where T

Return the direct sum of morphisms of vector spaces.
"""
function direct_sum(f::VectorSpaceMorphism,g::VectorSpaceMorphism)
    F = base_ring(domain(f))
    mf,nf = size(f.m)
    mg,ng = size(g.m)
    z1 = zero(matrix_space(F,mf,ng))
    z2 = zero(matrix_space(F,mg,nf))
    m = vcat(hcat(f.m,z1), hcat(z2,g.m))
    return VSMorphism(m,direct_sum(domain(f),domain(g))[1],direct_sum(codomain(f),codomain(g))[1])
end

#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

function kernel(f::VSMorphism)
    F = base_ring(domain(f))
    k = kernel(f.m)
    K = VectorSpaceObject(parent(domain(f)), number_of_rows(k))
    return K, morphism(K,domain(f), k)
end

function cokernel(f::VSMorphism)
    F = base_ring(domain(f))
    k = kernel(f.m, side = :right)

    K = VectorSpaceObject(parent(domain(f)), number_of_columns(k))
    return K, morphism(codomain(f), K, k)
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
    return morphism(D,C,m)
end
#

function braiding(X::VectorSpaceObject, Y::VectorSpaceObject)
    F = base_ring(X)
    n,m = int_dim(X),int_dim(Y)
    map = zero(matrix_space(F,n*m,n*m))
    for i ∈ 1:n, j ∈ 1:m
        v1 = matrix(F,transpose([k == i ? 1 : 0 for k ∈ 1:n]))
        v2 = matrix(F,transpose([k == j ? 1 : 0 for k ∈ 1:m]))
        map[(j-1)*n + i, :] = kronecker_product(v1,v2)
    end
    return morphism(X⊗Y, Y⊗X, transpose(map))
end

#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::VectorSpaceMorphism...)
    if length(f) == 1
        return f[1]
    end
    @assert all([is_isomorphic(domain(f[i]), codomain(f[i-1]))[1] for i ∈ 2:length(f)])  "Morphisms not compatible"

    return morphism(domain(f[1]),codomain(f[end]),*([g.m for g ∈ f]...))
end


function ==(f::VectorSpaceMorphism, g::VectorSpaceMorphism)
    a = domain(f) == domain(g)
    b = codomain(f) == codomain(g)
    c = f.m == g.m
    return a && b && c
end

function +(f::VectorSpaceMorphism, g::VectorSpaceMorphism)
    @assert is_isomorphic(domain(f),domain(g))[1] && is_isomorphic(codomain(f),codomain(g))[1]
    return morphism(domain(f),codomain(f), f.m + g.m)
end

"""
    id(X::VectorSpaceObject{T}) where T

Return the identity on the vector space ``X``.
"""
function id(X::VectorSpaceObject)
    n = int_dim(X)
    m = matrix(base_ring(X), [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return morphism(X,X,m)
end

inv(f::VectorSpaceMorphism)= morphism(codomain(f), domain(f), inv(matrix(f)))

*(λ,f::VectorSpaceMorphism)  = morphism(domain(f),codomain(f),parent(domain(f)).base_ring(λ)*f.m)

is_invertible(f::VectorSpaceMorphism) = rank(f.m) == int_dim(domain(f)) == int_dim(codomain(f))

function left_inverse(f::VectorSpaceMorphism)
    k = matrix(f)
    d = rank(k)
    F = base_ring(f)
    
    k_inv = transpose(solve(transpose(k), one(matrix_space(F,d,d))))
    return morphism(codomain(f), domain(f), k_inv)
end

function right_inverse(f::VectorSpaceMorphism)
    k = matrix(f)
    d = rank(k)
    F = base_ring(f)
    c_inv = solve(k, one(matrix_space(F,d,d)))
    return morphism(codomain(f),domain(f), c_inv)
end

function extension_of_scalars(C::VectorSpaces, K::Ring)
    VectorSpaces(K)
end

function extension_of_scalars(V::VectorSpaceObject, K::Ring; parent::VectorSpaces = parent(V) ⊗ K)
    VSObject(basis(V), parent)
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
    n = *(int_dim.([X,Y,Z])...)
    F = base_ring(X)
    m = matrix(F, [i == j ? 1 : 0 for i ∈ 1:n, j ∈ 1:n])
    return morphism((X⊗Y)⊗Z, X⊗(Y⊗Z), m)
end


#----------------------------------------------------------------------------
#   Hom Spaces
#----------------------------------------------------------------------------

struct VSHomSpace <: AbstractHomSpace
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
    n1,n2 = (length(basis(X)), length(basis(Y)))
    mats = [matrix(base_ring(X), [i==k && j == l ? 1 : 0 for i ∈ 1:n1, j ∈ 1:n2]) for k ∈ 1:n1, l ∈ 1:n2]
    base = [[morphism(X,Y,m) for m ∈ mats]...]
    return VSHomSpace(X,Y,base,VectorSpaces(base_ring(X)))
end

basis(V::VSHomSpace) = V.basis

zero(V::VSHomSpace) = morphism(V.X,V.Y,matrix(base_ring(V.X), [0 for i ∈ 1:int_dim(V.X), j ∈ 1:int_dim(V.Y)]))

zero_morphism(V::VectorSpaceObject,W::VectorSpaceObject) = morphism(V,W, zero(matrix_space(base_ring(V), int_dim(V), int_dim(W))))

function express_in_basis(f::VectorSpaceMorphism, B::Vector{<:VectorSpaceMorphism})
    F = base_ring(f)

    if typeof(F) <: Union{AcbField, ComplexField, ArbField}
        return express_in_basis_numeric(f,B)
    end

    B_mat = matrix(F,hcat([[x for x ∈ b.m][:] for b ∈ B]...))
    f_mat = matrix(F, 1, *(size(f.m)...), [x for x ∈ f.m][:])
    
    return [x for x ∈ solve(transpose(B_mat),f_mat, side = :left)][:]
end

function express_in_basis_numeric(f::VectorSpaceMorphism, B::Vector{<:VectorSpaceMorphism})
    F = base_ring(f)
    #m = Int(floor(precision(F)/2))
    
    B_mat = matrix(F,hcat([[x for x ∈ b.m][:] for b ∈ B]...))
    f_mat = matrix(F, *(size(f.m)...), 1, [x for x ∈ f.m][:])

    m = minimum([Int(floor(minimum([a for a in Oscar.accuracy_bits.([B_mat f_mat]) if a > 0], init = precision(F)))), Int(floor(precision(F)))])

    coeffs = complex_lindep(collect([f_mat B_mat]), m)

    n = length(B)
    lead = F(coeffs[1])
    overlaps(lead,F(0)) && return zeros(F,n)

    return  -1 .* coeffs[2:end]
end

(F::Field)(f::VectorSpaceMorphism) = F(matrix(f)[1,1])
(F::QQBarField)(f::VectorSpaceMorphism) = F(matrix(f)[1,1])
(F::CalciumField)(f::VectorSpaceMorphism) = F(matrix(f)[1,1])