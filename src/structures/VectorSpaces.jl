
"""
    VectorSpaces{T}(K::S) where T <: FieldElem

The category of finite dimensional vector spaces over K.
"""
struct VectorSpaces{T<:FieldElem} <: TensorCategory{T}
    base_ring::Field
end

"""
    VectorSpaceObject{T}

An object in the category of finite dimensional vector spaces.
"""
struct VectorSpaceObject{T<:FieldElem} <: Object
    basis::Vector{S} where S
    basis_strings::Vector{String}
    parent::VectorSpaces{T}
end

struct VectorSpaceMorphism{T<:FieldElem} <: Morphism
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


function VectorSpaceObject(Vec::VectorSpaces{T}, n::Int, basis::Vector{String} = String[]) where T
    if basis == String[]
        basis = ["v$i" for i ∈ 1:n]
    elseif length(basis) != n
        throw(ErrorException("Mismatching dimensions"))
    end
    return VectorSpaceObject{T}(basis,basis,Vec)
end


function VectorSpaceObject(K::F,n::Int, basis::Vector{String} = String[]) where {F<:Field}
    Vec = VectorSpaces(K)
    return VectorSpaceObject(Vec,n,basis)
end

function VectorSpaceObject(Vec::VectorSpaces{T}, basis::Vector{S},
        bstring::Vector{String} = String[]) where {T,S}

        if bstring == String[]
            bstring = ["v$i" for i ∈ 1:length(basis)]
        elseif length(basis) != length(bstring)
            throw(ErrorException("Mismatching dimensions"))
        end
        return VectorSpaceObject{T}(basis, bstring, Vec)
end

function VectorSpaceObject(K::F, basis::Vector{S},
        bstring::Vector{String} = String[]) where {F <: Field,S}

    Vec = VectorSpaces(K)
    return VectorSpaceObject(Vec,basis,bstring)
end

function VectorSpaceMorphism(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}, m::MatElem{T}) where {T}
    if parent(X) != parent(Y)
        throw(ErrorException("Missmatching parents."))
    elseif size(m) != (dim(X),dim(Y))
        throw(ErrorException("Mismatching dimensions"))
    else
        return VectorSpaceMorphism{T}(m,X,Y)
    end
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

parent(X::VectorSpaceObject) = X.parent

base_ring(V::VectorSpaceObject) = V.parent.base_ring
base_ring(Vec::VectorSpaces) = Vec.base_ring

dim(V::VectorSpaceObject) = length(V.basis)

basis(V::VectorSpaceObject) = V.basis

simples(Vec::VectorSpaces) = [VectorSpaceObject(base_ring(Vec),1)]

one(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec),1)

zero(Vec::VectorSpaces) = VectorSpaceObject(base_ring(Vec), 0)

==(V::VectorSpaces{T},W::VectorSpaces{T}) where T = V.base_ring == W.base_ring

function ==(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where T
    a = X.basis == Y.basis
    return a
end

#-----------------------------------------------------------------
#   Functionality: Direct Sum
#-----------------------------------------------------------------

"""
    dsum(X::VectorSpaceObject{T,S}...) where {T,S <: FreeModule}

Direct sum space of X... together with the embedding morphisms.
"""
function dsum(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T}
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end
    F = base_ring(X)
    b = [(1,x) for x in basis(X)] ∪ [(2,y) for y in basis(Y)]
    bstring = ["($x,0)" for x in X.basis_strings] ∪ ["(0,$y)" for y in Y.basis_strings]
    V = VectorSpaceObject(parent(X),b,bstring)
    ix = VectorSpaceMorphism(X,V, matrix(F,[i == j ? 1 : 0 for i ∈ 1:dim(X), j ∈ 1:dim(V)]))
    iy = VectorSpaceMorphism(Y,V, matrix(F,[i == j - dim(X) for i ∈ 1:dim(Y), j ∈ 1:dim(V)]))
    return V,[ix,iy]
end

function dsum(f::VectorSpaceMorphism{T},g::VectorSpaceMorphism{T}) where T
    F = base_ring(domain(f))
    mf,nf = size(f.m)
    mg,ng = size(g.m)
    z1 = zero(MatrixSpace(F,mf,ng))
    z2 = zero(MatrixSpace(F,mg,nf))
    m = vcat(hcat(f.m,z1), hcat(z2,g.m))
    return VectorSpaceMorphism{T}(m,dsum(domain(f),domain(g))[1],dsum(codomain(f),codomain(g))[1])
end



#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

function tensor_product(X::VectorSpaceObject{T}, Y::VectorSpaceObject{T}) where {T,S1,S2}
    if parent(X) != parent(Y)
        throw(ErrorException("Mismatching parents."))
    end
    b = [[(x,y) for x ∈ basis(X), y ∈ basis(Y)]...]
    bstring = [["$x⊗$y" for x ∈ X.basis_strings, y ∈ Y.basis_strings]...]
    return VectorSpaceObject(parent(X),b,bstring)
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
⊗(f::VectorSpaceMorphism, g::VectorSpaceMorphism) where {T} = tensor_product(f,g)
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
    return VectorSpaceMorphism{T}(*([g.m for g ∈ f]...),domain(f[1]),codomain(f[end]))
end

∘(f::VectorSpaceMorphism{T}, g::VectorSpaceMorphism) where T = compose(g,f)

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

inv(f::VectorSpaceMorphism{T}) where T = VectorSpaceMorphism(domain(f), codomain(f), inv(f.m))

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
