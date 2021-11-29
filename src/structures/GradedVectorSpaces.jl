
"""
    GradedVectorSpaces{T}(K::S) where T <: FieldElem

The category of finite dimensional vector spaces over K.
"""
struct GradedVectorSpaces{T,G} <: TensorCategory{T}
    base_ring::Field
    base_group::GAPGroup
end

struct GradedVectorSpaceObject{T,G} <: Object
    V::Dict{G,VectorSpaceObject{T}}
    parent::GradedVectorSpaces{T,G}
end

struct GradedVectorSpaceMorphism{T,G} <: Morphism
    m::Dict{G,VectorSpaceMorphism{T}}
    domain::GradedVectorSpaceObject{T,G}
    codomain::GradedVectorSpaceObject{T,G}
end
#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------

function GradedVectorSpaces(K::T,G::S) where {T <: Field, S <: GAPGroup}
    F = elem_type(K)
    G2 = elem_type(G)
    return GradedVectorSpaces{F,G2}(K,G)
end

# function (Vec::GradedVectorSpaces{T})(V::FreeModule{T}) where T <: FieldElem
#     return GradedVectorSpaceObject(V,Vec)
# end

function GradedVectorSpaceObject(V::Dict{G,VectorSpaceObject{T}}) where {T,G}
    Vec = GradedVectorSpaces{T,G}(base_ring(V[V.keys[1]]), parent(V.keys[1]))
    return GradedVectorSpaceObject{T,G}(V,Vec)
end

function GradedVectorSpaceObject(V::Pair{G,VectorSpaceObject{T}}...) where {G,T}
    return GradedVectorSpaceObject(Dict([v for v in V]...))
end

function GradedVectorSpaceMorphism(D::GradedVectorSpaceObject{T,G},
        C::GradedVectorSpaceObject{T,G}, m::Dict{G,VectorSpaceMorphism{T}}) where {T,G}
    if parent(D) != parent(C)
        throw(ErrorException("Missmatching parents."))
    elseif !(*([size(m[g].m) == (dim(D[g]),dim(C[g])) for g in keys(m)]...))
        throw(ErrorException("Mismatching dimensions"))
    else
        return GradedVectorSpaceMorphism{T,G}(m,D,C)
    end
end

function GradedVectorSpaceMorphism(D::GradedVectorSpaceObject{T,G}, C::GradedVectorSpaceObject{T,G},
        m::Pair{G,VectorSpaceMorphism{T}}) where {T,G}
    return GradedVectorSpaceMorphism(D,C, Dict(m...))
end
#-----------------------------------------------------------------
#   Pretty Printing
#-----------------------------------------------------------------

function Base.show(io::IO, Vec::GradedVectorSpaces)
    println(io, """
    Category of finite dimensional graded vector spaces over $(base_ring(Vec))
    with simple objects in $(Vec.base_group).""")
end

function Base.show(io::IO, V::GradedVectorSpaceObject)
    print(io, """
    Graded Vector Space of dimension $(dim(V)) over $(base_ring(V)). Direct sum of \n""")
    for (k,v) in V.V
        println(io, "   $k ⋅ $v")
    end
end

#-----------------------------------------------------------------
#   Functionality
#-----------------------------------------------------------------

base_ring(Vec::GradedVectorSpaces) = Vec.base_ring
base_ring(V:: GradedVectorSpaceObject) = V.parent.base_ring

dim(V::GradedVectorSpaceObject) = sum([dim(v) for (g,v) in V.V])

function simples(Vec::GradedVectorSpaces)
    F = base_ring(Vec)
    S = VectorSpaceObject(F,1)
    return [GradedVectorSpaceObject(g => S) for g in Vec.base_group]
end

function one(Vec::GradedVectorSpaces)
    F = base_ring(Vec)
    S = VectorSpaceObject(F,1)
    return GradedVectorSpaceObject(one(Vec.base_group) => S)
end

parent(V::GradedVectorSpaceObject) = V.parent

getindex(V::GradedVectorSpaceObject, x) = x ∈ keys(V.V) ? V.V[x] : VectorSpaceObject(base_ring(V),0)

function getindex(f::GradedVectorSpaceMorphism, x)
    if x ∈ keys(f.m)
        return f.m[x]
    end
    D = domain(f)[x]
    C = codomain(f)[x]
    return VectorSpaceMorphism(D,C,matrix(base_ring(D),zeros(Int64,dim(D),dim(C))))
end

#-----------------------------------------------------------------
#   Functionality: Direct Sums
#-----------------------------------------------------------------

function dsum(V::GradedVectorSpaceObject{T,G}, W::GradedVectorSpaceObject{T,G}) where {T,G}
    Wm,Vm = Dict{G,VectorSpaceMorphism{T}}(),Dict{G,VectorSpaceMorphism{T}}()
    Z = Dict{G,VectorSpaceObject{T}}()

    for g ∈ setdiff(keys(V.V), keys(W.V))
        Z[g] = V[g]
        Vm[g] = id(V[g])
    end
    for g ∈ keys(V.V) ∩ keys(W.V)
        S,i = dsum(V[g], W[g])
        Z[g] = S
        Vm[g] = i[1]
        Wm[g] = i[2]
    end
    for g ∈ setdiff(keys(W.V),keys(V.V))
        Z[g] = W[g]
        Wm[g] = id(W[g])
    end

    VZ = GradedVectorSpaceObject(Z)
    mV = GradedVectorSpaceMorphism(V,VZ, Vm)
    mW = GradedVectorSpaceMorphism(W,VZ, Wm)
    return VZ, [mV,mW]
end

⊕(V::GradedVectorSpaceObject{T,G}, W::GradedVectorSpaceObject{T,G}) where {T,G} = dsum(V,W)

function dsum(f::GradedVectorSpaceMorphism{T,G}, g::GradedVectorSpaceMorphism{T,G}) where {T,G}
    D,_ = domain(f) ⊕ domain(g)
    C,_ = codomain(f) ⊕ codomain(g)

    m = Dict{G,VectorSpaceMorphism{T}}()

    for x ∈ keys(f.m) ∩ keys(g.m)
        m[x] = f[x] ⊕ g[x]
    end
    for x ∈ setdiff(keys(f.m), keys(g.m))
        m[x] = f[x]
    end
    for x ∈ setdiff(keys(g.m), keys(f.m))
        m[x] = g[x]
    end

    return GradedVectorSpaceMorphism(D,C,m)
end


#-----------------------------------------------------------------
#   Functionality: Tensor Products
#-----------------------------------------------------------------

function tensor_product(V::GradedVectorSpaceObject{T,G}, W::GradedVectorSpaceObject{T,G}) where {T,G}
    if parent(V) != parent(W)
        throw(ErrorException("Mismatching parents."))
    end
    Z = Dict{G,VectorSpaceObject{T}}()

    for g ∈ keys(V.V), h ∈ keys(W.V)
        if g*h ∈ keys(Z)
            Z[g*h] = dsum(Z[g*h],V[g]⊗W[h])[1]
        else
            Z[g*h] = V[g]⊗W[h]
        end
    end
    return GradedVectorSpaceObject(Z)
end

⊗(X::GradedVectorSpaceObject{T}, Y::GradedVectorSpaceObject{T}) where {T} = tensor_product(X,Y)
#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::GradedVectorSpaceMorphism{T,G}...) where {T,G}
    if [domain(f[i]) == codomain(f[i-1]) for i ∈ 2:length(f)] != trues(length(f)-1)
        throw(ErrorException("Morphisms not compatible"))
    end
    m = Dict(g => compose([f[i][g] for i ∈ 1:length(f)]...) for g ∈ ∪([keys(f[i].m) for i ∈ 1:length(f)]...))
    return GradedVectorSpaceMorphism(domain(f[1]), codomain(f[end]), m)
end

∘(f::GradedVectorSpaceMorphism...) = compose(reverse(f)...)

function ==(f::GradedVectorSpaceMorphism, g::GradedVectorSpaceMorphism)
    a = domain(f) == domain(g)
    b = codomain(f) == codomain(g)
    c = f.m == g.m
    return a && b && c
end

function id(X::GradedVectorSpaceObject{T,G}) where {T,G}
    n = [dim(X[g]) for g ∈ keys(X.V)]
    m = Dict(g => id(X[g]) for g ∈ keys(X.V))
    return GradedVectorSpaceMorphism(X,X, m)
end


#------------------------------------------------------------------------------
# Associators
#------------------------------------------------------------------------------


function associator(X::GradedVectorSpaceObject{T,G}, Y::GradedVectorSpaceObject{T,G}, Z::GradedVectorSpaceObject{T,G}) where T,G
    #todo
end
