
"""
    GradedVectorSpaces{T}(K::S) where T <: FieldElem

The category of finite dimensional vector spaces over K.
"""
struct GradedVectorSpaces{T,G} <: TensorCategory{T}
    base_ring::Field
    base_group::GAPGroup
    twist::Cocycle{3}
end

struct GVSObject{T,G} <: VectorSpaceObject{T}
    V::Dict{G,VectorSpaceObject{T}}
    parent::GradedVectorSpaces{T,G}
end

struct GVSMorphism{T,G} <: VectorSpaceMorphism{T}
    m::Dict{G,VectorSpaceMorphism{T}}
    domain::GVSObject{T,G}
    codomain::GVSObject{T,G}
end
#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------

function GradedVectorSpaces(K::T,G::S,twist::Cocycle{3} = trivial_3_cocycle(G)) where {T <: Field, S <: GAPGroup}
    F = elem_type(K)
    G2 = elem_type(G)
    return GradedVectorSpaces{F,G2}(K,G,twist)
end

# function (Vec::GradedVectorSpaces{T})(V::FreeModule{T}) where T <: FieldElem
#     return GradedVectorSpaceObject(V,Vec)
# end

function VectorSpaceObject(V::Dict{G,S}) where {T, S <: VectorSpaceObject{T},G}
    Vec = GradedVectorSpaces(base_ring(V[V.keys[1]]), parent(V.keys[1]))
    return GVSObject{T,G}(V,Vec)
end

function VectorSpaceObject(V::Pair{G,T}...) where {G,T <: VectorSpaceObject}
    return VectorSpaceObject(Dict([v for v in V]...))
end

function VectorSpaceMorphism(D::GVSObject{T,G},
        C::GVSObject{T,G}, m::Dict{G,S}) where {T,S<:VectorSpaceMorphism,G}
    if parent(D) != parent(C)
        throw(ErrorException("Missmatching parents."))
    elseif !(*([size(m[g].m) == (dim(D[g]),dim(C[g])) for g in keys(m)]...))
        throw(ErrorException("Mismatching dimensions"))
    else
        return GVSMorphism{T,G}(m,D,C)
    end
end

function VectorSpaceMorphism(D::GVSObject{T,G}, C::GVSObject{T,G},
        m::Pair{G,S}...) where {T,G, S <: VectorSpaceMorphism{T}}
    return VectorSpaceMorphism(D,C, Dict(m...))
end
#-----------------------------------------------------------------
#   Pretty Printing
#-----------------------------------------------------------------

function Base.show(io::IO, Vec::GradedVectorSpaces)
    println(io, """
    Category of finite dimensional graded vector spaces over $(base_ring(Vec))
    with simple objects in $(Vec.base_group).""")
end

function Base.show(io::IO, V::GVSObject)
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
base_ring(V:: GVSObject) = V.parent.base_ring

issemisimple(::GradedVectorSpaces) = true

dim(V::GVSObject) = sum([dim(v) for (g,v) in V.V])

function simples(Vec::GradedVectorSpaces)
    F = base_ring(Vec)
    S = VectorSpaceObject(F,1)
    return [VectorSpaceObject(g => S) for g in Vec.base_group]
end

function one(Vec::GradedVectorSpaces)
    F = base_ring(Vec)
    S = VectorSpaceObject(F,1)
    return VectorSpaceObject(one(Vec.base_group) => S)
end

function decompose(V::GVSObject)
    F = base_ring(V)
    return [(VectorSpaceObject(g => one(VectorSpaces(F))),dim(v)) for (g,v) in V.V]
end

parent(V::GVSObject) = V.parent

getindex(V::GVSObject, x) = x ∈ keys(V.V) ? V.V[x] : VectorSpaceObject(base_ring(V),0)

function getindex(f::GVSMorphism, x)
    if x ∈ keys(f.m)
        return f.m[x]
    end
    D = domain(f)[x]
    C = codomain(f)[x]
    return VectorSpaceMorphism(D,C,matrix(base_ring(D),zeros(Int64,dim(D),dim(C))))
end

function ==(V::GVSObject, W::GVSObject)
    if keys(V.V) != keys(W.V) return false end
    for g ∈ keys(V.V)
        if V[g] != W[g] return false end
    end
    return true
end

function isisomorphic(V::GVSObject{T,G}, W::GVSObject{T,G}) where {T,G}
    if keys(V.V) != keys(W.V) return false, nothing end
    m = []
    for g ∈ keys(V.V)
        b,n = isisomorphic(V[g],W[g])
        if !b return false, nothing end
        m = [m;g => n]
    end
    return true, Morphism(V,W, m...)
end

#-----------------------------------------------------------------
#   Functionality: Direct Sums
#-----------------------------------------------------------------

function dsum(V::GVSObject{T,G}, W::GVSObject{T,G}, morphisms = false) where {T,G}
    incl_W, incl_V = Dict{G,VectorSpaceMorphism{T}}(),Dict{G,VectorSpaceMorphism{T}}()
    proj_W, proj_V = Dict{G,VectorSpaceMorphism{T}}(),Dict{G,VectorSpaceMorphism{T}}()
    Z = Dict{G,VectorSpaceObject{T}}()

    for g ∈ setdiff(keys(V.V), keys(W.V))
        Z[g] = V[g]
        incl_V[g] = id(V[g])
        proj_V[g] = id(V[g])
    end
    for g ∈ keys(V.V) ∩ keys(W.V)
        S,i,p = dsum(V[g], W[g])
        Z[g] = S
        incl_V[g] = i[1]
        incl_W[g] = i[2]
        proj_V[g] = p[1]
        proj_W[g] = p[2]
    end
    for g ∈ setdiff(keys(W.V),keys(V.V))
        Z[g] = W[g]
        incl_W[g] = id(W[g])
        proj_W[g] = id(W[g])
    end

    if !morphisms return VZ end

    VZ = VectorSpaceObject(Z)
    mor_incl_V = VectorSpaceMorphism(V,VZ, incl_V)
    mor_incl_W = VectorSpaceMorphism(W,VZ, incl_W)
    mor_proj_V = VectorSpaceMorphism(VZ,V, proj_V)
    mor_proj_W = VectorSpaceMorphism(VZ,W, proj_W)
    return VZ, [mor_incl_V,mor_incl_W], [mor_proj_V, mor_proj_W]
end

function dsum(f::GVSMorphism{T,G}, g::GVSMorphism{T,G}) where {T,G}
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

    return VectorSpaceMorphism(D,C,m)
end


#-----------------------------------------------------------------
#   Functionality: Tensor Products
#-----------------------------------------------------------------

function tensor_product(V::GVSObject{T,G}, W::GVSObject{T,G}) where {T,G}
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
    return VectorSpaceObject(Z)
end

#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::GVSMorphism{T,G}...) where {T,G}
    if [domain(f[i]) == codomain(f[i-1]) for i ∈ 2:length(f)] != trues(length(f)-1)
        throw(ErrorException("Morphisms not compatible"))
    end
    m = Dict(g => compose([f[i][g] for i ∈ 1:length(f)]...) for g ∈ ∪([keys(f[i].m) for i ∈ 1:length(f)]...))
    return VectorSpaceMorphism(domain(f[1]), codomain(f[end]), m)
end


function ==(f::GVSMorphism, g::GVSMorphism)
    a = domain(f) == domain(g)
    b = codomain(f) == codomain(g)
    c = f.m == g.m
    return a && b && c
end

function id(X::GVSObject{T,G}) where {T,G}
    n = [dim(X[g]) for g ∈ keys(X.V)]
    m = Dict(g => id(X[g]) for g ∈ keys(X.V))
    return GradedVectorSpaceMorphism(X,X, m)
end


#------------------------------------------------------------------------------
# Associators
#------------------------------------------------------------------------------


function associator(X::GVSObject{T,G}, Y::GVSObject{T,G}, Z::GVSObject{T,G}) where {T,G}
    D = (X⊗Y)⊗Z
    C = X⊗(Y⊗Z)
    mors = Vector{Pair{G,VSMorphism{T}}}() #Dict{G,VSMorphism{T}}()
    Vec = parent(X)
    for g1 ∈ X.V.keys
        for g2 ∈ Y.V.keys
            for g3 ∈ Z.V.keys
                a = g1*g2*g3
                # if false #(a = g1*g2*g3) ∈ keys(mors)
                #     mors[a] = mors[a] ⊕ Vec.twist(g1,g2,g3)*associator(X[g1],Y[g2],Z[g3])
                # else
                #     @show a
                #     push!(mors,a => Vec.twist(g1,g2,g3)*associator(X[g1],Y[g2],Z[g3]))
                # end
            end
        end
    end
    return VectorSpaceMorphism(D,C,mors)
end


#----------------------------------------------------------------------------
#   Hom Spaces
#----------------------------------------------------------------------------

struct GVSHomSpace{T,G} <: HomSpace{T}
    X::GVSObject{T,G}
    Y::GVSObject{T,G}
    basis::Vector{GVSMorphism{T,G}}
    parent::VectorSpaces{T}
end

function Hom(X::GVSObject{T,G}, Y::GVSObject{T,G}) where {T,G}
    key_elems = keys(X.V)∩keys(Y.V)
    bases = Dict(g => Hom(X[g],Y[g]).basis for g ∈ key_elems)
    basis = Vector{GVSMorphism{T,G}}()

    for g ∈ key_elems
        for b ∈ bases[g]
            push!(basis, VectorSpaceMorphism(X,Y,g => b))
        end
    end
    return GVSHomSpace{T,G}(X,Y,basis,VectorSpaces(base_ring(X)))
end

basis(V::GVSHomSpace) = V.basis
