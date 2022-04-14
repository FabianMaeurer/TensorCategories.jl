
"""
    GradedVectorSpaces{T}(K::S) where T <: FieldElem

The category of finite dimensional vector spaces over K.
"""
struct GradedVectorSpaces<: Category
    base_ring::Field
    base_group::GAPGroup
    twist::Cocycle{3}
end

struct GVSObject <: VectorSpaceObject
    V::Dict{G,VectorSpaceObject} where G <: GroupElem
    parent::GradedVectorSpaces
end

struct GVSMorphism <: VectorSpaceMorphism
    m::Dict{G,VectorSpaceMorphism} where G <: GroupElem
    domain::GVSObject
    codomain::GVSObject
end
#-----------------------------------------------------------------
#   Constructors
#-----------------------------------------------------------------

function GradedVectorSpaces(K::Field,G::GAPGroup)
    return GradedVectorSpaces(K,G,trivial_3_cocycle(G))
end

# function (Vec::GradedVectorSpaces{T})(V::FreeModule{T}) where T <: FieldElem
#     return GradedVectorSpaceObject(V,Vec)
# end

function VectorSpaceObject(V::Dict{G,S}) where {G <: GroupElem, S <: VectorSpaceObject}
    Vec = GradedVectorSpaces(base_ring(V[V.keys[1]]), parent(V.keys[1]))
    return GVSObject{T,G}(V,Vec)
end

function VectorSpaceObject(V::Pair{G,T}...) where {G,T <: VectorSpaceObject}
    return VectorSpaceObject(Dict([v for v in V]...))
end

VectorSpaceObject(x::Base.Generator{<:Dict{<:GroupElem, <:VectorSpaceObject}, <:Function}) = VectorSpaceObject(Dict(x))

"""
    Morphism(D::GVSObject, C::GVSObject, m::Dict)

Morphism between graded vector spaces defined by pairs of elements and morphisms
g => (Dg -> Cg).
"""
function Morphism(D::GVSObject,
        C::GVSObject, m::Dict{G,S}) where {S<:VectorSpaceMorphism,G}
    if parent(D) != parent(C)
        throw(ErrorException("Missmatching parents."))
    # elseif length(m) != 0 && !(*([size(m[g].m) == (dim(D[g]),dim(C[g])) for g in keys(m)]...))
    #     throw(ErrorException("Mismatching dimensions"))
    else
        return GVSMorphism(m,D,C)
    end
end

"""
    Morphism(D::GVSObject, C::GVSObject, m::Pair...)

Morphism between graded vector spaces defined by pairs of elements and morphisms
g => (Dg -> Cg).
"""
function Morphism(D::GVSObject, C::GVSObject,
        m::Pair{G,S}...) where {T,G, S <: VectorSpaceMorphism}
    return Morphism(D,C, Dict(m...))
end

function Morphism(m::Pair{G,S}...) where {G,S}
    dom = VectorSpaceObject([g => domain(f) for (g,f) in m]...)
    cod = VectorSpaceObject([g => codomain(f) for (g,f) in m]...)
    Morphism(dom,cod,m...)
end
#-----------------------------------------------------------------
#   Pretty Printing
#-----------------------------------------------------------------

function Base.show(io::IO, Vec::GradedVectorSpaces)
    print(io, """
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

base_group(Vec::GradedVectorSpaces) = Vec.base_group
base_group(V::GVSObject) = V.parent.base_group

issemisimple(::GradedVectorSpaces) = true

dim(V::GVSObject) = length(V.V) > 0 ? sum([dim(v) for (g,v) in V.V]) : 0

basis(V::GVSObject) = hcat([[(g,v) for v ∈ basis(dual(Vg))] for (g,Vg) ∈ V.V]...)[:]

function matrices(f::GVSMorphism)
    return Dict(g => matrix(m) for (g,m) ∈ f.m)
end

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
    return zero_morphism(D,C)
end

function ==(V::GVSObject, W::GVSObject)
    if keys(V.V) != keys(W.V) return false end
    for g ∈ keys(V.V)
        if V[g] != W[g] return false end
    end
    return true
end

function isisomorphic(V::GVSObject, W::GVSObject)
    if keys(V.V) != keys(W.V) return false, nothing end
    m = []
    for g ∈ keys(V.V)
        b,n = isisomorphic(V[g],W[g])
        if !b return false, nothing end
        m = [m;g => n]
    end
    return true, Morphism(V,W, m...)
end

dual(V::GVSObject) = VectorSpaceObject([g => dual(V[inv(g)]) for g ∈ base_group(V) if inv(g) ∈ keys(V.V)]...)

function ev(V::GVSObject)
    dom = dual(V)⊗V
    cod = one(parent(V))
    o = one(base_group(V))
    ω = parent(V).twist
    m = Morphism(dom[o], cod[o], vcat([ev(V[g]).m for g ∈ keys(V.V)]))
    return Morphism(dom,cod, o => m)
end

function coev(V::GVSObject)
    dom = one(parent(V))
    cod = V ⊗ dual(V)
    o = one(base_group(V))
    m = Morphism(dom[o], cod[o], hcat([coev(V[g]).m for g ∈  keys(V.V)]))
    return Morphism(dom,cod, o => m)
end

spherical(X::GVSObject) = Morphism(X,dual(dual(X)), [g => spherical(V) for (g,V) ∈ X.V]...)

#-----------------------------------------------------------------
#   Functionality: Direct Sums
#-----------------------------------------------------------------

function dsum(V::GVSObject, W::GVSObject, morphisms::Bool = false)
    incl_W, incl_V = Dict{G,VectorSpaceMorphism}(),Dict{G,VectorSpaceMorphism}()
    proj_W, proj_V = Dict{G,VectorSpaceMorphism}(),Dict{G,VectorSpaceMorphism}()
    Z = Dict{G,VectorSpaceObject}()

    for g ∈ setdiff(keys(V.V), keys(W.V))
        Z[g] = V[g]
        incl_V[g] = id(V[g])
        proj_V[g] = id(V[g])
    end
    for g ∈ keys(V.V) ∩ keys(W.V)
        S,i,p = dsum(V[g], W[g],true)
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


    VZ = GVSObject{T,G}(Z,parent(V))

    if !morphisms return VZ end

    mor_incl_V = Morphism(V,VZ, incl_V)
    mor_incl_W = Morphism(W,VZ, incl_W)
    mor_proj_V = Morphism(VZ,V, proj_V)
    mor_proj_W = Morphism(VZ,W, proj_W)
    return VZ, [mor_incl_V,mor_incl_W], [mor_proj_V, mor_proj_W]
end

function dsum(f::GVSMorphism, g::GVSMorphism) where {T,G}
    D = domain(f) ⊕ domain(g)
    C = codomain(f) ⊕ codomain(g)

    m = Dict{G,VectorSpaceMorphism}()

    for x ∈ keys(f.m) ∩ keys(g.m)
        m[x] = f[x] ⊕ g[x]
    end
    for x ∈ setdiff(keys(f.m), keys(g.m))
        m[x] = f[x]
    end
    for x ∈ setdiff(keys(g.m), keys(f.m))
        m[x] = g[x]
    end

    return Morphism(D,C,m)
end

zero(C::GradedVectorSpaces) where {T,G} = GVSObject(Dict{G,VectorSpaceObject}(),C)
#-----------------------------------------------------------------
#   Functionality: Tensor Products
#-----------------------------------------------------------------

function tensor_product(V::GVSObject, W::GVSObject)
    if parent(V) != parent(W)
        throw(ErrorException("Mismatching parents."))
    end
    Z = Dict{G,VectorSpaceObject}()
    grp = base_group(V)

    for g ∈ grp, h ∈ grp
        if g*h ∈ keys(Z)
            Z[g*h] = dsum(Z[g*h],V[g]⊗W[h])
        else
            Z[g*h] = V[g]⊗W[h]
        end
    end
    return GVSObject(Z,parent(V))
end

function tensor_product(f::GVSMorphism, g::GVSMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)
    mors = Dict{G,VSMorphism}()

    grp = base_group(domain(f))
    for x ∈ grp, y ∈ grp
        @show matrix.([f[x],g[x]])
        if (a = x*y) ∈ keys(mors)
            mors[a] = mors[a] ⊕ (f[x] ⊗ g[y])
        else
            mors[a] = (f[x] ⊗ g[y])
        end
    end
    # for g ∈ keys(mors)
    #     dom_left = dim(dom[g]) - dim(domain(mors[g]))
    #     cod_left = dim(cod[g]) - dim(codomain(mors[g]))
    #     Vd = VectorSpaceObject(base_ring(dom), dom_left)
    #     Vc = VectorSpaceObject(base_ring(cod), cod_left)
    #     mors[g] = mors[g] ⊕ zero_morphism(Vd,Vc)
    # end

    return Morphism(dom,cod, [x => mors[x] for x ∈ keys(mors)]...)
end

#-----------------------------------------------------------------
#   Functionality: Morphisms
#-----------------------------------------------------------------

function compose(f::GVSMorphism,g::GVSMorphism)
    if !isisomorphic(domain(g), codomain(f))[1]
        throw(ErrorException("Morphisms not compatible"))
    end
    if length(keys(f.m)∩keys(g.m)) == 0 return zero_morphism(domain(f),codomain(g)) end

    m = Dict(x => compose(f[x],g[x]) for x ∈ keys(f.m)∩keys(g.m))

    return Morphism(domain(f), codomain(g), m)
end


function ==(f::GVSMorphism, g::GVSMorphism)
    a = domain(f) == domain(g)
    b = codomain(f) == codomain(g)
    c = f.m == g.m
    return a && b && c
end

function id(X::GVSObject)
    n = [dim(X[g]) for g ∈ keys(X.V)]
    m = Dict(g => id(X[g]) for g ∈ keys(X.V))
    return Morphism(X,X, m)
end

zero_morphism(X::GVSObject, Y::GVSObject) = Morphism(X,Y,Dict{G,VSMorphism}())


function +(f::GVSMorphism, g::GVSMorphism)
    @assert isisomorphic(domain(f),domain(g))[1] && isisomorphic(codomain(f), codomain(g))[1]

    return Morphism(domain(f),codomain(f), Dict(x => f[x] + g[x] for x ∈ keys(f.m)∪keys(g.m)))
end

function *(x, f::GVSMorphism)
    return Morphism(domain(f),codomain(f), Dict(g => x*m for (g,m) ∈ f.m))
end

function inv(f::GVSMorphism)
    m = Dict(g => inv(m) for (g,m) ∈ f.m)
    return Morphism(codomain(f), domain(f), m)
end

function express_in_basis(f::GVSMorphism, basis::Vector{GVSMorphism})
    F = base_ring(domain(f))
    A = Array{T,2}(undef,length(basis),0)
    grp = base_group(domain(f))
    for g ∈ basis
        B = []
        for x ∈ grp
            if x ∈ keys(g.m)
                B = [B ; [n for n ∈ g[x].m][:]]
            else
                B = [B ; zeros(F,dim(domain(f)[x]), dim(codomain(f)[x]))[:]]
            end
        end
        A = [A B]
    end
    b = []
    for x ∈ grp
        if x ∈ keys(f.m)
            b = [b; [n for n ∈ f[x].m][:]]
        else
            b = [b; zeros(F,dim(domain(f)[x]), dim(codomain(f)[x]))[:]]
        end
    end
    return [i for  i ∈ solve_left(transpose(matrix(F,A)), MatrixSpace(F,1,length(b))(F.(b)))][:]
end
#------------------------------------------------------------------------------
# Associators
#------------------------------------------------------------------------------


function associator(X::GVSObject, Y::GVSObject, Z::GVSObject)
    D = (X⊗Y)⊗Z
    C = X⊗(Y⊗Z)
    Vec = parent(X)
    ω = Vec.twist

    mors = Dict{G,VSMorphism}()

    for (g,Vg) ∈ X.V, (h,Vh) ∈ Y.V, (k,Vk) ∈ Z.V
        if (a = g*h*k) ∈ keys(mors)
            mors[a] = mors[a] ⊕ ω(g,h,k)*associator(Vg,Vh,Vk)
        else
            mors[a] = ω(g,h,k)*associator(Vg,Vh,Vk)
        end
    end
    mors = Dict(g => Morphism(D[g],C[g], mors[g].m) for g in keys(mors))
    return Morphism(D,C,mors)
end


#----------------------------------------------------------------------------
#   Hom Spaces
#----------------------------------------------------------------------------

struct GVSHomSpace <: HomSpace
    X::GVSObject
    Y::GVSObject
    basis::Vector{GVSMorphism}
    parent::VectorSpaces
end

function Hom(X::GVSObject, Y::GVSObject)
    key_elems = keys(X.V)∩keys(Y.V)
    bases = Dict(g => Hom(X[g],Y[g]).basis for g ∈ key_elems)
    basis = Vector{GVSMorphism}()

    for g ∈ key_elems
        for b ∈ bases[g]
            push!(basis, Morphism(X,Y,g => b))
        end
    end
    return GVSHomSpace{T,G}(X,Y,basis,VectorSpaces(base_ring(X)))
end

basis(V::GVSHomSpace) = V.basis

zero(H::GVSHomSpace) = zero_morphism(H.X,H.Y)
