
struct GradedVectorSpaces <: Category
    base_ring::Field
    base_group::GAPGroup
    twist::Cocycle{3}
end

struct GVSObject <: VectorSpaceObject
    parent::GradedVectorSpaces
    V::VectorSpaceObject
    grading::Vector{<:GroupElem}
end

struct GVSMorphism <: VectorSpaceMorphism
    domain::GVSObject
    codomain::GVSObject
    m::MatElem
end

function GradedVectorSpaces(F::Field, G::GAPGroup)
    elems = elements(G)
    GradedVectorSpaces(F,G,trivial_3_cocycle(G,F))
end

function VectorSpaceObject(V::Pair{<:GroupElem, <:VectorSpaceObject}...)
    W = dsum([v for (_,v) ∈ V])
    G = parent(V[1][1])
    elems = elements(G)
    grading = vcat([[g for _ ∈ 1:dim(v)] for (g,v) ∈ V]...)
    C = GradedVectorSpaces(base_ring(W), G)
    return GVSObject(C, W, grading)
end

isfusion(C::GradedVectorSpaces) = true

dim(X::GVSObject) = dim(X.V)
basis(X::GVSObject) = basis(X.V)

function Morphism(X::GVSObject, Y::GVSObject, m::MatElem)
    if !isgraded(X,Y,m)
        throw(ErrorException("Matrix does not define graded morphism"))
    end
    return GVSMorphism(X,Y,m)
end

one(C::GradedVectorSpaces) = GVSObject(C,VectorSpaceObject(base_ring(C),1), [one(base_group(C))])
zero(C::GradedVectorSpaces) = GVSObject(C,VectorSpaceObject(base_ring(C),0), elem_type(base_group(C))[])

function isisomorphic(X::GVSObject, Y::GVSObject)
    b,f = isisomorphic(X.V,Y.V)

    return b && Set(X.grading) == Set(Y.grading) ? (true,Morphism(X,Y,matrix(f))) : (false,nothing)
end
#-----------------------------------------------------------------
#   Functionality: Direct Sums
#-----------------------------------------------------------------

function dsum(X::GVSObject, Y::GVSObject, morphisms::Bool = false)
    W,(ix,iy),(px,py) = dsum(X.V, Y.V, true)
    m,n = dim(X), dim(Y)
    F = base_ring(X)
    grading = [X.grading; Y.grading]
    j = 1

    Z = GVSObject(parent(X), W, grading)

    if morphisms
        ix = Morphism(X,Z,matrix(ix))
        iy = Morphism(Y,Z,matrix(iy))

        px = Morphism(Z,X,matrix(px))
        py = Morphism(Z,Y,matrix(py))

        return Z, [ix,iy], [px,py]
    end

    return Z
end

# function dsum(f::GVSMorphism, g::GVSMorphism)
#     dom = domain(f)⊕domain(g)
#     cod = codomain(f)⊕codomain(g)
#     F = base_ring(f)
#     m1,n1 = size(f.m)
#     m2,n2 = size(g.m)
#     m = [f.m zero(MatrixSpace(F,m1,n2)); zero(MatrixSpace(F,m2,n1)) g.m]
# end

#-----------------------------------------------------------------
#   Functionality: Tensor Product
#-----------------------------------------------------------------

function tensor_product(X::GVSObject, Y::GVSObject)
    W = X.V ⊗ Y.V
    G = base_group(X)
    elems = elements(G)
    grading = vcat([[i*j for i ∈ X.grading] for j ∈ Y.grading]...)
    return GVSObject(parent(X), W, length(grading) == 0 ? elem_type(G)[] : grading)
end

#-----------------------------------------------------------------
#   Functionality: Simple Objects
#-----------------------------------------------------------------

function simples(C::GradedVectorSpaces)
    K = VectorSpaceObject(base_ring(C),1)
    G = base_group(C)
    n = Int(order(G))
    return [GVSObject(C,K,[g]) for g ∈ G]
end

function decompose(V::GVSObject)
    simpls = simples(parent(V))
    return filter(e -> e[2] > 0, [(s, dim(Hom(s,V))) for s ∈ simpls])
end
#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

function kernel(f::GVSMorphism)
    F = base_ring(f)
    G = base_group(domain(f))
    X,Y = domain(f),codomain(f)
    n = dim(X) - rank(f.m)
    m = zero(MatrixSpace(F, n, dim(X)))
    l = 1
    grading = elem_type(G)[]

    for x ∈ unique(domain(f).grading)
        i = findall(e -> e == x, X.grading)
        j = findall(e -> e == x, Y.grading)

        if length(i) == 0 continue end
        if length(j) == 0
            grading = [grading; [x for _ ∈ i]]
            for k in i
                m[l,k] = F(1)
                l = l +1
            end
            continue
        end

        mx = f.m[i,j]

        d,k = kernel(mx, side = :left)
        k = k[1:d,:]

        m[l:l+d-1 ,i] = k
        l = l+d
        grading = [grading; [x for _ ∈ 1:d]]
    end
    K = GVSObject(parent(X), VectorSpaceObject(F,n), grading)
    return K, GVSMorphism(K,domain(f), m)
end

function cokernel(f::GVSMorphism)
    g = GVSMorphism(codomain(f),domain(f),transpose(f.m))
    C,c = kernel(g)
    return C, GVSMorphism(codomain(f),C, transpose(c.m))
end

#-----------------------------------------------------------------
#   Functionality: Associators
#-----------------------------------------------------------------

function associator(X::GVSObject, Y::GVSObject, Z::GVSObject)
    C = parent(X)
    twist = C.twist
    elems = elements(base_group(C))

    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)

    m = one(MatrixSpace(base_ring(X),dim(dom),dim(cod)))

    j = 1
    for x ∈ X.grading, y ∈ Y.grading, z ∈ Z.grading
        m[j,j] = twist(x,y,z)
        j = j+1
    end
    return Morphism(dom,cod,m)
end

#-----------------------------------------------------------------
#   Functionality: Duals
#-----------------------------------------------------------------

function dual(V::GVSObject)
    W = dual(V.V)
    G = base_group(V)
    grading = [inv(j) for j ∈ V.grading]
    return GVSObject(parent(V), W, grading)
end

function ev(V::GVSObject)
    dom = dual(V)⊗V
    cod = one(parent(V))
    elems = elements(base_group(V))
    twist = parent(V).twist
    m = [i == j ? inv(twist(g,inv(g),g)) : 0 for (i,g) ∈ zip(1:dim(V), V.grading), j ∈ 1:dim(V)][:]
    Morphism(dom,cod, matrix(base_ring(V), reshape(m,dim(dom),1)))
end

#-----------------------------------------------------------------
#   Functionality: Hom-Spaces
#-----------------------------------------------------------------

struct GVSHomSpace <: HomSpace
    X::GVSObject
    Y::GVSObject
    basis::Vector{VectorSpaceMorphism}
    parent::VectorSpaces
end

function Hom(V::GVSObject, W::GVSObject)
    G = base_group(V)
    B = VSMorphism[]

    zero_M = MatrixSpace(base_ring(V), dim(V), dim(W))

    for x ∈ unique(V.grading)
        V_grading = findall(e -> e == x, V.grading)
        W_grading = findall(e -> e == x, W.grading)


        for i ∈ V_grading, j ∈ W_grading
            m = zero(zero_M)
            m[i,j] = 1
            B = [B; GVSMorphism(V,W,m)]
        end
    end
    return GVSHomSpace(V,W,B,VectorSpaces(base_ring(V)))
end

function isgraded(X::GVSObject, Y::GVSObject, m::MatElem)
    G = base_group(X)
    for k ∈ 1:order(G)
        i = findall(e -> e == k, X.grading)
        j = findall(e -> e == k, Y.grading)
        for t ∈ i, s ∈ 1:length(j)
            if m[t,s] != 0 && !(s ∈ j)
                return false
            end
        end
    end
    true
end

id(X::GVSObject) = Morphism(X,X,one(MatrixSpace(base_ring(X),dim(X),dim(X))))
#-----------------------------------------------------------------
#   Pretty Printing
#-----------------------------------------------------------------

function show(io::IO, C::GradedVectorSpaces)
    print(io, "Category of G-graded vector spaces over $(base_ring(C)) where G is $(base_group(C))")
end
function show(io::IO, V::GVSObject)
    elems = elements(base_group(V))
    print(io, "Graded vector space of dimension $(dim(V)) with grading\n$(V.grading)")
end
