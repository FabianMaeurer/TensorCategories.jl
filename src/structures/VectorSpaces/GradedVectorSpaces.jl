
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

"""
    GradedVectorSpaces(F::Field, G::GAPGroup)

The category of ```G```-graded vector spaces.
"""
function GradedVectorSpaces(F::Field, G::GAPGroup)
    elems = elements(G)
    GradedVectorSpaces(F,G,trivial_3_cocycle(G,F))
end

function GradedVectorSpaces(G::GAPGroup)
    GradedVectorSpaces(QQBar, G)
end
"""
    VectorSpaceObject(V::Pair{<:GroupElem, <:VectorSpaceObject}...)

TBW
"""
function VectorSpaceObject(V::Pair{<:GroupElem, <:VectorSpaceObject}...)
    W,_,_ = direct_sum([v for (_,v) ∈ V])
    G = parent(V[1][1])
    grading = vcat([[g for _ ∈ 1:int_dim(v)] for (g,v) ∈ V]...)
    C = GradedVectorSpaces(base_ring(W), G)
    return GVSObject(C, W, grading)
end

is_fusion(C::GradedVectorSpaces) = true

basis(X::GVSObject) = basis(X.V)

"""
    function Morphism(V::GVSObject, Y::GVSObject, m::MatElem)

Return the morphism ``V → W``defined by ``m``.
"""
function Morphism(X::GVSObject, Y::GVSObject, m::MatElem)
    if !isgraded(X,Y,m)
        throw(ErrorException("Matrix does not define graded morphism"))
    end
    return GVSMorphism(X,Y,m)
end

"""
    function one(C::GradedVectorSpaces)

Return ``k`` as the one dimensional graded vector space.
"""
one(C::GradedVectorSpaces) = GVSObject(C,VectorSpaceObject(base_ring(C),1), [one(base_group(C))])

"""
    function zero(C::GradedVectorSpaces)

Return the zero diemsnional graded vector space.
"""
zero(C::GradedVectorSpaces) = GVSObject(C,VectorSpaceObject(base_ring(C),0), elem_type(base_group(C))[])

"""
    function is_isomorphic(V::GVSObject, W::GVSObject)

Check whether ``V``and ``W``are isomorphic as ``G``-graded vector spaces and return an
isomorphism in the positive case.
"""
function is_isomorphic(X::GVSObject, Y::GVSObject)
    if Set(X.grading) != Set(Y.grading) || !is_isomorphic(X.V,Y.V)[1]
        return false, nothing
    else
        m = zero(MatrixSpace(base_ring(X),int_dim(X),int_dim(Y)))
        for g ∈ X.grading
            i = findall(h -> h == g, X.grading)
            j = findall(h -> h == g, Y.grading)
            for (l,k) ∈ zip(i,j)
                m[l,k] = 1
            end
        end
        return true, Morphism(X,Y,m)
    end
end

function ==(V::GVSObject, W::GVSObject)
    if parent(V) == parent(W) && V.grading == W.grading
        return true
    end
    return false
end

function ==(f::GVSMorphism, g::GVSMorphism)
    return domain(f) == domain(g) && codomain(f) == codomain(g) && matrix(f) == matrix(g)
end

dim(V::GVSObject) = base_ring(V)(tr(id(V)))

function spherical(V::GVSObject) 
    DDV = dual(dual(V))
    dims = filter!(e -> e != 0, collect(matrix(ev(V)))[:])
    m = diagonal_matrix([d for d ∈ dims])
    Morphism(V,DDV,m)
end
#-----------------------------------------------------------------
#   Functionality: Direct Sums
#-----------------------------------------------------------------

"""
    function direct_sum(V::GVSObject, W::GVSObject)

Return the direct sum object ``V⊕W``.
"""
function direct_sum(X::GVSObject, Y::GVSObject)
    W,(ix,iy),(px,py) = direct_sum(X.V, Y.V)
    m,n = int_dim(X), int_dim(Y)
    F = base_ring(X)
    grading = [X.grading; Y.grading]
    j = 1

    Z = GVSObject(parent(X), W, grading)

    ix = Morphism(X,Z,matrix(ix))
    iy = Morphism(Y,Z,matrix(iy))

    px = Morphism(Z,X,matrix(px))
    py = Morphism(Z,Y,matrix(py))

    return Z, [ix,iy], [px,py]
end

# function direct_sum(f::GVSMorphism, g::GVSMorphism)
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

"""
    function tensor_product(V::GVSObject, W::GVSObject)

Return the tensor product ``V⊗W``.
"""
function tensor_product(X::GVSObject, Y::GVSObject)
    W = X.V ⊗ Y.V
    G = base_group(X)
    elems = elements(G)
    grading = vcat([[i*j for i ∈ Y.grading] for j ∈ X.grading]...)
    return GVSObject(parent(X), W, length(grading) == 0 ? elem_type(G)[] : grading)
end

#-----------------------------------------------------------------
#   Functionality: Simple Objects
#-----------------------------------------------------------------

"""
    function simples(C::GradedVectorSpaces)

Return a vector containing the simple objects of ``C``.
"""
function simples(C::GradedVectorSpaces)
    K = VectorSpaceObject(base_ring(C),1)
    G = base_group(C)
    n = Int(order(G))
    return [GVSObject(C,K,[g]) for g ∈ G]
end

"""
    function decompose(V::GVSObject)

Return a vector with the simple objects together with their multiplicities ``[V:Xi]``.
"""
function decompose(V::GVSObject)
    simpls = simples(parent(V))
    return filter(e -> e[2] > 0, [(s, int_dim(Hom(s,V))) for s ∈ simpls])
end
#-----------------------------------------------------------------
#   Functionality: (Co)Kernel
#-----------------------------------------------------------------

"""
    function kernel(f::GVSMorphism)

Return the graded vector space kernel of ``f``.
"""
function kernel(f::GVSMorphism)
    F = base_ring(f)
    G = base_group(domain(f))
    X,Y = domain(f),codomain(f)
    n = int_dim(X) - rank(f.m)
    m = zero(MatrixSpace(F, n, int_dim(X)))
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

"""
    function cokernel(f::GVSMorphism)

Return the graded vector space cokernel of ``f``.
"""
function cokernel(f::GVSMorphism)
    g = GVSMorphism(codomain(f),domain(f),transpose(f.m))
    C,c = kernel(g)
    return C, GVSMorphism(codomain(f),C, transpose(c.m))
end

#-----------------------------------------------------------------
#   Functionality: Associators
#-----------------------------------------------------------------

"""
    function associator(U::GVSObject, V::GVSObject, W::GVSObject)

return the associator isomorphism ``(U⊗V)⊗W → U⊗(V⊗W)``.
"""
function associator(X::GVSObject, Y::GVSObject, Z::GVSObject)
    C = parent(X)
    twist = C.twist
    elems = elements(base_group(C))

    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)

    m = one(MatrixSpace(base_ring(X),int_dim(dom),int_dim(cod)))

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

"""
    function dual(V::GVSObject)

Return the graded dual vector space of ``V``.
"""
function dual(V::GVSObject)
    W = dual(V.V)
    G = base_group(V)
    grading = [inv(j) for j ∈ V.grading]
    return GVSObject(parent(V), W, grading)
end

"""
    function ev(V::GVSObject)

Return the evaluation map ``V*⊗V → 𝟙``.
"""
function ev(V::GVSObject)
    dom = dual(V)⊗V
    cod = one(parent(V))
    elems = elements(base_group(V))
    twist = parent(V).twist
    m = [i == j ? inv(twist(g,inv(g),g)) : 0 for (i,g) ∈ zip(1:int_dim(V), V.grading), j ∈ 1:int_dim(V)][:]
    Morphism(dom,cod, matrix(base_ring(V), reshape(m,int_dim(dom),1)))
end

#-----------------------------------------------------------------
#   Functionality: Hom-Spaces
#-----------------------------------------------------------------

struct GVSCategoryHomSpace <: AbstractCategoryHomSpace
    X::GVSObject
    Y::GVSObject
    basis::Vector{VectorSpaceMorphism}
    parent::VectorSpaces
end

"""
    function Hom(V::GVSObject, W::GVSObject)

Return the space of morphisms between graded vector spaces ``V`` and ``W``.
"""
function Hom(V::GVSObject, W::GVSObject)
    G = base_group(V)
    B = VSMorphism[]

    zero_M = MatrixSpace(base_ring(V), int_dim(V), int_dim(W))

    for x ∈ unique(V.grading)
        V_grading = findall(e -> e == x, V.grading)
        W_grading = findall(e -> e == x, W.grading)


        for i ∈ V_grading, j ∈ W_grading
            m = zero(zero_M)
            m[i,j] = 1
            B = [B; GVSMorphism(V,W,m)]
        end
    end
    return GVSCategoryHomSpace(V,W,B,VectorSpaces(base_ring(V)))
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

is_simple(V::VectorSpaceObject) = dim(V) == 1


# """
#     function id(V::GVSObject)

# description
# """
# id(X::GVSObject) = Morphism(X,X,one(MatrixSpace(base_ring(X),dim(X),dim(X))))
#-----------------------------------------------------------------
#   Pretty Printing
#-----------------------------------------------------------------

function show(io::IO, C::GradedVectorSpaces)
    print(io, "Category of G-graded vector spaces over $(base_ring(C)) where G is $(base_group(C))")
end
function show(io::IO, V::GVSObject)
    elems = elements(base_group(V))
    print(io, "Graded vector space of dimension $(int_dim(V)) with grading\n$(V.grading)")
end
