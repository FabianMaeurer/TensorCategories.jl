#------------------------------------------------------------------------
#   Structs for categories
#------------------------------------------------------------------------

abstract type Category end

abstract type Object end

abstract type Morphism end


"""
    VectorSpaceObject

An object in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceObject <: Object end

"""
    VectorSpaceMorphism

A morphism in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceMorphism <: Morphism end

abstract type HomSet end

abstract type AbstractHomSpace <: VectorSpaceObject end

struct HomSpace <: AbstractHomSpace
    X::Object
    Y::Object
    basis::Vector{<:Morphism}
    parent
end

domain(m::Morphism) = m.domain
codomain(m::Morphism) = m.codomain

"""
    parent(X::Object)

Return the parent category of the object X.
"""
parent(X::Object) = X.parent

"""
    function parent(f::Morphism)

Return the parent category of ``f``.
"""
parent(f::Morphism) = parent(domain(f))

"""
    base_ring(X::Object)

Return the base ring ```k``` of the ```k```-linear parent category of ```X```.
"""
base_ring(X::Object) = base_ring(parent(X))
base_ring(X::Morphism) = parent(domain(X)).base_ring

"""
    base_ring(C::Category)

Return the base ring ```k```of the ```k```-linear category ```C```.
"""
base_ring(C::Category) = C.base_ring

base_group(C::Category) = C.base_group
base_group(X::Object) = parent(X).base_group

#---------------------------------------------------------
#   Direct Sums, Products, Coproducts
#---------------------------------------------------------

function ⊕(T::Tuple{S,Vector{R},Vector{R2}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism, R2 <: Morphism}
    Z,ix,px = dsum(T[1],X,true)
    incl = vcat([ix[1] ∘ t for t in T[2]], ix[2:2])
    proj = vcat([t ∘ px[1] for t in T[3]], px[2:2])
    return Z, incl, proj
end

⊕(X::S1,T::Tuple{S,Vector{R}, Vector{R2}}) where {S <: Object,S1 <: Object, R <: Morphism, R2 <: Morphism} = ⊕(T,X)

function dsum(X::Object...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y ∈ X[2:end]
        Z = dsum(Z,Y)
    end
    return Z
end

function dsum_with_morphisms(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])],[id(X[1])]
    end
    Z,ix,px = dsum(X[1],X[2],true)
    for Y in X[3:end]
        Z,ix,px = ⊕((Z,ix,px),Y)
    end
    return Z,ix,px
end

function dsum(f::Morphism...)
    g = f[1]

    for h ∈ f[2:end]
        g = g ⊕ h
    end
    return g
end

function ×(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = product(T[1],X)
    m = vcat([t ∘ px[1] for t in T[2]], px[2])
    return Z, m
end

×(X::S1,T::Tuple{S,Vector{R}}) where {S <: Object,S1 <: Object, R <: Morphism} = ×(T,X)

function product(X::Object...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y ∈ X[2:end]
        Z = product(Z,Y)
    end
    return Z
end

function product_with_morphisms(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,px = product(X[1],X[2], true)
    for Y in X[3:end]
        Z,px = ×((Z,px),Y)
    end
    return Z,px
end

function ∐(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = coproduct(T[1],X)
    m = vcat([px[1] ∘ t for t in T[2]], px[2])
    return Z, m
end

∐(X::S1,T::Tuple{S,Vector{R}}) where {S <: Object,S1 <: Object, R <: Morphism} = ∐(T,X)

function coproduct(X::Object...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y in X[2:end]
        Z = coproduct(Z,Y)
    end
    return Z
end

function coproduct_with_morphisms(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,ix = coproduct(X[1],X[2])
    for Y in X[3:end]
        Z,ix = ∐((Z,ix),Y)
    end
    return Z,ix
end

"""
    ×(X::Object...)

Return the product Object and an array containing the projection morphisms.
"""
×(X::Object...) = product(X...)

"""
    ∐(X::Object...)

Return the coproduct Object and an array containing the injection morphisms.
"""
∐(X::Object...) = coproduct(X...)

"""
    ⊕(X::Object...)

Return the direct sum Object and arrays containing the injection and projection
morphisms.
"""

⊕(X::Object...) = dsum(X...)

⊕(X::Morphism...) = dsum(X...)

"""
    ⊗(X::Object...)

Return the tensor product object.
"""
⊗(X::Object...) = tensor_product(X...)

"""
    ^(X::Object, n::Integer)

Return the n-fold product object ```X^n```.
"""
^(X::Object,n::Integer) = n == 0 ? zero(parent(X)) : product([X for i in 1:n]...)

^(X::Morphism,n::Integer) = n == 0 ? zero_morphism(zero(parent(domain(X))), zero(parent(domain(X)))) : dsum([X for i in 1:n]...)
"""
    ⊗(f::Morphism, g::Morphism)

Return the tensor product morphism of ```f```and ```g```.
"""
⊗(f::Morphism, g::Morphism) = tensor_product(f,g)


dsum(X::T) where T <: Union{Vector,Tuple} = dsum(X...)
product(X::T) where T <: Union{Vector,Tuple} = product(X...)
coproduct(X::T) where T <: Union{Vector,Tuple} = coproduct(X...)

product(X::Object,Y::Object) = dsum(X,Y)
coproduct(X::Object, Y::Object) = dsum(X,Y)

#---------------------------------------------------------
#   Horizontal and Vertical direct sums
#---------------------------------------------------------

"""
    function horizontal_dsum(f::Morphism, g::Morphism)

Return the sum of ``f:X → Z``, ``g:Y → Z`` as ``f+g:X⊕Y → Z.
"""
function horizontal_dsum(f::Morphism, g::Morphism)
    #@assert codomain(f) == codomain(g) "Codomains do not coincide"
    sum = f ⊕ g
    _,_,(p1,p2) = dsum_with_morphisms(codomain(f),codomain(g))
    return p1∘sum + p2∘sum
end

function horizontal_dsum(f::Vector{M}) where M <: Morphism
    #@assert codomain(f) == codomain(g) "Codomains do not coincide"
    f_sum = dsum(f...)
    _,_,p = dsum_with_morphisms([codomain(fi) for fi ∈ f]...)
    return sum([p1∘f_sum for p1 ∈ p])
end

"""
    function vertical_dsum(f::Morphism, g::Morphism)

Return the sum of ``f:X → Y``, ``g:X → Z`` as ``f+g: X → Y⊕Z.
"""
function vertical_dsum(f::Morphism, g::Morphism)
    #@assert domain(f) == domain(g) "Domains do not coincide"

    sum = f ⊕ g
    _,(i1,i2),_ = dsum_with_morphisms(domain(f), domain(g))
    return sum∘i1 + sum∘i2
end

function vertical_dsum(f::Vector{M}) where M <: Morphism
    f_sum = dsum(f...)
    _,i,_ = dsum_with_morphisms([domain(fi) for fi ∈ f]...)
    return sum([f_sum∘ix for ix ∈ i])

end

#---------------------------------------------------------
#   tensor_product
#---------------------------------------------------------


function tensor_product(X::Object...)
    if length(X) == 1 return X end

    Z = X[1]
    for Y ∈ X[2:end]
        Z = Z⊗Y
    end
    return Z
end

tensor_product(X::T) where T <: Union{Vector,Tuple} = tensor_product(X...)


"""
    distribute_left(X::RingCatObject, Y::RingCatObject, Z::RingCatObject)

Return the canonical isomorphism ```(X⊕Y)⊗Z → (X⊗Z)⊕(Y⊗Z)```.
"""
function distribute_left(X::Object, Y::Object, Z::Object)
    XY,(ix,iy),(px,py) = dsum_with_morphisms(X,Y)
    return  vertical_dsum(px⊗id(Z), py⊗id(Z))
end

"""
    distribute_left(X::Vector{O}, Z::O) where O <: Object

Return the canonical isomorphism ```(⨁Xi)⊗Z → ⨁(Xi⊗Z)```.
"""
function distribute_left(X::Vector{O}, Z::O) where O <: Object
    XY,ix,px = dsum_with_morphisms(X...)
    return vertical_dsum([pi⊗id(Z) for pi ∈ px])
end


"""
    distribute_right(X::RingCatObject, Y::RingCatObject, Z::RingCatObject)

Return the canonical isomorphism ```X⊗(Y⊕Z) → (X⊗Y)⊕(X⊗Z)````
"""
function distribute_right(X::Object, Y::Object, Z::Object)
    XY,(iy,iz),(py,pz) = dsum_with_morphisms(Y,Z)
    return  vertical_dsum(id(X)⊗py, id(X)⊗pz)
end

"""
    distribute_left(X::O, Z::Vector{O}) where O <: Object

Return the canonical isomorphism ```Z⊗(⨁Xi) → ⨁(Z⊗Xi)```.
"""
function distribute_right(X::O, Z::Vector{O}) where O <: Object
    XY,ix,px = dsum_with_morphisms(Z...)
    return vertical_dsum([id(X)⊗pi for pi ∈ px])
end

function distribute_left_to_right(X::Vector{T}, Y::Vector{T}) where T <: Object
    X_sum,ix,px = dsum_with_morphisms(X...)
    Y_sum,iy,py = dsum_with_morphisms(Y...)
    Z_sum,iz,pz = dsum_with_morphisms(Z...)
    dsum([(pxk ⊗ pyj ⊗ pzi) ∘ (ixk ⊗ iyj ⊗ izi) for (izi, pzi) ∈ zip(iz,pz), (iyj,pyj) ∈ zip(iy,py), (ixk,pxk) ∈ zip(ix,px)][:]...)
end

function distribute_right_to_left(X::Vector{T}, Y::Vector{T}, Z::Vector{T}) where T <: Object
    X_sum,ix,px = dsum_with_morphisms(X...)
    Y_sum,iy,py = dsum_with_morphisms(Y...)
    Z_sum,iz,pz = dsum_with_morphisms(Z...)
    dsum([(pxk ⊗ (pyj ⊗ pzi)) ∘ (ixk ⊗ (iyj ⊗ izi)) for (izi, pzi) ∈ zip(iz,pz), (iyj,pyj) ∈ zip(iy,py), (ixk,pxk) ∈ zip(ix,px)][:]...)
end

#------------------------------------------------------
#   Abstract Methods
#------------------------------------------------------
isfusion(C::Category) = false
ismultifusion(C::Category) = isfusion(C)

istensor(C::Category) = isfusion(C)
ismultitensor(C::Category) = ismultifusion(C) || istensor(C)

isring(C::Category) = istensor(C)
ismultiring(C::Category) = ismultitensor(C)

ismonoidal(C::Category) = ismultitensor(C)

isabelian(C::Category) = ismultitensor(C)

isadditive(C::Category) = isabelian(C)

islinear(C::Category) = isabelian(C)

issemisimple(C::Category) = ismultitensor(C)

function image(f::Morphism)
    C,c = cokernel(f)
    return kernel(c)
end

∘(f::Morphism...) = compose(reverse(f)...)

-(f::Morphism, g::Morphism) = f + (-1)*g
-(f::Morphism) = (-1)*f

getindex(C::Category, x::Int) = simples(C)[x]

#-------------------------------------------------------
# Hom Spaces
#-------------------------------------------------------

dim(V::HomSpace) = length(basis(V))

End(X::Object) = Hom(X,X)

zero_morphism(C::Category) = zero_morphism(zero(C), zero(C))

Base.iterate(H::AbstractHomSpace, state = 1) = state > int_dim(H) ? nothing : (basis(H)[state], state + 1)
Base.length(H::AbstractHomSpace) = int_dim(H)
Base.eltype(::Type{T}) where T <: AbstractHomSpace = Morphism 

function (F::Field)(f::Morphism)
    m = matrix(f)
    b,c = is_scalar_multiple(m, matrix(id(domain(f))))
    if b 
        return c
    end
    m = collect(m)[m .!= 0]
    if size(m) == (1,)
        return F(m[1,1])
    end
    @show size(m)
    throw(ErrorException("Cannot convert to element of $F"))
end

function is_scalar_multiple(M::MatElem,N::MatElem)
    n,m = size(M)
    (i,j) = Tuple(findfirst(e -> M[e...] != 0 && M[e...] != 0, [(i,j) for i ∈ 1:n, j ∈ 1:m]))
    if i === nothing return false, nothing end
    k = M[i,j] * inv(N[i,j])
    for (a,b) ∈ zip(M,N)
        if a == b == 0 
            continue
        elseif a == 0 || b == 0 
            return false, nothing
        elseif a * inv(b) != k
            return false, nothing
        end
    end
    return true,k
end

function express_in_basis(f::T, B::Vector{T}) where T <: Morphism
    F = base_ring(f)
    B_mat = matrix(F,hcat([[x for x ∈ matrix(b)][:] for b ∈ B]...))
    f_mat = matrix(F, 1, *(size(matrix(f))...), [x for x ∈ matrix(f)][:])

    return [x for x ∈ solve_left(transpose(B_mat),f_mat)][:]
end
#-------------------------------------------------------
# Duals
#-------------------------------------------------------

left_dual(X::Object) = dual(X)
right_dual(X::Object) = dual(X)

dual(f::Morphism) = left_dual(f)

function left_dual(f::Morphism)
    X = domain(f)
    Y = codomain(f)
    a = ev(Y)⊗id(dual(X))
    b = (id(dual(Y))⊗f)⊗id(dual(X))
    c = inv(associator(dual(Y),X,dual(X)))
    d = id(dual(Y))⊗coev(X)
    (a)∘(b)∘(c)∘(d)
end

tr(f::Morphism) = left_trace(f)

function left_trace(f::Morphism)
    V = domain(f)
    W = codomain(f)
    C = parent(V)
    if V == zero(C) || W == zero(C) return zero_morphism(one(C),one(C)) end

    if V == W
        return ev(left_dual(V)) ∘ ((spherical(V)∘f) ⊗ id(left_dual(V))) ∘ coev(V)
    end
    return ev(left_dual(V)) ∘ (f ⊗ id(left_dual(V))) ∘ coev(V)
end

function right_trace(f::Morphism)
    V = domain(f)
    W = codomain(f)
    dV = right_dual(V)
    _,i = isisomorphic(left_dual(dV),V)
    _,j = isisomorphic(right_dual(V), left_dual(right_dual(dV)))
    return (ev(right_dual(dV))) ∘ (j⊗(f∘i)) ∘ coev(right_dual(V))
end

#-------------------------------------------------------
# Spherical structure
#-------------------------------------------------------

function drinfeld_morphism(X::Object)
     (ev(X)⊗id(dual(dual(X)))) ∘ (braiding(X,dual(X))⊗id(dual(dual(X)))) ∘ (id(X)⊗coev(dual(X)))
 end

dim(X::Object) = base_ring(X)(tr(spherical(X)))

dim(C::Category) = sum(dim(s)^2 for s ∈ simples(C))
#-------------------------------------------------------
# S-Matrix
#-------------------------------------------------------

function smatrix(C::Category, simples = simples(C))
    @assert issemisimple(C) "Category has to be semisimple"
    F = base_ring(C)
    m = [tr(braiding(s,t)∘braiding(t,s)) for s ∈ simples, t ∈ simples]
    try
        return matrix(F,[F(n) for n ∈ m])
    catch
    end
    return matrix(F,m)
end

function tmatrix(C::Category, simples = simples(C))
    F=base_ring(C)
    T=[1//dim(S)*F(tr(braiding(S,dual(S)))) for S in simples]
    return diagonal_matrix(T)
end

#-------------------------------------------------------
# decomposition morphism
#-------------------------------------------------------

function decompose(X::Object, S = simples(parent(X)))
    C = parent(X)
    @assert issemisimple(C) "Category not semisimple"
    dimensions = [int_dim(Hom(X,s)) for s ∈ S]
    return [(s,d) for (s,d) ∈ zip(S,dimensions) if d > 0]
end

function decompose_morphism(X::Object, S = simples(parent(X)))
    C = parent(X)
    @assert issemisimple(C) "Semisimplicity required"
    
    if X == zero(C) return id(X), [], [] end

    components = decompose(X,S)
    Z, incl, proj = dsum_with_morphisms(vcat([[s for _ ∈ 1:d] for (s,d) ∈ components]...)...)

    # temporary solution!
    iso = isisomorphic(X,Z)[2]
    return iso, [inv(iso)∘i for i ∈ incl], [p∘iso for p ∈ proj]

    #----------------------------------
    f = zero_morphism(X,Z)

    for (p,i) ∈ zip(proj, incl)
        g = i∘p
        f = f + g
    end
    return f, incl, proj
end




#-------------------------------------------------------
# Semisimple: Subobjects
#-------------------------------------------------------

function eigenspaces(f::Morphism)
    @assert domain(f) == codomain(f) "Not an endomorphism"

    values = collect(keys(eigenspaces(matrix(f))))

    return Dict(λ => kernel(f-λ*id(domain(f)))[1] for λ ∈ values)
end

function simple_subobjects(X::Object)
    B = basis(End(X))

    if length(B) == 1 return [X] end

    for f ∈ B
        eig_spaces = eigenspaces(f)

        if length(eig_spaces) == 1 
            continue
        end

        simple_subs = vcat([simple_subobjects(K) for (_,K) ∈ eig_spaces]...)

        return unique_simples(simple_subs)
    end
    return [X]
end

function unique_simples(simples::Vector{<:Object})
    unique_simples = simples[1:1]
    for s ∈ simples[2:end]
        if sum([dim(Hom(s,u)) for u ∈ unique_simples]) == 0
            unique_simples = [unique_simples; s]
        end
    end
    return unique_simples
end


#=-------------------------------------------------
    Duals in Fusion Categories
-------------------------------------------------=#

function coev(X::Object)
    if is_simple(X)
        return simple_objects_coev(X)
    end

    C = parent(X)
    𝟙 = one(C)

    summands = vcat([[x for _ ∈ 1:k] for (x,k) ∈ decompose(X)]...)
    dual_summands = dual.(summands)
    d = length(summands)

    c = vertical_dsum([i == j ? coev(summands[i]) : zero_morphism(𝟙, summands[j]⊗dual_summands[i]) for j ∈ 1:d, i ∈ 1:d][:])

    distr = dsum([distribute_right(x,dual_summands) for x ∈ summands]) ∘ distribute_left(summands, dual(X))

    return distr ∘ c
end

function ev(X::Object)
    if is_simple(X)
        return simple_objects_ev(X)
    end
    C = parent(X)
    𝟙 = one(C)

    summands = vcat([[x for _ ∈ 1:k] for (x,k) ∈ decompose(X)]...)
    dual_summands = dual.(summands)
    d = length(summands)

    e = horizontal_dsum([i == j ? ev(summands[i]) : zero_morphism(dual_summands[j]⊗summands[i], 𝟙)  for j ∈ 1:d, i ∈ 1:d][:])

    distr = dsum([distribute_right(x,summands) for x ∈ dual_summands]) ∘ distribute_left(dual_summands, X)

    return e ∘ inv(distr) 
end

function simple_objects_coev(X::Object)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    cod = X ⊗ DX

    if X == zero(C) return zero_morphism(one(C), X) end

    return basis(Hom(one(C), cod))[1]
end

function simple_objects_ev(X::Object)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    dom = DX ⊗ X

    if X == zero(C) return zero_morphism(X,one(C)) end

    unscaled_ev = basis(Hom(dom,one(C)))[1]

    factor = F((id(X)⊗unscaled_ev)∘associator(X,DX,X)∘(coev(X)⊗id(X)))

    return inv(factor) * unscaled_ev
end


#-------------------------------------------------------
# Misc
#-------------------------------------------------------

*(f::Morphism, x) = x*f
