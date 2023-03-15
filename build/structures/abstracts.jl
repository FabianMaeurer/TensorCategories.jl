#------------------------------------------------------------------------
#   Structs for categories
#------------------------------------------------------------------------

abstract type Category end

abstract type CategoryObject end

abstract type CategoryMorphism end


"""
    VectorSpaceCategoryObject

An object in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceCategoryObject <: CategoryObject end

"""
    VectorSpaceCategoryMorphism

A morphism in the category of finite dimensional vector spaces.
"""
abstract type VectorSpaceCategoryMorphism <: CategoryMorphism end

abstract type CategoryHomSet end

abstract type AbstractCategoryHomSpace <: VectorSpaceCategoryObject end

struct CategoryHomSpace <: AbstractCategoryHomSpace
    X::CategoryObject
    Y::CategoryObject
    basis::Vector{<:CategoryMorphism}
    parent
end

domain(m::CategoryMorphism) = m.domain
codomain(m::CategoryMorphism) = m.codomain

"""
    parent(X::CategoryObject)

Return the parent category of the object X.
"""
parent(X::CategoryObject) = X.parent

"""
    function parent(f::CategoryMorphism)

Return the parent category of ``f``.
"""
parent(f::CategoryMorphism) = parent(domain(f))

"""
    base_ring(X::CategoryObject)

Return the base ring ```k``` of the ```k```-linear parent category of ```X```.
"""
base_ring(X::CategoryObject) = base_ring(parent(X))
base_ring(X::CategoryMorphism) = base_ring(parent(domain(X)))

"""
    base_ring(C::Category)

Return the base ring ```k```of the ```k```-linear category ```C```.
"""
base_ring(C::Category) = C.base_ring

base_group(C::Category) = C.base_group
base_group(X::CategoryObject) = parent(X).base_group

#---------------------------------------------------------
#   Direct Sums, Products, Coproducts
#---------------------------------------------------------

function ⊕(T::Tuple{S,Vector{R},Vector{R2}},X::S1) where {S <: CategoryObject,S1 <: CategoryObject, R <: CategoryMorphism, R2 <: CategoryMorphism}
    Z,ix,px = direct_sum(T[1],X)
    incl = vcat([ix[1] ∘ t for t in T[2]], ix[2:2])
    proj = vcat([t ∘ px[1] for t in T[3]], px[2:2])
    return Z, incl, proj
end

⊕(X::S1,T::Tuple{S,Vector{R}, Vector{R2}}) where {S <: CategoryObject,S1 <: CategoryObject, R <: CategoryMorphism, R2 <: CategoryMorphism} = ⊕(T,X)

function direct_sum(X::CategoryObject...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y ∈ X[2:end]
        Z = direct_sum(Z,Y)
    end
    return Z
end

function direct_sum(X::CategoryObject...)
    if length(X) == 1
        return X[1], [id(X[1])],[id(X[1])]
    end
    Z,ix,px = direct_sum(X[1],X[2])
    for Y in X[3:end]
        Z,ix,px = ⊕((Z,ix,px),Y)
    end
    return Z,ix,px
end

function direct_sum(f::CategoryMorphism...)
    g = f[1]

    for h ∈ f[2:end]
        g = g ⊕ h
    end
    return g
end

function ×(T::Tuple{S,Vector{R}},X::S1) where {S <: CategoryObject,S1 <: CategoryObject, R <: CategoryMorphism}
    Z,px = product(T[1],X)
    m = vcat([t ∘ px[1] for t in T[2]], px[2])
    return Z, m
end

×(X::S1,T::Tuple{S,Vector{R}}) where {S <: CategoryObject,S1 <: CategoryObject, R <: CategoryMorphism} = ×(T,X)

function product(X::CategoryObject...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y ∈ X[2:end]
        Z = product(Z,Y)
    end
    return Z
end

function product(X::CategoryObject...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,px = product(X[1],X[2], true)
    for Y in X[3:end]
        Z,px = ×((Z,px),Y)
    end
    return Z,px
end

function ∐(T::Tuple{S,Vector{R}},X::S1) where {S <: CategoryObject,S1 <: CategoryObject, R <: CategoryMorphism}
    Z,px = coproduct(T[1],X)
    m = vcat([px[1] ∘ t for t in T[2]], px[2])
    return Z, m
end

∐(X::S1,T::Tuple{S,Vector{R}}) where {S <: CategoryObject,S1 <: CategoryObject, R <: CategoryMorphism} = ∐(T,X)

function coproduct(X::CategoryObject...)
    if length(X) == 0 return nothing end
    Z = X[1]
    for Y in X[2:end]
        Z = coproduct(Z,Y)
    end
    return Z
end

function coproduct(X::CategoryObject...)
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
    ×(X::CategoryObject...)

Return the product CategoryObject and an array containing the projection morphisms.
"""
×(X::CategoryObject...) = product(X...)

"""
    ∐(X::CategoryObject...)

Return the coproduct CategoryObject and an array containing the injection morphisms.
"""
∐(X::CategoryObject...) = coproduct(X...)

"""
    ⊕(X::CategoryObject...)

Return the direct sum CategoryObject and arrays containing the injection and projection
morphisms.
"""

⊕(X::CategoryObject...) = direct_sum(X...)

⊕(X::CategoryMorphism...) = direct_sum(X...)

"""
    ⊗(X::CategoryObject...)

Return the tensor product object.
"""
⊗(X::CategoryObject...) = tensor_product(X...)

"""
    ^(X::CategoryObject, n::Integer)

Return the n-fold product object ```X^n```.
"""
^(X::CategoryObject,n::Integer) = n == 0 ? zero(parent(X)) : product([X for i in 1:n]...)

^(X::CategoryMorphism,n::Integer) = n == 0 ? zero_morphism(zero(parent(domain(X))), zero(parent(domain(X)))) : direct_sum([X for i in 1:n]...)
"""
    ⊗(f::CategoryMorphism, g::CategoryMorphism)

Return the tensor product morphism of ```f```and ```g```.
"""
⊗(f::CategoryMorphism, g::CategoryMorphism) = tensor_product(f,g)


direct_sum(X::T) where T <: Union{Vector,Tuple} = direct_sum(X...)
product(X::T) where T <: Union{Vector,Tuple} = product(X...)
coproduct(X::T) where T <: Union{Vector,Tuple} = coproduct(X...)

product(X::CategoryObject,Y::CategoryObject) = direct_sum(X,Y)
coproduct(X::CategoryObject, Y::CategoryObject) = direct_sum(X,Y)

#---------------------------------------------------------
#   Horizontal and Vertical direct sums
#---------------------------------------------------------

"""
    function horizontal_direct_sum(f::CategoryMorphism, g::CategoryMorphism)

Return the sum of ``f:X → Z``, ``g:Y → Z`` as ``f+g:X⊕Y → Z.
"""
function horizontal_direct_sum(f::CategoryMorphism, g::CategoryMorphism)
    #@assert codomain(f) == codomain(g) "Codomains do not coincide"
    sum = f ⊕ g
    _,_,(p1,p2) = direct_sum(codomain(f),codomain(g))
    return p1∘sum + p2∘sum
end

function horizontal_direct_sum(f::Vector{M}) where M <: CategoryMorphism
    #@assert codomain(f) == codomain(g) "Codomains do not coincide"
    f_sum = direct_sum(f...)
    _,_,p = direct_sum([codomain(fi) for fi ∈ f]...)
    return sum([p1∘f_sum for p1 ∈ p])
end

"""
    function vertical_direct_sum(f::CategoryMorphism, g::CategoryMorphism)

Return the sum of ``f:X → Y``, ``g:X → Z`` as ``f+g: X → Y⊕Z.
"""
function vertical_direct_sum(f::CategoryMorphism, g::CategoryMorphism)
    #@assert domain(f) == domain(g) "Domains do not coincide"

    sum = f ⊕ g
    _,(i1,i2),_ = direct_sum(domain(f), domain(g))
    return sum∘i1 + sum∘i2
end

function vertical_direct_sum(f::Vector{M}) where M <: CategoryMorphism
    f_sum = direct_sum(f...)
    _,i,_ = direct_sum([domain(fi) for fi ∈ f]...)
    return sum([f_sum∘ix for ix ∈ i])

end

is_simple(X::CategoryObject) = sum([dim(Hom(X,s)) for s ∈ simples(parent(X))]) == 1
#---------------------------------------------------------
#   tensor_product
#---------------------------------------------------------


function tensor_product(X::CategoryObject...)
    if length(X) == 1 return X end

    Z = X[1]
    for Y ∈ X[2:end]
        Z = Z⊗Y
    end
    return Z
end

tensor_product(X::T) where T <: Union{Vector,Tuple} = tensor_product(X...)


"""
    distribute_left(X::RingCategoryObject, Y::RingCategoryObject, Z::RingCategoryObject)

Return the canonical isomorphism ```(X⊕Y)⊗Z → (X⊗Z)⊕(Y⊗Z)```.
"""
function distribute_left(X::CategoryObject, Y::CategoryObject, Z::CategoryObject)
    XY,(ix,iy),(px,py) = direct_sum(X,Y)
    return  vertical_direct_sum(px⊗id(Z), py⊗id(Z))
end

"""
    distribute_left(X::Vector{O}, Z::O) where O <: CategoryObject

Return the canonical isomorphism ```(⨁Xi)⊗Z → ⨁(Xi⊗Z)```.
"""
function distribute_left(X::Vector{O}, Z::O) where O <: CategoryObject
    XY,ix,px = direct_sum(X...)
    return vertical_direct_sum([pi⊗id(Z) for pi ∈ px])
end


"""
    distribute_right(X::RingCategoryObject, Y::RingCategoryObject, Z::RingCategoryObject)

Return the canonical isomorphism ```X⊗(Y⊕Z) → (X⊗Y)⊕(X⊗Z)````
"""
function distribute_right(X::CategoryObject, Y::CategoryObject, Z::CategoryObject)
    XY,(iy,iz),(py,pz) = direct_sum(Y,Z)
    return  vertical_direct_sum(id(X)⊗py, id(X)⊗pz)
end

"""
    distribute_left(X::O, Z::Vector{O}) where O <: CategoryObject

Return the canonical isomorphism ```Z⊗(⨁Xi) → ⨁(Z⊗Xi)```.
"""
function distribute_right(X::O, Z::Vector{O}) where O <: CategoryObject
    XY,ix,px = direct_sum(Z...)
    return vertical_direct_sum([id(X)⊗pi for pi ∈ px])
end

function distribute_left_to_right(X::Vector{T}, Y::Vector{T}) where T <: CategoryObject
    X_sum,ix,px = direct_sum(X...)
    Y_sum,iy,py = direct_sum(Y...)
    Z_sum,iz,pz = direct_sum(Z...)
    direct_sum([(pxk ⊗ pyj ⊗ pzi) ∘ (ixk ⊗ iyj ⊗ izi) for (izi, pzi) ∈ zip(iz,pz), (iyj,pyj) ∈ zip(iy,py), (ixk,pxk) ∈ zip(ix,px)][:]...)
end

function distribute_right_to_left(X::Vector{T}, Y::Vector{T}, Z::Vector{T}) where T <: CategoryObject
    X_sum,ix,px = direct_sum(X...)
    Y_sum,iy,py = direct_sum(Y...)
    Z_sum,iz,pz = direct_sum(Z...)
    direct_sum([(pxk ⊗ (pyj ⊗ pzi)) ∘ (ixk ⊗ (iyj ⊗ izi)) for (izi, pzi) ∈ zip(iz,pz), (iyj,pyj) ∈ zip(iy,py), (ixk,pxk) ∈ zip(ix,px)][:]...)
end

inv_associator(X::CategoryObject, Y::CategoryObject, Z::CategoryObject) = inv(associator(X,Y,Z))


#------------------------------------------------------
#   Abstract Methods
#------------------------------------------------------


function image(f::CategoryMorphism)
    C,c = cokernel(f)
    return kernel(c)
end

∘(f::CategoryMorphism...) = compose(reverse(f)...)

-(f::CategoryMorphism, g::CategoryMorphism) = f + (-1)*g
-(f::CategoryMorphism) = (-1)*f

getindex(C::Category, x::Int) = simples(C)[x]

#=-------------------------------------------------
    Multifusion Categories 
-------------------------------------------------=#

function decompose(C::Category)
    @assert is_multitensor(C)
    one_components = [o for (o,_) in decompose(one(C), simples(C))] 

    if length(one_components) == 1
        return [C]
    end
    S = simples(C)
    structure = [length(filter!(e -> e != zero(C), [𝟙ᵢ⊗s⊗𝟙ⱼ for s ∈ S])) for 𝟙ⱼ ∈ one_components, 𝟙ᵢ ∈ one_components]

    components = []
    comp = [1]
    while Set(vcat(components...)) != Set([i for i ∈ 1:length(one_components)])
        js = findall(e -> e != 0, filter(e -> e != structure[comp[end],comp[end]], structure[:,comp[end]]))
        if length(js) == 0
            components = [components; [comp]]
            k = findfirst(e -> !(e ∈ vcat(components...)), 1:length(one_components))
            if k === nothing
                continue
            end
            comp = [k]
        end
        comp = [comp; js]
    end
    return [RingSubcategory(C,c) for c ∈ components]
end

#-------------------------------------------------------
# Hom Spaces
#-------------------------------------------------------

dim(V::CategoryHomSpace) = length(basis(V))

End(X::CategoryObject) = Hom(X,X)

zero_morphism(C::Category) = zero_morphism(zero(C), zero(C))

Base.iterate(H::AbstractCategoryHomSpace, state = 1) = state > int_dim(H) ? nothing : (basis(H)[state], state + 1)
Base.length(H::AbstractCategoryHomSpace) = int_dim(H)
Base.eltype(::Type{T}) where T <: AbstractCategoryHomSpace = CategoryMorphism 

function (F::Field)(f::CategoryMorphism)
    m = matrix(f)
    if m == zero(parent(m))
        return zero(F)
    end
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
    ind = findfirst(e -> M[e...] != 0 && M[e...] != 0, [(i,j) for i ∈ 1:n, j ∈ 1:m])
    if ind === nothing return false, nothing end
    i,j = Tuple(ind)
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

function express_in_basis(f::T, B::Vector{T}) where T <: CategoryMorphism
    F = base_ring(f)
    B_mat = matrix(F,hcat([[x for x ∈ matrix(b)][:] for b ∈ B]...))
    f_mat = matrix(F, 1, *(size(matrix(f))...), [x for x ∈ matrix(f)][:])

    return [x for x ∈ solve_left(transpose(B_mat),f_mat)][:]
end
#-------------------------------------------------------
# Duals
#-------------------------------------------------------

left_dual(X::CategoryObject) = dual(X)
right_dual(X::CategoryObject) = dual(X)

dual(f::CategoryMorphism) = left_dual(f)

function left_dual(f::CategoryMorphism)
    X = domain(f)
    Y = codomain(f)
    a = ev(Y)⊗id(dual(X))
    b = (id(dual(Y))⊗f)⊗id(dual(X))
    c = inv(associator(dual(Y),X,dual(X)))
    d = id(dual(Y))⊗coev(X)
    (a)∘(b)∘(c)∘(d)
end

tr(f::CategoryMorphism) = left_trace(f)

function left_trace(f::CategoryMorphism)
    V = domain(f)
    W = codomain(f)
    C = parent(V)
    if V == zero(C) || W == zero(C) return zero_morphism(one(C),one(C)) end

    if V == W
        return ev(left_dual(V)) ∘ ((spherical(V)∘f) ⊗ id(left_dual(V))) ∘ coev(V)
    end
    return ev(left_dual(V)) ∘ (f ⊗ id(left_dual(V))) ∘ coev(V)
end

function right_trace(f::CategoryMorphism)
    V = domain(f)
    W = codomain(f)
    dV = right_dual(V)
    _,i = is_isomorphic(left_dual(dV),V)
    _,j = is_isomorphic(right_dual(V), left_dual(right_dual(dV)))
    return (ev(right_dual(dV))) ∘ (j⊗(f∘i)) ∘ coev(right_dual(V))
end

#-------------------------------------------------------
# Spherical structure
#-------------------------------------------------------

function drinfeld_morphism(X::CategoryObject)
     (ev(X)⊗id(dual(dual(X)))) ∘ (braiding(X,dual(X))⊗id(dual(dual(X)))) ∘ (id(X)⊗coev(dual(X)))
 end

dim(X::CategoryObject) = base_ring(X)(tr(spherical(X)))

dim(C::Category) = sum(dim(s)^2 for s ∈ simples(C))
#-------------------------------------------------------
# S-Matrix
#-------------------------------------------------------

function smatrix(C::Category, simples = simples(C))
    @assert is_semisimple(C) "Category has to be semisimple"
    F = base_ring(C)
    m = [tr(braiding(s,t)∘braiding(t,s)) for s ∈ simples, t ∈ simples]
    try
        return matrix(F,[F(n) for n ∈ m])
    catch
    end
    return matrix(F,m)
end

function normalized_smatrix(C::Category, simples = simples(C))
    d = inv(sqrt(dim(C)))
    K = base_ring(C)
    if characteristic(K) == 0
        f = complex_embeddings(K)[1]
        if real(f(d)) < 0
            d = -d
        end
    end
    return d * smatrix(C)
end

function tmatrix(C::Category, simples = simples(C))
    F=base_ring(C)
    T=[1//dim(S)*F(tr(braiding(S,dual(S)))) for S in simples]
    return diagonal_matrix(T)
end

#-------------------------------------------------------
# decomposition morphism
#-------------------------------------------------------

function decompose(X::CategoryObject, S = simples(parent(X)))
    C = parent(X)
    @assert is_semisimple(C) "Category not semisimple"
    dimensions = [int_dim(Hom(s,X)) for s ∈ S]
    return [(s,d) for (s,d) ∈ zip(S,dimensions) if d > 0]
end

function decompose_morphism(X::CategoryObject, S = simples(parent(X)))
    C = parent(X)
    @assert is_semisimple(C) "Semisimplicity required"
    
    if X == zero(C) return id(X), [], [] end

    components = decompose(X,S)
    Z, incl, proj = direct_sum(vcat([[s for _ ∈ 1:d] for (s,d) ∈ components]...)...)

    # temporary solution!
    iso = is_isomorphic(X,Z)[2]
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

function eigenspaces(f::CategoryMorphism)
    @assert domain(f) == codomain(f) "Not an endomorphism"

    #@show factor(minpoly(matrix(f)))
    if base_ring(f) == QQBar
        values = eigenvalues(matrix(f))
    else
        values = collect(keys(eigenspaces(matrix(f))))
    end

    return Dict(λ => kernel(f-λ*id(domain(f)))[1] for λ ∈ values)
end

function simple_subobjects(X::CategoryObject)
    B = basis(End(X))

    if length(B) == 1 return [X] end

    for f ∈ B
        eig_spaces = eigenspaces(f)

        if length(eig_spaces) == 1 && collect(values(eig_spaces))[1] == X
            continue
        end

        simple_subs = vcat([simple_subobjects(K) for (_,K) ∈ eig_spaces]...)

        return unique_simples(simple_subs)
    end
    return [X]
end

function unique_simples(simples::Vector{<:CategoryObject})
    unique_simples = simples[1:1]
    for s ∈ simples[2:end]
        if sum([dim(Hom(s,u)) for u ∈ unique_simples]) == 0
            unique_simples = [unique_simples; s]
        end
    end
    return unique_simples
end

function simples_names(C::Category) 
    @assert is_semisimple(C)
    return ["X$i" for i ∈ 1:length(simples(C))]
end

#=-------------------------------------------------
    Duals in Fusion Categories
-------------------------------------------------=#

function coev(X::CategoryObject)
    if is_simple(X)
        return simple_objects_coev(X)
    end

    C = parent(X)
    𝟙 = one(C)

    summands = vcat([[x for _ ∈ 1:k] for (x,k) ∈ decompose(X)]...)
    dual_summands = dual.(summands)
    d = length(summands)

    c = vertical_direct_sum([i == j ? coev(summands[i]) : zero_morphism(𝟙, summands[j]⊗dual_summands[i]) for j ∈ 1:d, i ∈ 1:d][:])

    distr = direct_sum([distribute_right(x,dual_summands) for x ∈ summands]) ∘ distribute_left(summands, dual(X))

    return distr ∘ c
end

function ev(X::CategoryObject)
    if is_simple(X)
        return simple_objects_ev(X)
    end
    C = parent(X)
    𝟙 = one(C)

    summands = vcat([[x for _ ∈ 1:k] for (x,k) ∈ decompose(X)]...)
    dual_summands = dual.(summands)
    d = length(summands)

    e = horizontal_direct_sum([i == j ? ev(summands[i]) : zero_morphism(dual_summands[j]⊗summands[i], 𝟙)  for j ∈ 1:d, i ∈ 1:d][:])

    distr = direct_sum([distribute_right(x,summands) for x ∈ dual_summands]) ∘ distribute_left(dual_summands, X)

    return e ∘ inv(distr) 
end

function simple_objects_coev(X::CategoryObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    cod = X ⊗ DX

    if X == zero(C) return zero_morphism(one(C), X) end

    return basis(Hom(one(C), cod))[1]
end

function simple_objects_ev(X::CategoryObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    dom = DX ⊗ X

    if X == zero(C) return zero_morphism(X,one(C)) end

    unscaled_ev = basis(Hom(dom,one(C)))[1]

    factor = F((id(X)⊗unscaled_ev)∘associator(X,DX,X)∘(coev(X)⊗id(X)))

    return inv(factor) * unscaled_ev
end

#=-------------------------------------------------
    Frobenius Perron dimension 
-------------------------------------------------=#

function fpdim(X::CategoryObject)
    @assert is_fusion(parent(X))
    S = simples(parent(X))
    n = length(S)

    K = base_ring(X)

 
    A = Array{Int,2}(undef,n,n)
    for i ∈ 1:n
        Y = S[i]
        A[:,i] = [length(basis(Hom(X⊗Y,S[j]))) for j ∈ 1:n]
    end

    if characteristic(K) != 0 || K == QQBar
        K = QQBar
        λ = eigenvalues(matrix(QQ,A),K)
        filter!(e -> isreal(e), λ)
        return findmax(e -> abs(e), λ)
    end


    f = complex_embeddings(K)[1]

    λ = [k for (k,_) ∈ eigenspaces(matrix(K,A))]
    
    filter!(e -> real(f(e)) > 0, λ)

    _,i = findmax(e -> abs(f(e)), λ)
    return λ[i]
end

function fpdim(C::Category)
    @assert is_fusion(C)
    sum(fpdim.(simples(C)).^2)
end


#-------------------------------------------------------
# Misc
#-------------------------------------------------------

*(f::CategoryMorphism, x) = x*f


function is_subobject(X::CategoryObject, Y::CategoryObject)
    @assert parent(X) == parent(Y)
    S = simples(parent(X))

    incl = zero_morphism(X,Y)

    for s ∈ S
        X_s = basis(Hom(X,s))
        s_Y = basis(Hom(s,Y))

        if length(X_s) > length(s_Y) 
            return false, nothing
        elseif length(X_s) > 0
            incl = incl + sum([f∘g for (f,g) ∈ zip(s_Y,X_s)])
        end
    end

    return true,incl
end
