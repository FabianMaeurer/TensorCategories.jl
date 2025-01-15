
isequal(f::Morphism, g::Morphism) = f == g

#=----------------------------------------------------------
    Constructors 
----------------------------------------------------------=##

function morphism(X::T, Y::T, m...) where T <: Object 
    morphism_type(parent(X))(X,Y,m...)
end

#---------------------------------------------------------
#   Direct Sums, Products, Coproducts
#---------------------------------------------------------

direct_sum(X::Vector{<:Object}) = direct_sum(X...)

function ⊕(X::T,Y::T) where {T <: Object}
    return direct_sum(X,Y)[1]
end


function ⊕(T::Tuple{S,Vector{R},Vector{R2}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism, R2 <: Morphism}
    Z,ix,px = direct_sum(T[1],X)
    incl = vcat([ix[1] ∘ t for t in T[2]], ix[2:2])
    proj = vcat([t ∘ px[1] for t in T[3]], px[2:2])
    return Z, incl, proj
end

function direct_sum(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])],[id(X[1])]
    end
    Z,ix,px = direct_sum(X[1],X[2])
    for Y in X[3:end]
        Z,ix,px = ⊕((Z,ix,px),Y)
    end
    return Z,ix,px
end

function direct_sum(f::Morphism...)
    g = f[1]

    for h ∈ f[2:end]
        g = g ⊕ h
    end
    return g
end

function ×(X::T,Y::T) where {T <: Object}
    return product(X,Y)[1] 
end

function ×(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = product(T[1],X)
    m = vcat([t ∘ px[1] for t in T[2]], px[2])
    return Z, m
end

function product(X::Object...)
    if length(X) == 1
        return X[1], [id(X[1])]
    end
    Z,px = product(X[1],X[2])
    for Y in X[3:end]
        Z,px = ×((Z,px),Y)
    end
    return Z,px
end

function ∐(X::T,Y::T) where {T <: Object}
    return coproduct(T[1],X)[1]
end

function ∐(T::Tuple{S,Vector{R}},X::S1) where {S <: Object,S1 <: Object, R <: Morphism}
    Z,px = coproduct(T[1],X)
    m = vcat([px[1] ∘ t for t in T[2]], px[2])
    return Z, m
end
function coproduct(X::Object...)
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
×(X::Object...) = product(X...)[1]

"""
    ∐(X::Object...)

Return the coproduct Object and an array containing the injection morphisms.
"""
∐(X::Object...) = coproduct(X...)[1]

"""
    ⊕(X::Object...)

Return the direct sum Object and arrays containing the injection and projection
morphisms.
"""

⊕(X::Vector{<:Object}) = direct_sum(X...)[1]
⊕(X::Object...) = direct_sum(X...)[1]

⊕(X::Morphism...) = direct_sum(X...)

"""
    ⊗(X::Object...)

Return the tensor product object.
"""
⊗(X::T...) where T <: Object = tensor_product(X...)
⊗(X::Object, Y::Object) = tensor_product(X,Y)

⊗(C::Category, K::Field) = extension_of_scalars(C,K)
⊗(X::Object, K::Field) = extension_of_scalars(X,K)
⊗(f::Morphism, K::Field) = extension_of_scalars(f,K)

"""
    ^(X::Object, n::Integer)

Return the n-fold product object ```X^n```.
"""
^(X::Object,n) = n == 0 ? zero(parent(X)) : product([X for i in 1:n]...)[1]

^(X::Morphism,n) = n == 0 ? zero_morphism(zero(parent(domain(X))), zero(parent(domain(X)))) : direct_sum([X for i in 1:n]...)
"""
    ⊗(f::Morphism, g::Morphism)

Return the tensor product morphism of ```f```and ```g```.
"""
⊗(f::Morphism, g::Morphism) = tensor_product(f,g)


direct_sum(X::T) where T <: Union{Vector,Tuple} = direct_sum(X...)
product(X::T) where T <: Union{Vector,Tuple} = product(X...)
coproduct(X::T) where T <: Union{Vector,Tuple} = coproduct(X...)

product(X::Object,Y::Object) = direct_sum(X,Y)[[1,3]]
coproduct(X::Object, Y::Object) = direct_sum(X,Y)[[1,2]]



""" 

    initial_object(C::Category)

Return the initial object of `C`.
"""
function initial_object(C::Category)
    @assert is_additive(C) "No initial object known"
    return zero(C)
end

""" 

    terminal_object(C::Category)

Return the terminal object of C.
"""
function terminal_object(C::Category)
    @assert is_additive(C) "No terminal object known"
    return zero(C)
end

compose(f::T...) where T <: Morphism = reduce(compose, f)
∘(f::Morphism...) = compose(reverse(f)...)
∘(f::AbstractFunctor...) = compose(reverse(f)...)


-(f::Morphism, g::Morphism) = f + (-1)*g
-(f::Morphism) = (-1)*f

getindex(C::Category, x::Int) = simples(C)[x]
getindex(C::Category, x::Vector{Int}) = [C[i] for i ∈ x]

function getindex(C::Category, x::Int...) 
    S = simples(C)
    direct_sum([S[i] for i ∈ x])[1]
end

#-------------------------------------------------------
# Hom Spaces
#-------------------------------------------------------

dim(V::HomSpace) = length(basis(V))

End(X::Object) = Hom(X,X)


Base.iterate(H::AbstractHomSpace, state = 1) = state > int_dim(H) ? nothing : (basis(H)[state], state + 1)
getindex(H::AbstractHomSpace, k) = getindex(basis(H),k)
Base.length(H::AbstractHomSpace) = int_dim(H)
Base.eltype(H::Type{T}) where T <: AbstractHomSpace = fieldtype(H, :basis).parameters[1] 


(R::Ring)(f::Morphism) = morphism_to_scalar(R,f)
(R::QQBarField)(f::Morphism) = morphism_to_scalar(R,f)

function morphism_to_scalar(R::Ring, f::Morphism)  
    if is_zero(f)
        return zero(R)
    end
    B = basis(Hom(domain(f), codomain(f)))
    if length(B) == 0 
        return 0
    elseif length(B) == 1
        if domain(f) == codomain(f)
            return express_in_basis(f,[id(domain(f))])[1]
        else
            return express_in_basis(f,B)[1]
        end
    end

    try 
        m = matrix(f)
        if m == zero(parent(m))
            return zero(F)
        end
        b,c = is_scalar_multiple(m, matrix(id(domain(f))))
        if b 
            return c
        end
    catch end
    # m = collect(m)[m .!= 0]
    # if size(m) == (1,)
    #     return F(m[1,1])
    # end
    throw(ErrorException("Cannot convert to element of $R"))
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

*(f::Morphism, x) = x*f

#=----------------------------------------------------------
    Wrapper type fallbacks 
----------------------------------------------------------=#


*(k, f::Morphism) = morphism(domain(f), codomain(f), k*morphism(f))

+(f::T, g::T) where T <: Morphism = morphism(domain(f), codomain(g), morphism(f) + morphism(g))

compose(f::T, g::T) where T <: Morphism = morphism(domain(f), codomain(g), compose(morphism(f), morphism(g)))

function matrix(f::Morphism) 
    try 
        return matrix(morphism(f))
    catch
        error("matrix(::$(typeof(f))) not implemented")
    end
end