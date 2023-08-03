
abstract type SixJCategory <: Category end
abstract type SixJCategoryObject <: CategoryObject end
abstract type SixJCategoryMorphism <: CategoryMorphism end
#=----------------------------------------------------------
    Category defined by 6j-Symbols with finitely many
    simple objects 
----------------------------------------------------------=#
mutable struct FinSixJCategory <: SixJCategory
    base_ring::Ring
    simples::Int64
    simples_names::Vector{String}
    ass::Array{MatElem,4}
    braiding::Array{MatElem,3}
    tensor_product::Array{Int,3}
    spherical::Vector
    twist::Vector
    one::Vector{Int}
    name::String

    function FinSixJCategory(F::Ring, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1,1,:])])
        C = new(F, length(mult[1,1,:]), names)
        set_tensor_product!(C,mult)
        set_spherical!(C, [F(1) for _ ∈ names])
    
        #C.ass = [id(⊗(X,Y,Z)) for X ∈ simples(C), Y ∈ simples(C), Z ∈ simples(C)]
        #C.dims = [1 for i ∈ 1:length(names)]
        return C
    end

    function FinSixJCategory(F::Ring, names::Vector{String})
        C = new(F,length(names), names)
        #C.dims = [1 for i ∈ 1:length(names)]
        set_spherical!(C, [F(1) for _ ∈ names])
        (C)
        return C
    end

    function SixJCategory()
        new()
    end

end

SixJCategory(F::Ring, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1,1,:])]) = FinSixJCategory(F,mult,names)

SixJCategory(F::Ring, names::Vector{String}) = FinSixJCategory(F,names)

struct FinSixJCategoryObject <: SixJCategoryObject
    parent::SixJCategory
    components::Vector{Int}
end

SixJCategoryObject(C::FinSixJCategory, c::Vector) = FinSixJCategoryObject(C,c)

struct FinSixJCategoryMorphism <: SixJCategoryMorphism
    domain::SixJCategoryObject
    codomain::SixJCategoryObject
    m::Vector{<:MatElem}
end

SixJCategoryMorphism(domain::SixJCategoryObject,
    codomain::SixJCategoryObject,
    m::Vector{<:MatElem}) = FinSixJCategoryMorphism(domain,codomain,m)

#=----------------------------------------------------------
    Category from 6j-Symbols with countably many simple
    objects 
----------------------------------------------------------=#

mutable struct InfSixJCategory <: SixJCategory
    base_ring::Ring
    simples_names::Function
    ass::Function
    dual::Function
    braiding::Function
    tensor_product::Function
    spherical::Function
    twist::Vector
    one::Vector
    name::String

    function SixJCategory(F::Ring, mult::Function, names::Function = n::ZZRingElem -> "X$n")
        C = new(F,names)
        set_tensor_product!(C,mult)
        return C
    end

    function SixJCategory(F::Ring, mult::Function)
        C = new(F,mult)
    end

end

struct InfSixJCategoryObject <: SixJCategoryObject
    parent::SixJCategory
    components::Dict{ZZRingElem,ZZRingElem}
end

SixJCategoryObject(C::InfSixJCategory, c::Dict) = FinSixJCategoryObject(C,c)

struct InfSixJCategoryMorphism <: SixJCategoryMorphism
    domain::SixJCategoryObject
    codomain::SixJCategoryObject
    m::Dict{ZZRingElem, <:MatElem}
end

SixJCategoryMorphism(domain::InfSixJCategoryObject,
codomain::InfSixJCategoryObject,
m::Dict{ZZRingElem, <:MatElem}) = InfSixJCategoryMorphism(domain,codomain,m)

#-------------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------------

#SixJCategory(x...) = SixJCategory(x...)

Morphism(X::SixJCategoryObject, Y::SixJCategoryObject, m) = SixJCategoryMorphism(X,Y,m)

#-------------------------------------------------------------------------------
#   Setters/Getters
#-------------------------------------------------------------------------------

function set_tensor_product!(F::SixJCategory, tensor)
    F.tensor_product = tensor
    n = size(tensor,1)
    if typeof(F) == FinSixJCategory
        ass = Array{MatElem,4}(undef,n,n,n,n)
        for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
            ass[i,j,k,:] = matrices(id(F[i]⊗F[j]⊗F[k]))
        end
        F.ass = ass
    end
end

function set_braiding!(F::SixJCategory, braiding)
    F.braiding = braiding
end

set_associator!(F::SixJCategory, ass) = F.ass = ass
function set_associator!(F::FinSixJCategory, i::Int, j::Int, k::Int, ass::Vector{<:MatElem})
    F.ass[i,j,k,:] = ass
end

function set_associator!(F::FinSixJCategory, i::Int, j::Int, k::Int, l::Int, ass::MatElem)
    F.ass[i,j,k,l] = ass
end

function set_associator!(F::FinSixJCategory, i::Int, j::Int, k::Int, l::Int, ass::Array{T,N}) where {T,N}
    F.ass[i,j,k,l] = matrix(base_ring(F), (N > 1 ? size(ass) : (1,1))..., ass)
end

function set_dual!(C::InfSixJCategory, dual::Function)
    C.dual = dual
end

function set_spherical!(F::SixJCategory, sp)
    F.spherical = sp
end

function set_canonical_spherical!(C::SixJCategory)
    @assert is_fusion(C)
    set_spherical!(C, [fpdim(s)*inv(dim(s)) for s ∈ simples(C)])
end

function set_one!(F::SixJCategory, v) 
    F.one = v
end 

function set_ribbon!(F::SixJCategory, r)
    F.ribbon = r
end

function set_twist!(F::SixJCategory, t)
    F.twist = t
end

function set_name!(F::SixJCategory, name)
    F.name = name
end

function set_simples_name!(F::SixJCategory, names::Vector{String})
    F.simples_names = names
end

simples_names(C::SixJCategory) = C.simples_names

dim(X::SixJCategoryObject) = base_ring(X)(tr(id(X)))

#(::Type{Int})(x::fmpq) = Int(numerator(x))

function braiding(X::FinSixJCategoryObject, Y::FinSixJCategoryObject) 
    if is_simple(X) && is_simple(Y)
        i = findfirst(e -> e != 0, X.components)
        j = findfirst(e -> e != 0, Y.components)
        return Morphism(X⊗Y,Y⊗X, parent(X).braiding[i,j,:])
    end

    simple_objects = simples(parent(X))

    X_summands = vcat([[s for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Y_summands = vcat([[s for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)

    braid = direct_sum([braiding(x,y) for x ∈ X_summands, y ∈ Y_summands][:])

    distr_before = direct_sum([distribute_right(x,Y_summands) for x ∈ X_summands]) ∘ distr_left(X_summands,Y) 
    distr_after = direct_sum([distribute_left(y, X_summands) for y ∈ Y_summands]) ∘ distribute_right(Y,X_summands)
    
    return inv(distr_after) ∘ braid ∘ distr_before
end

associator(C::SixJCategory) = C.ass


"""
    associator(X::SixJCategoryObject, Y::SixJCategoryObject, Z::SixJCategoryObject)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
 function associator(X::SixJCategoryObject, Y::SixJCategoryObject, Z::SixJCategoryObject)
    @assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"

    C = parent(X)

    if zero(C) == X ⊗ Y ⊗ Z
        return zero_morphism(zero(C),zero(C))
    end
    F = base_ring(C)
    n = C.simples
    dom = X⊗Y⊗Z

    C_associator = C.ass

    #---------------------------------
    # associators on simple objects
    #---------------------------------
    if is_simple(X) && is_simple(Y) && is_simple(Z)
        if C isa FinSixJCategory
            i = findfirst(e -> e ≠ 0, X.components)
            j = findfirst(e -> e ≠ 0, Y.components)
            k = findfirst(e -> e ≠ 0, Z.components)
            return Morphism(dom,dom, C_associator[i,j,k,:])
        elseif C isa InfSixJCategory
            return Morphism(dom,dom, C_associator(i,j,k))
        end
    end

    #---------------------------------
    # associators for arbitrary objects
    #---------------------------------
    simple_objects = simples(parent(X))

    X_summands = vcat([[s for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Y_summands = vcat([[s for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Z_summands = vcat([[s for l ∈ 1:Z.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)

    #=-------------------------------------------------
        Distribution 
    -------------------------------------------------=#

    # Before
    distr_before = distribute_left(X_summands, Y) ⊗ id(Z)
    distr_before = (direct_sum([distribute_right(Xᵢ,Y_summands) for Xᵢ ∈ X_summands]...)⊗id(Z)) ∘ distr_before
    distr_before = distribute_left([Xᵢ⊗Yⱼ for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:], Z) ∘ distr_before
    distr_before = direct_sum([distribute_right(Xᵢ⊗Yⱼ,Z_summands) for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:]...) ∘ distr_before
    
    # After
    distr_after = id(X)⊗distribute_left(Y_summands, Z)
    distr_after = (id(X)⊗direct_sum([distribute_right(Yⱼ,Z_summands) for Yⱼ ∈ Y_summands]...)) ∘ distr_after
    distr_after = distribute_left(X_summands, Y⊗Z) ∘ distr_after
    YZ_arr = [Yⱼ⊗Zₖ for  Zₖ ∈ Z_summands, Yⱼ ∈ Y_summands][:]
    distr_after = direct_sum([distribute_right(Xᵢ, YZ_arr) for Xᵢ ∈ X_summands]) ∘ distr_after

    #-----------------------------------
    # Associator morphism
    #-----------------------------------
    m = zero_morphism(zero(C),zero(C))
    for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands
        m = m ⊕ associator(x,y,z)
    end

    return inv(distr_after) ∘ m ∘ distr_before
end

@memoize Dict function inv_associator(X::SixJCategoryObject, Y::SixJCategoryObject, Z::SixJCategoryObject)
    @assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"

    C = parent(X)

    if zero(C) == X ⊗ Y ⊗ Z
        return zero_morphism(zero(C),zero(C))
    end
    F = base_ring(C)
    n = C.simples
    dom = X⊗Y⊗Z

    C_associator = C.ass

    #---------------------------------
    # associators on simple objects
    #---------------------------------
    if is_simple(X) && is_simple(Y) && is_simple(Z)
        if C isa FinSixJCategory
            i = findfirst(e -> e ≠ 0, X.components)
            j = findfirst(e -> e ≠ 0, Y.components)
            k = findfirst(e -> e ≠ 0, Z.components)
            return inv(Morphism(dom,dom, C_associator[i,j,k,:]))
        elseif C isa InfSixJCategory
            return inv(Morphism(dom,dom, C_associator(i,j,k)))
        end
    end

    #---------------------------------
    # associators for arbitrary objects
    #---------------------------------
    simple_objects = simples(parent(X))

    X_summands = vcat([[s for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Y_summands = vcat([[s for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Z_summands = vcat([[s for l ∈ 1:Z.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)

    #=-------------------------------------------------
        Distribution 
    -------------------------------------------------=#

    # Before
    distr_before = distribute_left(X_summands, Y) ⊗ id(Z)
    distr_before = (direct_sum([distribute_right(Xᵢ,Y_summands) for Xᵢ ∈ X_summands]...)⊗id(Z)) ∘ distr_before
    distr_before = distribute_left([Xᵢ⊗Yⱼ for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:], Z) ∘ distr_before
    distr_before = direct_sum([distribute_right(Xᵢ⊗Yⱼ,Z_summands) for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:]...) ∘ distr_before
    
    # After
    distr_after = id(X)⊗distribute_left(Y_summands, Z)
    distr_after = (id(X)⊗direct_sum([distribute_right(Yⱼ,Z_summands) for Yⱼ ∈ Y_summands]...)) ∘ distr_after
    distr_after = distribute_left(X_summands, Y⊗Z) ∘ distr_after
    YZ_arr = [Yⱼ⊗Zₖ for  Zₖ ∈ Z_summands, Yⱼ ∈ Y_summands][:]
    distr_after = direct_sum([distribute_right(Xᵢ, YZ_arr) for Xᵢ ∈ X_summands]) ∘ distr_after

    #-----------------------------------
    # Associator morphism
    #-----------------------------------
    m = zero_morphism(zero(C),zero(C))
    for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands
        m = m ⊕ inv(associator(x,y,z))
    end

    return inv(distr_before) ∘ m ∘ distr_after
end

function vector_permutation(A::Vector,B::Vector)
    perm = Int[]
    for a ∈ A
        i = findall(e -> e == a, B)
        j = filter(e -> !(e ∈ perm), i)[1]
        perm = [perm; j]
    end
    return perm
end


#-------------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------------
is_semisimple(::SixJCategory) = true
is_multiring(::SixJCategory) = true

function is_multifusion(C::SixJCategory)
    try 
        dual.(simples(C))
    catch 
        return false
    end
    true
end

function is_fusion(C::SixJCategory)
    is_multifusion(C) && (sum(one(C).components) == 1)
end



is_simple(X::FinSixJCategoryObject) = sum(X.components) == 1
is_simple(X::InfSixJCategoryObject) = sum(values(X.components)) == 1

==(X::SixJCategoryObject, Y::SixJCategoryObject) = parent(X) == parent(Y) && X.components == Y.components
==(f::SixJCategoryMorphism, g::SixJCategoryMorphism) = domain(f) == domain(g) && codomain(f) == codomain(g) && f.m == g.m


decompose(X::FinSixJCategoryObject) = [(x,k) for (x,k) ∈ zip(simples(parent(X)), X.components) if k != 0]
decompose(X::InfSixJCategoryObject) = [(parent(X)[k],v) for (k,v) ∈ X.components]

inv(f::FinSixJCategoryMorphism) = SixJCategoryMorphism(codomain(f),domain(f), inv.(f.m))

inv(f::InfSixJCategoryMorphism) = Morphism(codomain(f), domain(f), Dict(k => inv(m) for (k,m) ∈ f.m))

id(X::FinSixJCategoryObject) = SixJCategoryMorphism(X,X, [one(MatrixSpace(base_ring(X),d,d)) for d ∈ X.components])

id(X::InfSixJCategoryObject) = Morphism(X,X, Dict(k => one(MatrixSpace(base_ring(X),d,d)) for (k,d) ∈ X.components))

function compose(f::FinSixJCategoryMorphism, g::FinSixJCategoryMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"

    return SixJCategoryMorphism(domain(f), codomain(g), [m*n for (m,n) ∈ zip(f.m,g.m)])
end

function compose(f::InfSixJCategoryMorphism, g::InfSixJCategoryMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"

    return Morphism(domain(f), codomain(g), Dict(k => m*n for (k,m,n) ∈ zip(keys(f.m),values(f.m), values(g.m))))
end


function +(f::FinSixJCategoryMorphism, g::FinSixJCategoryMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    SixJCategoryMorphism(domain(f), codomain(f), [m + n for (m,n) ∈ zip(f.m,g.m)])
end

function +(f::InfSixJCategoryMorphism, g::InfSixJCategoryMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    SixJCategoryMorphism(domain(f), codomain(f), merge(+, f.m, g.m))
end
"""
    dual(X::SixJCategoryObject)

Return the dual object of ``X``. An error is thrown if ``X`` is not rigid.
"""
function dual(X::FinSixJCategoryObject)
    C = parent(X)

    # Dual of simple CategoryObject
    if is_simple(X)
        # Check for rigidity
        i = findfirst(e -> e == 1, X.components)
        j = []
        for k ∈ 1:C.simples 
            if C.one[k] == 1
                j = [j; findall(e -> C.tensor_product[i,e,k] >= 1, 1:C.simples)]
            end
        end
        if length(j) != 1
            throw(ErrorException("CategoryObject not rigid."))
        end
        return SixJCategoryObject(C,[i == j[1] ? 1 : 0 for i ∈ 1:C.simples])
    end

    # Build dual from simple objects
    return direct_sum([dual(Y)^(X.components[i]) for (Y,i) ∈ zip(simples(C), 1:C.simples)])[1]
end

function dual(X::InfSixJCategoryObject)
    C = parent(C)

    if is_simple(X) 
        return simple_dual(X)
    end

    return direct_sum([dual(C[Y])^(i) for (Y,i) ∈ X.components])[1]
end

function coev(X::SixJCategoryObject)
    if X == zero(parent(X))
        return zero_morphism(one(parent(X)),X)
    end

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

function ev(X::SixJCategoryObject)
    if X == zero(parent(X))
        return zero_morphism(X,one(parent(X)))
    end
    if is_simple(X)
        return simple_objects_ev(X)
    end
    C = parent(X)
    𝟙 = one(C)

    summands = vcat([[x for _ ∈ 1:k] for (x,k) ∈ decompose(X)]...)
    dual_summands = dual.(summands)
    d = length(summands)

    e = horizontal_direct_sum(SixJCategoryMorphism[i == j ? ev(summands[i]) : zero_morphism(dual_summands[j]⊗summands[i], 𝟙)  for j ∈ 1:d, i ∈ 1:d][:])

    distr = direct_sum([distribute_right(x,summands) for x ∈ dual_summands]) ∘ distribute_left(dual_summands, X)

    return e ∘ inv(distr) 
end

function simple_objects_coev(X::SixJCategoryObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    cod = X ⊗ DX

    if sum(X.components) == 0 return zero_morphism(one(C), X) end

    mats = [diagonal_matrix(F(1),n,m) for (n,m) ∈ zip(C.one, cod.components)]

    return Morphism(one(C), cod, mats)
end

function simple_objects_ev(X::SixJCategoryObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    dom = DX ⊗ X

    if sum(X.components) == 0 return zero_morphism(X,one(C)) end

    mats = [diagonal_matrix(F(1),n,m) for (n,m) ∈ zip(dom.components, C.one)]

    unscaled_ev = Morphism(dom, one(C), mats)

    factor = F((id(X)⊗unscaled_ev)∘associator(X,DX,X)∘(coev(X)⊗id(X)))


    return inv(factor) * unscaled_ev
end

function spherical(X::SixJCategoryObject)
    C = parent(X)
    F = base_ring(C)
    sp = C.spherical
    mats = [diagonal_matrix(θ, k) for (θ,k) ∈ zip(sp, X.components)]
    return Morphism(X,X,mats)
end


*(λ,f::SixJCategoryMorphism) = SixJCategoryMorphism(domain(f), codomain(f), λ .*f.m)


function getindex(f::FinSixJCategoryMorphism, i)
    simple = simples(parent(domain(f)))
    dom = simple[i]^domain(f).components[i]
    cod = simple[i]^codomain(f).components[i]
    m = zero_morphism(dom,cod).m
    m[i] = f.m[i]
    return SixJCategoryMorphism(dom,cod,m)
end

function getindex(f::InfSixJCategoryMorphism, i)
    X = domain(f)
    Y = codomain(f)
    C = parent(X)
    if ZZ(i) ∈ keys(f.m)
        return Morphism(X[i], Y[i], f.m[i])
    else
        return zero_morphism(X[i],Y[i])
    end
end

getindex(X::FinSixJCategoryObject, i::Int64) = X.components[i]

getindex(X::InfSixJCategoryObject, i::Int64) = ZZ(i) ∈ keys(X.components) ? X.components[i] : zero(parent(X))

getindex(C::InfSixJCategory, i) = SixJCategoryObject(C,Dict(ZZ(i),ZZ(1)))

function matrices(f::FinSixJCategoryMorphism)
    f.m
end

function matrices(f::InfSixJCategoryMorphism)
    collect(values(f.m))
end


function matrix(f::FinSixJCategoryMorphism)
    diagonal_matrix(f.m)
end

function matrix(f::InfSixJCategoryMorphism)
    diagonal_matrix(collect(values(f.m)))
end

# function (F::Field)(f::SixJCategoryMorphism)
#     if !(domain(f) == codomain(f) && is_simple(domain(f)))
#         throw(ErrorException("Cannot convert Morphism to $F"))
#     end
#     i = findfirst(e -> e == 1, domain(f).components)
#     return F(f.m[i][1,1])
# end

#-------------------------------------------------------------------------------
#   Tensor Product
#-------------------------------------------------------------------------------



function tensor_product(X::FinSixJCategoryObject, Y::FinSixJCategoryObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    C = parent(X)
    n = C.simples
    T = [0 for i ∈ 1:n]

    Xc = X.components
    Yc = Y.components

    for (i,j) ∈ Base.product(1:n, 1:n)
        if (c = Xc[i]) != 0 && (d = Yc[j]) != 0
            coeffs = C.tensor_product[i,j,:]
            T = T .+ ((c*d) .* coeffs)
        end
    end

    return SixJCategoryObject(C,T)
end

function tensor_product(X::InfSixJCategoryObject, Y::InfSixJCategoryObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    C = parent(X)
    
    Z = zero(C)
    T = C.tensor_product

    for (x,k) ∈ X.components, (y,j) ∈ Y.components
        Z = Z ⊕ SixJCategoryObject(C, Dict(ZZ(z) => k*j*ZZ(l) for (z,l) ∈ T(x,y)))
    end

    return SixJCategoryObject(C,Z)
end

function tensor_product(f::FinSixJCategoryMorphism, g::FinSixJCategoryMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)
    C = parent(dom)

    h = zero_morphism(zero(C), zero(C))

    table = C.tensor_product
    simpl = simples(C)

    for i ∈ 1:C.simples, j ∈ 1:C.simples
        A = kronecker_product(f.m[i],g.m[j])
        d1,d2 = size(A)
        #if d1*d2 == 0 continue end
        for k ∈ 1:C.simples
            if (c = table[i,j,k]) > 0
                m = zero_morphism(simpl[k]^(c*d1),simpl[k]^(c*d2)).m
                m[k] = kronecker_product(identity_matrix(base_ring(C),c), A)

                h = h ⊕ SixJCategoryMorphism(simpl[k]^(c*d1),simpl[k]^(c*d2), m)
                
            end
        end
    end
    #dom_left = dom.components - domain(h).components
    #cod_left = cod.components - codomain(h).components
    return h #⊕ zero_morphism(SixJCategoryObject(C,dom_left), SixJCategoryObject(C,cod_left))
end

function tensor_product(f::InfSixJCategoryMorphism, g::InfSixJCategoryMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)
    C = parent(dom)

    h = zero_morphism(zero(C), zero(C))

    table = C.tensor_product


    for i ∈ keys(f.m), j ∈ keys(g.m)
        A = kronecker_product(f.m[i],g.m[j])
        d1,d2 = size(A)
        #if d1*d2 == 0 continue end
        Ks = keys(table(i,j))
        for k ∈ Ks
            if (c = Ks[k]) > 0
                m = zero_morphism(C[k]^(c*d1),C[k]^(c*d2)).m
                m[k] = kronecker_product(identity_matrix(base_ring(C),c), A)

                h = h ⊕ SixJCategoryMorphism(C[k]^(c*d1),C[k]^(c*d2), m)
                
            end
        end
    end
    #dom_left = dom.components - domain(h).components
    #cod_left = cod.components - codomain(h).components
    return h #⊕ zero_morphism(SixJCategoryObject(C,dom_left), SixJCategoryObject(C,cod_left))
end


function one(C::SixJCategory) 
    if !isdefined(C, :one) 
        throw(ErrorException("There is no unit object defined"))
    end
    SixJCategoryObject(C,C.one)
end

function one(C::InfSixJCategory)
    if !isdefined(C, :one) 
        throw(ErrorException("There is no unit object defined"))
    end
    SixJCategoryObject(C,Dict(ZZ(v) => ZZ(1) for v ∈ C.one))
end
#-------------------------------------------------------------------------------
#   Direct sum
#-------------------------------------------------------------------------------

function direct_sum(X::SixJCategoryObject, Y::SixJCategoryObject)
    S = SixJCategoryObject(parent(X), X.components .+ Y.components)
    ix_mats = matrices(zero_morphism(X,S))
    iy_mats = matrices(zero_morphism(Y,S))
    px_mats = matrices(zero_morphism(S,X))
    py_mats = matrices(zero_morphism(S,Y))

    for i ∈ 1:parent(X).simples
        (x,y) = X.components[i], Y.components[i]
        for j ∈ 1:x 
            ix_mats[i][j,j] = 1
            px_mats[i][j,j] = 1
        end
        for j ∈ 1:y 
            iy_mats[i][j,j+x] = 1
            py_mats[i][j+x,j] = 1
        end
    end

    ix = Morphism(X,S, ix_mats)
    px = Morphism(S,X, px_mats)
    iy = Morphism(Y,S, iy_mats)
    py = Morphism(S,Y, py_mats)

    return S,[ix,iy],[px,py]
end

function direct_sum(X::InfSixJCategoryObject, Y::InfSixJCategoryObject)
    S = SixJCategoryObject(parent(X), sort(merge(+, X.components, Y.components)))
    ix_mats = matrices(zero_morphism(X,S))
    iy_mats = matrices(zero_morphism(Y,S))
    px_mats = matrices(zero_morphism(S,X))
    py_mats = matrices(zero_morphism(S,Y))

    for i ∈ keys(S.components)
        x = X.components[i]
        y = Y.components[i]
        for j ∈ 1:x 
            ix_mats[i][j,j] = 1
            px_mats[i][j,j] = 1
        end
        for j ∈ 1:y 
            iy_mats[i][j,j+x] = 1
            py_mats[i][j+x,j] = 1
        end
    end

    ix = Morphism(X,S, ix_mats)
    px = Morphism(S,X, px_mats)
    iy = Morphism(Y,S, iy_mats)
    py = Morphism(S,Y, py_mats)

    return S,[ix,iy],[px,py]
end


function direct_sum(f::SixJCategoryMorphism, g::SixJCategoryMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    F = base_ring(dom)
    m = zero_morphism(dom,cod).m
    for i ∈ 1:parent(dom).simples
        mf,nf = size(f.m[i])
        mg,ng = size(g.m[i])
        z1 = zero(MatrixSpace(F,mf,ng))
        z2 = zero(MatrixSpace(F,mg,nf))
        m[i] = [f.m[i] z1; z2 g.m[i]]
    end

    return SixJCategoryMorphism(dom,cod, m)
end


zero(C::SixJCategory) = SixJCategoryObject(C,[0 for i ∈ 1:C.simples])

function zero_morphism(X::SixJCategoryObject, Y::SixJCategoryObject)
    return SixJCategoryMorphism(X,Y,[zero(MatrixSpace(base_ring(X), cX, cY)) for (cX,cY) ∈ zip(X.components, Y.components)])
end

function is_isomorphic(X::SixJCategoryObject, Y::SixJCategoryObject)
    if X != Y
        return false, nothing
    else
        return true, id(X)
    end
end
#-------------------------------------------------------------------------------
#   Simple CategoryObjects
#-------------------------------------------------------------------------------

function simples(C::SixJCategory)
    n = C.simples
    [SixJCategoryObject(C, [i == j ? 1 : 0 for j ∈ 1:n]) for i ∈ 1:n]
end

#-------------------------------------------------------------------------------
#   Kernel and Cokernel
#-------------------------------------------------------------------------------

function kernel(f::SixJCategoryMorphism)
    C = parent(domain(f))
    kernels = [kernel(Morphism(m)) for m ∈ f.m]
    mats = [matrix(m) for (_,m) ∈ kernels]
    ker = SixJCategoryObject(C,[int_dim(k) for (k,m) ∈ kernels])

    return ker, Morphism(ker, domain(f), mats)
end

function cokernel(f::SixJCategoryMorphism)
    C = parent(domain(f))
    cokernels = [cokernel(Morphism(m)) for m ∈ f.m]
    mats = [matrix(m) for (_,m) ∈ cokernels]
    coker = SixJCategoryObject(C,[int_dim(k) for (k,m) ∈ cokernels])

    return coker, Morphism(codomain(f),coker, mats)
end


function left_inverse(f::SixJCategoryMorphism)
    inverses = [left_inverse(Morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return Morphism(codomain(f), domain(f), mats)
end

function right_inverse(f::SixJCategoryMorphism)
    inverses = [right_inverse(Morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return Morphism(codomain(f), domain(f), mats)
end



#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct SixJCategoryHomSpace<: AbstractCategoryHomSpace
    X::SixJCategoryObject
    Y::SixJCategoryObject
    basis::Vector{SixJCategoryMorphism}
    parent::VectorSpaces
end

function Hom(X::SixJCategoryObject, Y::SixJCategoryObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Xi, Yi = X.components, Y.components
    F = base_ring(X)

    d = sum([x*y for (x,y) ∈ zip(Xi,Yi)])

    if d == 0 return SixJCategoryHomSpace(X,Y,SixJCategoryMorphism[], VectorSpaces(F)) end

    basis = [zero_morphism(X,Y).m for i ∈ 1:d]
    next = 1
    for k ∈ 1:parent(X).simples

        for i ∈ 1:Xi[k], j ∈ 1:Yi[k]
            basis[next][k][i,j] = 1
            next = next + 1
        end
    end
    basis_mors = [SixJCategoryMorphism(X,Y,m) for m ∈ basis]
    return SixJCategoryHomSpace(X,Y,basis_mors, VectorSpaces(F))
end

function express_in_basis(f::SixJCategoryMorphism, base::Vector{SixJCategoryMorphism})
    F = base_ring(domain(f))
    A = Array{elem_type(F),2}(undef,length(base),0)
    b = []
    for g ∈ base
        y = []
        for m ∈ g.m
            y = [y; [x for x ∈ m][:]]
        end
        A = [A y]
    end
    for m ∈ f.m
        b = [b; [x for x ∈ m][:]]
    end

    return [i for  i ∈ solve_left(transpose(matrix(F,A)), MatrixSpace(F,1,length(b))(F.(b)))][:]
end


#-------------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------------

function show(io::IO, C::SixJCategory)
    if isdefined(C,:name)
        print(io, "$(C.name)")
    else
      print(io, "Fusion Category with $(C.simples) simple objects")
    end
end

function show(io::IO, X::SixJCategoryObject)
    coeffs = X.components

    if sum(coeffs) == 0
        print(io,"0")
        return
    end

    strings = parent(X).simples_names
    non_zero_coeffs = coeffs[coeffs .> 0]
    non_zero_strings = strings[coeffs .> 0]

    disp = non_zero_coeffs[1] == 1 ? "$(non_zero_strings[1])" : "$(non_zero_coeffs[1])⋅$(non_zero_strings[1])"

    for (Y,d) ∈ zip(non_zero_strings[2:end], non_zero_coeffs[2:end])
        disp = d == 1 ? disp*" ⊕ $Y" : disp*" ⊕ $(d)⋅$Y"
    end
    print(io,disp)
end

function show(io::IO, f::SixJCategoryMorphism)
    print(io, """Morphism with
Domain: $(domain(f))
Codomain: $(codomain(f))
Matrices: """)
    print(io, join(["$(m)" for m ∈ f.m], ", "))
end

#-------------------------------------------------------------------------------
#   Utility
#-------------------------------------------------------------------------------
