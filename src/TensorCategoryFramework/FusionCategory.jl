

#=----------------------------------------------------------
    Category defined by 6j-Symbols with finitely many
    simple objects 
----------------------------------------------------------=#
mutable struct SixJCategory <: Category
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
    
    function SixJCategory(F::Ring, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1,1,:])])
        C = new(F, length(mult[1,1,:]), names)
        set_tensor_product!(C,mult)
        set_spherical!(C, [F(1) for _ ∈ names])
                #C.ass = [id(⊗(X,Y,Z)) for X ∈ simples(C), Y ∈ simples(C), Z ∈ simples(C)]
                return C
    end

    function SixJCategory(F::Ring, names::Vector{String})
        C = new(F,length(names), names)
                set_spherical!(C, [F(1) for _ ∈ names])
        
        (C)
        return C
    end

    function SixJCategory()
        new()
    end

end

struct SixJObject <: Object
    parent::SixJCategory
    components::Vector{Int}
end

struct SixJMorphism <: Morphism
    domain::SixJObject
    codomain::SixJObject
    m::Vector{<:MatElem}
end



#-------------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------------

#SixJCategory(x...) = SixJCategory(x...)

Morphism(X::SixJObject, Y::SixJObject, m) = SixJMorphism(X,Y,m)

#-------------------------------------------------------------------------------
#   Setters/Getters
#-------------------------------------------------------------------------------

function set_tensor_product!(F::SixJCategory, tensor)
    F.tensor_product = tensor
    n = size(tensor,1)

    ass = Array{MatElem,4}(undef,n,n,n,n)
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        ass[i,j,k,:] = matrices(id(F[i]⊗F[j]⊗F[k]))
    end
    F.ass = ass

end

function set_braiding!(F::SixJCategory, braiding)
    F.braiding = braiding
end


set_associator!(F::SixJCategory, ass) = F.ass = ass
function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, ass::Vector{<:MatElem})
    F.ass[i,j,k,:] = ass
end

function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int, ass::MatElem)
    F.ass[i,j,k,l] = ass
end

function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int, ass::Array{T,N}) where {T,N}
    F.ass[i,j,k,l] = matrix(base_ring(F), (N > 1 ? size(ass) : (1,1))..., ass)
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
indecomposables_names(C::SixJCategory) = C.simples_names

#(::Type{Int})(x::fmpq) = Int(numerator(x))

function braiding(X::SixJObject, Y::SixJObject) 
    if is_simple(X) && is_simple(Y)
        i = findfirst(e -> e != 0, X.components)
        j = findfirst(e -> e != 0, Y.components)
        return Morphism(X⊗Y,Y⊗X, parent(X).braiding[i,j,:])
    end

    simple_objects = simples(parent(X))
    n = length(simple_objects)

    X_summands = vcat([[s for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Y_summands = vcat([[s for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)

    braid = direct_sum([braiding(x,y) for x ∈ X_summands, y ∈ Y_summands][:])

    return braid 

    distr_before = direct_sum([distribute_right(x,Y_summands) for x ∈ X_summands]) ∘ distribute_left(X_summands,Y) 
    distr_after = direct_sum([distribute_left(X_summands,Y) for y ∈ Y_summands]) ∘ distribute_right(Y_summands,X)
    
    return inv(distr_after) ∘ braid ∘ distr_before
end

associator(C::SixJCategory) = C.ass


"""
    associator(X::SixJObject, Y::SixJObject, Z::SixJObject)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
 #= @memoize Dict =# function associator(X::SixJObject, Y::SixJObject, Z::SixJObject)
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
        i = findfirst(e -> e ≠ 0, X.components)
        j = findfirst(e -> e ≠ 0, Y.components)
        k = findfirst(e -> e ≠ 0, Z.components)
        return Morphism(dom,dom, C_associator[i,j,k,:])

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

#= @memoize Dict =# function inv_associator(X::SixJObject, Y::SixJObject, Z::SixJObject)
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

        i = findfirst(e -> e ≠ 0, X.components)
        j = findfirst(e -> e ≠ 0, Y.components)
        k = findfirst(e -> e ≠ 0, Z.components)
        return inv(Morphism(dom,dom, C_associator[i,j,k,:]))

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



is_simple(X::SixJObject) = sum(X.components) == 1

==(X::SixJObject, Y::SixJObject) = parent(X) == parent(Y) && X.components == Y.components
==(f::SixJMorphism, g::SixJMorphism) = domain(f) == domain(g) && codomain(f) == codomain(g) && f.m == g.m


decompose(X::SixJObject, simpls::Vector{SixJObject} = SixJObject[]) = [(x,k) for (x,k) ∈ zip(simples(parent(X)), X.components) if k != 0]


inv(f::SixJMorphism) = SixJMorphism(codomain(f),domain(f), inv.(f.m))


id(X::SixJObject) = SixJMorphism(X,X, [one(MatrixSpace(base_ring(X),d,d)) for d ∈ X.components])


function compose(f::SixJMorphism, g::SixJMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"

    return SixJMorphism(domain(f), codomain(g), [m*n for (m,n) ∈ zip(f.m,g.m)])
end

#function vertical_direct_sum(f::SixJMorphism, g::SixJMorphism)

function +(f::SixJMorphism, g::SixJMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    SixJMorphism(domain(f), codomain(f), [m + n for (m,n) ∈ zip(f.m,g.m)])
end

function tr(f::SixJMorphism)
    # Make use of the fact that the trace is invariant under basis transformation.

    return sum([left_trace(f[i]) for i ∈ 1:parent(f).simples])

end

"""
    dual(X::SixJObject)

Return the dual object of ``X``. An error is thrown if ``X`` is not rigid.
"""
function dual(X::SixJObject)
    C = parent(X)

    # Dual of simple Object
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
            throw(ErrorException("Object not rigid."))
        end
        return SixJObject(C,[i == j[1] ? 1 : 0 for i ∈ 1:C.simples])
    end

    # Build dual from simple objects
    return direct_sum([dual(Y)^(X.components[i]) for (Y,i) ∈ zip(simples(C), 1:C.simples)])[1]
end



function coev(X::SixJObject)
    C = parent(X)
    if X == zero(C)
        return zero_morphism(one(C),X)
    end
    𝟙 = one(parent(X))
    ks = findall(e -> e > 0, X.components)
    if length(ks) == 1
        c = simple_objects_coev(C[ks[1]])
        k = X.components[ks[1]]
        m = collect(identity_matrix(base_ring(X), k))[:]
        #m = collect(matrix(coev(VectorSpaceObject(base_ring(X),k))))[:]

        return vertical_direct_sum([i * c for i ∈ m]...)
        # cod = X ⊗ dual(X)
        # n = matrices(zero_morphism(𝟙, cod))
        # n[1] = m
        # return Morphism(𝟙, cod, n)
    end

    C = parent(X)


    summands = [x^k for (x,k) ∈ decompose(X)]
    dual_summands = dual.(summands)
    d = length(summands)

    c = vertical_direct_sum([i == j ? coev(summands[i]) : zero_morphism(𝟙, summands[j]⊗dual_summands[i]) for j ∈ 1:d, i ∈ 1:d][:])
    
    distr = direct_sum([distribute_right(x,dual_summands) for x ∈ summands]) ∘ distribute_left(summands, dual(X))

    return distr ∘ c
end

function ev(X::SixJObject)
    C = parent(X)
    if X == zero(C)
        return zero_morphism(X,one(C))
    end
    𝟙 = one(parent(X))
    ks = findall(e -> e > 0, X.components)
    if length(ks) == 1
        e = simple_objects_ev(C[ks[1]])
        k = X.components[ks[1]]
        m = collect(identity_matrix(base_ring(X), k))[:]
        #m = collect(matrix(coev(VectorSpaceObject(base_ring(X),k))))[:]

        return horizontal_direct_sum([i * e for i ∈ m]...)
        # dom = dual(X) ⊗ X
        # n = matrices(zero_morphism(dom, 𝟙))
        # n[1] = m
        # return Morphism(dom, 𝟙, n)
    end

    summands = [x^k for (x,k) ∈ decompose(X)]
    dual_summands = dual.(summands)
    d = length(summands)

    e = horizontal_direct_sum(SixJMorphism[i == j ? ev(summands[i]) : zero_morphism(dual_summands[j]⊗summands[i], 𝟙)  for j ∈ 1:d, i ∈ 1:d][:])

    distr = direct_sum([distribute_right(x,summands) for x ∈ dual_summands]) ∘ distribute_left(dual_summands, X)

    return e ∘ inv(distr) 
end

function simple_objects_coev(X::SixJObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    cod = X ⊗ DX

    if sum(X.components) == 0 return zero_morphism(one(C), X) end

    return basis(Hom(one(C), cod))[1]
    #mats = [diagonal_matrix(F(1),n,m) for (n,m) ∈ zip(C.one, cod.components)]
    #return Morphism(one(C), cod, mats)
end

function simple_objects_ev(X::SixJObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    dom = DX ⊗ X

    if sum(X.components) == 0 return zero_morphism(X,one(C)) end

    #mats = [diagonal_matrix(F(1),n,m) for (n,m) ∈ zip(dom.components, C.one)]
    #unscaled_ev = Morphism(dom, one(C), mats)
    unscaled_ev = basis(Hom(dom, one(C)))[1]

    factor = F((id(X)⊗unscaled_ev)∘associator(X,DX,X)∘(coev(X)⊗id(X)))


    return inv(factor) * unscaled_ev
end

function spherical(X::SixJObject)
    C = parent(X)
    F = base_ring(C)
    sp = C.spherical
    mats = [diagonal_matrix(θ, k) for (θ,k) ∈ zip(sp, X.components)]
    return Morphism(X,X,mats)
end


*(λ,f::SixJMorphism) = SixJMorphism(domain(f), codomain(f), λ .*f.m)


function getindex(f::SixJMorphism, i)
    simple = simples(parent(domain(f)))
    dom = simple[i]^domain(f).components[i]
    cod = simple[i]^codomain(f).components[i]
    m = zero_morphism(dom,cod).m
    m[i] = f.m[i]
    return SixJMorphism(dom,cod,m)
end



getindex(X::SixJObject, i::Int64) = X.components[i]


function matrices(f::SixJMorphism)
    f.m
end


function matrix(f::SixJMorphism)
    diagonal_matrix(f.m)
end


# function (F::Field)(f::SixJMorphism)
#     if !(domain(f) == codomain(f) && is_simple(domain(f)))
#         throw(ErrorException("Cannot convert Morphism to $F"))
#     end
#     i = findfirst(e -> e == 1, domain(f).components)
#     return F(f.m[i][1,1])
# end

function dim(X::SixJObject)
    if X == zero(parent(X))
        return base_ring(X)(0)
    end

    K = base_ring(X)

    if is_simple(X)
        return K(tr(id(X)))
    end

    return sum([k*dim(x) for (x,k) ∈ decompose(X)])
end

        

#-------------------------------------------------------------------------------
#   Tensor Product
#-------------------------------------------------------------------------------



function tensor_product(X::SixJObject, Y::SixJObject)
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

    return SixJObject(C,T)
end


function tensor_product(f::SixJMorphism, g::SixJMorphism)
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
                m[k] = kronecker_product(A, identity_matrix(base_ring(C),c))

                h = h ⊕ SixJMorphism(simpl[k]^(c*d1),simpl[k]^(c*d2), m)
                
            end
        end
    end
    #dom_left = dom.components - domain(h).components
    #cod_left = cod.components - codomain(h).components
    return h #⊕ zero_morphism(SixJObject(C,dom_left), SixJObject(C,cod_left))
end



function one(C::SixJCategory) 
    if !isdefined(C, :one) 
        throw(ErrorException("There is no unit object defined"))
    end
    SixJObject(C,C.one)
end


#-------------------------------------------------------------------------------
#   Direct sum
#-------------------------------------------------------------------------------

# function direct_sum(X::SixJObject, Y::SixJObject)
#     S = SixJObject(parent(X), X.components .+ Y.components)
#     ix_mats = matrices(zero_morphism(X,S))
#     iy_mats = matrices(zero_morphism(Y,S))
#     px_mats = matrices(zero_morphism(S,X))
#     py_mats = matrices(zero_morphism(S,Y))

#     for i ∈ 1:parent(X).simples
#         (x,y) = X.components[i], Y.components[i]
#         for j ∈ 1:x 
#             ix_mats[i][j,j] = 1
#             px_mats[i][j,j] = 1
#         end
#         for j ∈ 1:y 
#             iy_mats[i][j,j+x] = 1
#             py_mats[i][j+x,j] = 1
#         end
#     end

#     ix = Morphism(X,S, ix_mats)
#     px = Morphism(S,X, px_mats)
#     iy = Morphism(Y,S, iy_mats)
#     py = Morphism(S,Y, py_mats)

#     return S,[ix,iy],[px,py]
# end

function direct_sum(X::SixJObject...)
    if length(X) == 1
        return X...,[id(X...)], [id(X...)]
    end

    S = SixJObject(parent(X[1]), vec(sum(hcat([x.components for x in X]...), dims = 2)))

    inc = [matrices(zero_morphism(x,S)) for x ∈ X]
    proj = [matrices(zero_morphism(S,x)) for x ∈ X]

    for i ∈ 1:parent(X[1]).simples
        k = [x[i] for x ∈ X]
        for j ∈ 1:length(k)
            for l ∈ 1:k[j]
                shift = sum(k[1:j-1])
                inc[j][i][l, shift + l] = 1
                proj[j][i][shift + l, l] = 1
            end
        end
    end
    inc = [Morphism(x,S,i) for (x,i) ∈ zip(X,inc)]
    proj = [Morphism(S,x,p) for (x,p) ∈ zip(X,proj)]

    return S, inc, proj
end


function ⊕(X::SixJObject...) 
    SixJObject(parent(X[1]), vec(sum(hcat([x.components for x in X]...), dims = 2)))
end

function ^(X::SixJObject, k::Int)
    SixJObject(parent(X), k.*(X.components))
end

function direct_sum(f::SixJMorphism, g::SixJMorphism)
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

    return Morphism(dom,cod, m)
end

function vertical_direct_sum(f::SixJMorphism...)
    if length(f) == 1
        return f[1]
    end
    
    #@assert length(unique!([domain.(f)...])) == 1 "Not compatible"

    ms = matrices.(f)
    m = [hcat([n[i] for n ∈ ms]...) for i ∈ 1:parent(f[1]).simples]
    return Morphism(domain(f[1]), ⊕(codomain.(f)...), m)
end

function horizontal_direct_sum(f::SixJMorphism...)
    if length(f) == 1
        return f[1]
    end
    # @assert length(unique!([codomain.(f)...])) == 1 "Not compatible"

    ms = matrices.(f)
    m = [vcat([n[i] for n ∈ ms]...) for i ∈ 1:parent(f[1]).simples]
    return Morphism(⊕(domain.(f)...), codomain(f[1]), m)
end


zero(C::SixJCategory) = SixJObject(C,[0 for i ∈ 1:C.simples])

function zero_morphism(X::SixJObject, Y::SixJObject)
    return SixJMorphism(X,Y,[zero(MatrixSpace(base_ring(X), cX, cY)) for (cX,cY) ∈ zip(X.components, Y.components)])
end

function is_isomorphic(X::SixJObject, Y::SixJObject)
    if X != Y
        return false, nothing
    else
        return true, id(X)
    end
end
#-------------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------------

function simples(C::SixJCategory)
    n = C.simples
    [SixJObject(C, [i == j ? 1 : 0 for j ∈ 1:n]) for i ∈ 1:n]
end

#-------------------------------------------------------------------------------
#   Kernel and Cokernel
#-------------------------------------------------------------------------------

function kernel(f::SixJMorphism)
    C = parent(domain(f))
    kernels = [kernel(m) for m ∈ f.m]
    
    ker = SixJObject(C,[k for (k,m) ∈ kernels])

    return ker, Morphism(ker, domain(f), [m for (_,m) ∈ kernels])
end

function cokernel(f::SixJMorphism)
    C = parent(domain(f))
    cokernels = [kernel(transpose(m)) for m ∈ f.m]
    
    coker = SixJObject(C,[k for (k,m) ∈ cokernels])

    return coker, Morphism(codomain(f),coker, [m for (_,m) ∈ cokernels])
end


function left_inverse(f::SixJMorphism)
    inverses = [left_inverse(Morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return Morphism(codomain(f), domain(f), mats)
end

function right_inverse(f::SixJMorphism)
    inverses = [right_inverse(Morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return Morphism(codomain(f), domain(f), mats)
end



#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct SixJHomSpace<: AbstractHomSpace
    X::SixJObject
    Y::SixJObject
    basis::Vector{SixJMorphism}
end

function Hom(X::SixJObject, Y::SixJObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Xi, Yi = X.components, Y.components
    F = base_ring(X)

    d = sum([x*y for (x,y) ∈ zip(Xi,Yi)])

    if d == 0 return SixJHomSpace(X,Y,SixJMorphism[]) end

    basis = [zero_morphism(X,Y).m for i ∈ 1:d]
    next = 1
    for k ∈ 1:parent(X).simples

        for i ∈ 1:Xi[k], j ∈ 1:Yi[k]
            basis[next][k][i,j] = 1
            next = next + 1
        end
    end
    basis_mors = [SixJMorphism(X,Y,m) for m ∈ basis]
    return SixJHomSpace(X,Y,basis_mors)
end

function express_in_basis(f::SixJMorphism, base::Vector{SixJMorphism})
    F = base_ring(domain(f))
    A = Array{elem_type(F),2}(undef,int_dim(Hom(domain(f),codomain(f))),0)
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

function show(io::IO, X::SixJObject)
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

function show(io::IO, f::SixJMorphism)
    print(io, """Morphism with
Domain: $(domain(f))
Codomain: $(codomain(f))
Matrices: """)
    print(io, join(["$(m)" for m ∈ f.m], ", "))
end

#-------------------------------------------------------------------------------
#   Utility
#-------------------------------------------------------------------------------

""" 

    extension_of_scalars(C::SixJCategory, K::Field)

Return the category ``C⊗K``.
"""
function extension_of_scalars(C::SixJCategory, L::Field)
    K = base_ring(C)
    if K != QQ && characteristic(K) == 0 
        if K isa AnticNumberField && L isa NfRel
            f = L
        else
            if base_field(K) == base_field(L)
                _,f = is_subfield(K,L)
            else
                f = L
            end
        end
    else
        f = L
    end

    try
        D = SixJCategory(L, C.tensor_product, simples_names(C))

        if isdefined(C, :ass)
            D.ass = [matrix(L, size(a)..., f.(collect(a))) for a ∈ C.ass]
        end
        if isdefined(C, :one)
            D.one = C.one
        end
        if isdefined(C, :spherical)
            D.spherical = f.(C.spherical)
        end
        if  isdefined(C, :braiding)
            D.braiding = [matrix(L, size(a)..., f.(collect(a))) for a ∈ C.braiding]
        end
        if isdefined(C, :twist) 
            D.twist = f.(C.twist)
        end
        return D
    catch 
        error("Extension of scalars not possible")
    end
end

""" 

    extension_of_scalars(X::SixJObject, K::Field)

Return the object ``X`` as an object of the category ``C⊗K``.
"""
function extension_of_scalars(X::SixJObject, L::Field, CL = parent(X) ⊗ L)
    SixJObject(CL, X.components)
end

""" 

    extension_of_scalars(f::SixJmorphism, K::Field)

Return the category ``C⊗K``.
"""
function extension_of_scalars(m::SixJMorphism, L::Field)
    K = base_ring(m)
    if K != QQ && characteristic(K) == 0 
        if K isa AnticNumberField && L isa NfRel
            f = L
        else
            if base_field(K) == base_field(L)
                _,f = is_subfield(K,L)
            else
                f = L
            end
        end
    else
        f = L
    end



    try

        mats = [matrix(L, size(m)..., f.(collect(m))) for m ∈ matrices(m)] 
        g = Morphism(domain(m) ⊗ L, codomain(m) ⊗ L, mats)

        return g
    catch
        error("Extension of scalars not possible")
    end
end

#=----------------------------------------------------------
    Reverse braided  
----------------------------------------------------------=#

function reverse_braiding(C::SixJCategory)
    @assert is_braided(C)

    D = deepcopy(C)
    n = length(simples(D))
    for i ∈ 1:n, j ∈ 1:n
        D.braiding[i,j,:] = inv.(C.braiding[j,i,:])
    end

    set_name!(D, "$(C.name) with reversed braiding")
    D
end