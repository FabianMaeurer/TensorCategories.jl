

#=----------------------------------------------------------
    Category defined by 6j-Symbols with finitely many
    simple objects 
----------------------------------------------------------=#
@attributes mutable struct SixJCategory <: Category
    base_ring::Ring
    simples::Int64
    simples_names::Vector{String}
    ass::Array{MatElem,4}
    braiding::Array{MatElem,3}
    tensor_product::Array{Int,3}
    pivotal::Vector
    twist::Vector
    one::Vector{Int}
    name::String
    

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

==(X::SixJObject, Y::SixJObject) = base_ring(X) == base_ring(Y) && X.components == Y.components

function ==(C::SixJCategory, D::SixJCategory)
    base_ring(C) ≠ base_ring(D) && return false 
    multiplication_table(C) ≠ multiplication_table(D) && return false 
    true 
end

#-------------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------------

#six_j_category(x...) = six_j_category(x...)

morphism(X::SixJObject, Y::SixJObject, m) = SixJMorphism(X,Y,m)


function six_j_category(F::Ring, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1,1,:])])
    C = SixJCategory()
    C.base_ring = F
    C.simples = length(mult[1,1,:])
    C.simples_names = names
    set_tensor_product!(C,mult)
    set_pivotal!(C, [F(1) for _ ∈ names])
    #C.ass = [id(⊗(X,Y,Z)) for X ∈ simples(C), Y ∈ simples(C), Z ∈ simples(C)]
    return C
end

function six_j_category(F::Ring, names::Vector{String})
    C = SixJCategory()
    C.base_ring = F
    C.simples = length(names)
    C.simples_names = names
    set_pivotal!(C, [F(1) for _ ∈ names])
    
    return C
end
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

function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int, m::Int, n::Int, v::RingElem) 
    F.ass[i,j,k,l][m,n] = v
end

function set_pivotal!(F::SixJCategory, sp)
    F.pivotal = sp
end

function set_spherical!(F::SixJCategory, sp) 
    F.pivotal = sp
    @req is_spherical(F) "Not a spherical structure"

end

# function is_spherical(F::SixJCategory, sp) 
#     all([s*ev(dual(x))∘coev(x) == inv(s)*ev(dual(x))∘coev(x) for (s,x) ∈ zip(sp,simples(F))]) 
# end

function is_spherical(F::SixJCategory)
    get_attribute!(F, :is_spherical) do
        if isdefined(F, :pivotal) 
            return all([dim(x) == dim(dual(x)) for x in simples(F)])
        end
        false
    end
    
end

function set_canonical_spherical!(C::SixJCategory)
    @assert is_fusion(C)
    set_pivotal!(C, base_ring(C).([1 for _ ∈ 1:C.simples]))
    set_pivotal!(C, [fpdim(s)*inv(dim(s)) for s ∈ simples(C)])
    end

function set_one!(F::SixJCategory, v::Vector) 
    F.one = v
end 

function set_one!(F::SixJCategory, i::Int)
    F.one = [k == i for k ∈ 1:F.simples]
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

#(::Type{Int})(x::QQFieldElem) = Int(numerator(x))

function braiding(X::SixJObject, Y::SixJObject) 
    C = parent(X)
    if is_simple(X) && is_simple(Y)
        i = findfirst(e -> e != 0, X.components)
        j = findfirst(e -> e != 0, Y.components)

        if ! all(isassigned(C.braiding, i,j,k) for k ∈ 1:C.simples)
            r_symbol = get_attribute(C, :r_symbol)
            C.braiding[i,j,:] = [r_symbol(i,j,k) for k ∈ 1:C.simples]
        end
        return morphism(X⊗Y,Y⊗X, C.braiding[i,j,:])
    end

    simple_objects = simples(C)
    n = length(simple_objects)

    X_summands = vcat([[s for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Y_summands = vcat([[s for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)

    braid = direct_sum([braiding(x,y) for y ∈ Y_summands, x ∈ X_summands][:])

    distr_before = compose(
        distribute_left(X_summands,Y), 
        direct_sum([distribute_right(x,Y_summands) for x ∈ X_summands]) 
    )
    distr_after = compose( 
        distribute_right(Y,X_summands),
        direct_sum([distribute_left(Y_summands, x) for x ∈ X_summands])
    )

    return inv(distr_after) ∘ braid ∘ distr_before
end

associator(C::SixJCategory) = C.ass


"""
    associator(X::SixJObject, Y::SixJObject, Z::SixJObject)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
function associator(X::SixJObject, Y::SixJObject, Z::SixJObject)
    #@assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"

    C = parent(X)

    if zero(C) == X ⊗ Y ⊗ Z
        return zero_morphism(zero(C),zero(C))
    end
    F = base_ring(C)
    n = C.simples
    dom = X⊗Y⊗Z

    C_associator = C.ass

    if one(C) ∈ [X,Y,Z]
        return id(X ⊗ Y ⊗ Z)
    end

    #---------------------------------
    # associators on simple objects
    #---------------------------------
    if is_simple(X) && is_simple(Y) && is_simple(Z)
        i = findfirst(e -> e ≠ 0, X.components)
        j = findfirst(e -> e ≠ 0, Y.components)
        k = findfirst(e -> e ≠ 0, Z.components)

        
        return morphism(dom,dom, [six_j_symbol(C,i,j,k,l) for l ∈ 1:n])

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
    # distr_before = distribute_left(X_summands, Y) ⊗ id(Z)
    # distr_before = (direct_sum([distribute_right(Xᵢ,Y_summands) for Xᵢ ∈ X_summands]...)⊗id(Z)) ∘ distr_before
    # distr_before = distribute_left([Xᵢ⊗Yⱼ for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:], Z) ∘ distr_before
    # distr_before = direct_sum([distribute_right(Xᵢ⊗Yⱼ,Z_summands) for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:]...) ∘ distr_before

    # # After
    # distr_after = id(X)⊗distribute_left(Y_summands, Z)
    # distr_after = (id(X)⊗direct_sum([distribute_right(Yⱼ,Z_summands) for Yⱼ ∈ Y_summands]...)) ∘ distr_after
    # distr_after = distribute_left(X_summands, Y⊗Z) ∘ distr_after
    # YZ_arr = [Yⱼ⊗Zₖ for  Zₖ ∈ Z_summands, Yⱼ ∈ Y_summands][:]
    # distr_after = direct_sum([distribute_right(Xᵢ, YZ_arr) for Xᵢ ∈ X_summands]) ∘ distr_after


    #-----------------------------------
    # Associator morphism
    #-----------------------------------
    m = zero_morphism(zero(C),zero(C))
    # for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands
    #     m = m ⊕ associator(x,y,z)
    # end
    m = direct_sum([associator(x,y,z) for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands][:])

    # if length(X_summands) == 1 || length(Y_summands) == 1 || length(Z_summands) == 1
    #     return m 
    # end
    _,ix,px = direct_sum(X_summands)
    _,iy,py = direct_sum(Y_summands)
    _,iz,pz = direct_sum(Z_summands)

    distr_before = vertical_direct_sum([(f⊗g)⊗h for f ∈ px, g ∈ py, h ∈ pz][:])
    #distr_after = vertical_direct_sum([f⊗(g⊗h) for f ∈ px, g ∈ py, h ∈ pz][:])

    distr_after = horizontal_direct_sum([f⊗(g⊗h) for f ∈ ix, g ∈ iy, h ∈ iz][:])
    
    return compose(distr_before, m , distr_after)
end

function inv_associator(X::SixJObject, Y::SixJObject, Z::SixJObject)
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

       

        return morphism(dom,dom, [inv(six_j_symbol(C,i,j,k,l)) for l ∈ 1:n])
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
    # distr_before = distribute_left(X_summands, Y) ⊗ id(Z)
    # distr_before = (direct_sum([distribute_right(Xᵢ,Y_summands) for Xᵢ ∈ X_summands]...)⊗id(Z)) ∘ distr_before
    # distr_before = distribute_left([Xᵢ⊗Yⱼ for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:], Z) ∘ distr_before
    # distr_before = direct_sum([distribute_right(Xᵢ⊗Yⱼ,Z_summands) for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:]...) ∘ distr_before
    
    # # After
    # distr_after = id(X)⊗distribute_left(Y_summands, Z)
    # distr_after = (id(X)⊗direct_sum([distribute_right(Yⱼ,Z_summands) for Yⱼ ∈ Y_summands]...)) ∘ distr_after
    # distr_after = distribute_left(X_summands, Y⊗Z) ∘ distr_after
    # YZ_arr = [Yⱼ⊗Zₖ for  Zₖ ∈ Z_summands, Yⱼ ∈ Y_summands][:]
    # distr_after = direct_sum([distribute_right(Xᵢ, YZ_arr) for Xᵢ ∈ X_summands]) ∘ distr_after

    _,ix,px = direct_sum(X_summands)
    _,iy,py = direct_sum(Y_summands)
    _,iz,pz = direct_sum(Z_summands)

    distr_before= horizontal_direct_sum([(f⊗g)⊗h for f ∈ ix, g ∈ iy, h ∈ iz][:])
    distr_after = vertical_direct_sum([f⊗(g⊗h) for f ∈ px, g ∈ py, h ∈ pz][:])

    #distr_after = horizontal_direct_sum([f⊗(g⊗h) for f ∈ ix, g ∈ iy, h ∈ iz][:])
    #-----------------------------------
    # Associator morphism
    #-----------------------------------
    m = zero_morphism(zero(C),zero(C))
    # for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands
    #     m = m ⊕ inv(associator(x,y,z))
    # end

    m = direct_sum([inv_associator(x,y,z) for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands][:])

    return compose(distr_after, m, distr_before)
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

function six_j_symbol(C::SixJCategory, i::Int, j::Int, k::Int, l::Int)
    if ! isassigned(C.ass, i,j,k,l)
        six_j_symbol = get_attribute(C, :six_j_symbol)
        C.ass[i,j,k,l] = six_j_symbol(i,j,k,l) 
    end
    return C.ass[i,j,k,l]
end

function r_symbol(C::SixJCategory, i::Int, j::Int, k::Int)
    if ! isassigned(C.braiding, i,j,k)
        r_symbol = get_attribute(C, :r_symbol)
        C.braiding[i,j,k,l] = r_symbol(i,j,k) 
    end
    return C.braiding[i,j,k]
end


#-------------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------------
is_semisimple(::SixJCategory) = true
is_multiring(::SixJCategory) = true
is_braided(C::SixJCategory) = isdefined(C, :braiding)

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

# ==(X::SixJObject, Y::SixJObject) = parent(X) == parent(Y) && X.components == Y.components
==(f::SixJMorphism, g::SixJMorphism) = domain(f) == domain(g) && codomain(f) == codomain(g) && f.m == g.m


decompose(X::SixJObject, simpls::Vector{SixJObject} = SixJObject[]) = [(x,k) for (x,k) ∈ zip(simples(parent(X)), X.components) if k != 0]


inv(f::SixJMorphism) = SixJMorphism(codomain(f),domain(f), inv.(f.m))


id(X::SixJObject) = SixJMorphism(X,X, [one(matrix_space(base_ring(X),d,d)) for d ∈ X.components])


function compose(f::SixJMorphism...)
    if length(f) == 1 
        return f[1]
    end
    #@show [(codomain(f[i]), domain(f[i+1])) for i ∈ 1:length(f)-1]
    @assert all([codomain(f[i]) == domain(f[i+1]) for i ∈ 1:length(f)-1]) "Morphisms not compatible"

    return SixJMorphism(domain(f[1]), codomain(f[end]), [*(m...) for m ∈ zip(matrices.(f)...)])
end

function vertical_direct_sum(f::Vector{SixJMorphism})
    @assert all(domain(g) == domain(f[1]) for g ∈ f[2:end])

    C = parent(f[1])
    cod = SixJObject(C, reduce(.+, [X.components for X ∈ codomain.(f)]))

    mats = [hcat([g.m[i] for g ∈ f]...) for i ∈ 1:C.simples]

    morphism(domain(f[1]), cod, mats)
end

function horizontal_direct_sum(f::Vector{SixJMorphism})
    @assert all(codomain(g) == codomain(f[1]) for g ∈ f[2:end])

    C = parent(f[1])
    dom = SixJObject(C, reduce(.+, [X.components for X ∈ domain.(f)]))

    mats = [vcat([g.m[i] for g ∈ f]...) for i ∈ 1:C.simples]

    morphism(dom, codomain(f[1]), mats)
end

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

        return vertical_direct_sum([i * c for i ∈ m])
        # cod = X ⊗ dual(X)
        # n = matrices(zero_morphism(𝟙, cod))
        # n[1] = m
        # return morphism(𝟙, cod, n)
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

        return horizontal_direct_sum([i * e for i ∈ m])
        # dom = dual(X) ⊗ X
        # n = matrices(zero_morphism(dom, 𝟙))
        # n[1] = m
        # return morphism(dom, 𝟙, n)
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
    #return morphism(one(C), cod, mats)
end

function simple_objects_ev(X::SixJObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    dom = DX ⊗ X

    if sum(X.components) == 0 return zero_morphism(X,one(C)) end

    #mats = [diagonal_matrix(F(1),n,m) for (n,m) ∈ zip(dom.components, C.one)]
    #unscaled_ev = morphism(dom, one(C), mats)
    unscaled_ev = basis(Hom(dom, one(C)))[1]

    factor = F((id(X)⊗unscaled_ev)∘associator(X,DX,X)∘(coev(X)⊗id(X)))


    return inv(factor) * unscaled_ev
end

function spherical(X::SixJObject)
    @req is_spherical(parent(X)) "Not spherical"
    pivotal(X)
end

function is_pivotal(C::SixJCategory) 
    all([pivotal(x)⊗pivotal(y) == double_dual_monoidal_structure(x,y) ∘ pivotal(x ⊗ y) for x ∈ simples(C), y ∈ simples(C)])
end 

function pivotal(X::SixJObject)
    C = parent(X)
    F = base_ring(C)
    sp = C.pivotal
    mats = [diagonal_matrix(θ, k) for (θ,k) ∈ zip(sp, X.components)]
    return morphism(X,X,mats)
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
    diagonal_matrix(f.m...)
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
    #@assert parent(X) == parent(Y) "Mismatching parents"
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

#     ix = morphism(X,S, ix_mats)
#     px = morphism(S,X, px_mats)
#     iy = morphism(Y,S, iy_mats)
#     py = morphism(S,Y, py_mats)

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
    inc = [morphism(x,S,i) for (x,i) ∈ zip(X,inc)]
    proj = [morphism(S,x,p) for (x,p) ∈ zip(X,proj)]

    return S, inc, proj
end


function ⊕(X::SixJObject...) 
    SixJObject(parent(X[1]), vec(sum(hcat([x.components for x in X]...), dims = 2)))
end

function ^(X::SixJObject, k::Int)
    SixJObject(parent(X), k.*(X.components))
end

function direct_sum(f::SixJMorphism...)
    dom = ⊕(domain.(f)...)
    cod = ⊕(codomain.(f)...)
    F = base_ring(dom)

    mats = [diagonal_matrix([g.m[i] for g ∈ f]) for i ∈ 1:parent(dom).simples]

    return morphism(dom,cod, mats)
end

# function vertical_direct_sum(f::Vector{SixJMorphism})
#     if length(f) == 1
#         return f[1]
#     end
    
#     #@assert length(unique!([domain.(f)...])) == 1 "Not compatible"

#     ms = matrices.(f)
#     m = [hcat([n[i] for n ∈ ms]...) for i ∈ 1:parent(f[1]).simples]
#     return morphism(domain(f[1]), ⊕(codomain.(f)...), m)
# end

# function horizontal_direct_sum(f::Vector{SixJMorphism})
#     if length(f) == 1
#         return f[1]
#     end
#     # @assert length(unique!([codomain.(f)...])) == 1 "Not compatible"

#     ms = matrices.(f)
#     m = [vcat([n[i] for n ∈ ms]...) for i ∈ 1:parent(f[1]).simples]
#     return morphism(⊕(domain.(f)...), codomain(f[1]), m)
# end


zero(C::SixJCategory) = SixJObject(C,[0 for i ∈ 1:C.simples])

function zero_morphism(X::SixJObject, Y::SixJObject)
    return SixJMorphism(X,Y,[zero(matrix_space(base_ring(X), cX, cY)) for (cX,cY) ∈ zip(X.components, Y.components)])
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


function sort_simples!(C::SixJCategory, order::Vector{Int})
    C.tensor_product = [C.tensor_product[i,j,k] for i ∈ order, j ∈ order, k ∈ order]
    if has_attribute(C, :multiplication_table)
        set_attribute!(C, :multiplication_table, C.tensor_product)
    end
    n = C.simples
   

    C.ass = C.ass[order,order,order,order]
    # for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n, l ∈ 1:n 
    #     C.ass[order[[i,j,k,l]]...] = ass[i,j,k,l]
    # end
    
    C.simples_names = C.simples_names[order]

    isdefined(C, :one) && (C.one = C.one[order])
    isdefined(C, :pivotal) && (C.pivotal = C.pivotal[order])
    isdefined(C, :braiding) && (C.braiding = [C.braiding[i,j,k] for i ∈ order, j ∈ order, k ∈ order])
    isdefined(C, :twist) && (C.twist = C.twist[order])
    return C
end

#-------------------------------------------------------------------------------
#   Kernel and Cokernel
#-------------------------------------------------------------------------------

function kernel(f::SixJMorphism)
    C = parent(domain(f))
    kernels = [kernel(m, side = :left) for m ∈ f.m]
    
    ker = SixJObject(C,[number_of_rows(k) for k ∈ kernels])

    return ker, morphism(ker, domain(f), [m for m ∈ kernels])
end

function cokernel(f::SixJMorphism)
    C = parent(domain(f))

    cokernels = [kernel(m, side = :right) for m ∈ f.m]
    
    coker = SixJObject(C,[number_of_columns(k) for k ∈ cokernels])

    return coker, morphism(codomain(f),coker, [m for m ∈ cokernels])
end


function left_inverse(f::SixJMorphism)
     inverses = [left_inverse(morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return morphism(codomain(f), domain(f), mats)
end

function right_inverse(f::SixJMorphism)
    inverses = [right_inverse(morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return morphism(codomain(f), domain(f), mats)
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
    #@assert parent(X) == parent(Y) "Mismatching parents"
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

# function express_in_basis(f::SixJMorphism, base::Vector{SixJMorphism})
#     F = base_ring(domain(f))
#     A = Array{elem_type(F),2}(undef,length(base),0)
#     b = []
#     for g ∈ base
#         y = []
#         for m ∈ g.m
#             y = [y; [x for x ∈ m][:]]
#         end
#         A = [A y]
#     end
#     for m ∈ f.m
#         b = [b; [x for x ∈ m][:]]
#     end

#     return [i for  i ∈ solve(transpose(matrix(F,A)), matrix_space(F,1,length(b))(F.(b)))][:]
# end


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
function extension_of_scalars(C::SixJCategory, L::Ring; embedding = embedding(base_ring(C), L))


    try
        D = six_j_category(L, C.tensor_product, simples_names(C))

        set_name!(D, C.name)
        
        if isdefined(C, :ass)
            D.ass = [matrix(L, size(a)..., embedding.(collect(a))) for a ∈ C.ass]
        end
        if isdefined(C, :one)
            D.one = C.one
        end
        if isdefined(C, :pivotal)
            D.pivotal = embedding.(C.pivotal)
        end
        if  isdefined(C, :braiding)
            D.braiding = [matrix(L, size(a)..., embedding.(collect(a))) for a ∈ C.braiding]
        end
        if isdefined(C, :twist) 
            D.twist = f.(C.twist)
        end
        set_name!(D, C.name)

        return D
    catch 
        error("Extension of scalars not possible")
    end
end

function extension_of_scalars(C::SixJCategory, K::QQBarField)
    e = complex_embeddings(base_ring(C))[1]
    if base_ring(C) == QQ 
        to_qqbar = QQBarField()
    else
        to_qqbar = x -> guess(QQBarField(), e(x,1024), maximum([1,degree(x)]))
    end   
    extension_of_scalars(C,K, embedding = to_qqbar)
end

""" 

    extension_of_scalars(X::SixJObject, K::Field)

Return the object ``X`` as an object of the category ``C⊗K``.
"""
function extension_of_scalars(X::SixJObject, L::Ring, CL = parent(X) ⊗ L;embedding = nothing)
    SixJObject(CL, X.components)
end

""" 

    extension_of_scalars(f::SixJMorphism, K::Field)

Return the category ``C⊗K``.
"""
function extension_of_scalars(m::SixJMorphism, L::Ring, CL = parent(m) ⊗ L; embedding = embedding(base_ring(m), L))
    try

        mats = [matrix(L, size(m)..., embedding.(collect(m))) for m ∈ matrices(m)] 
        g = morphism(extension_of_scalars(domain(m), L, CL, embedding = embedding),
                    extension_of_scalars(codomain(m), L, CL, embedding = embedding), 
                    mats)

        return g
    catch
        error("Extension of scalars not possible")
    end
end

function restriction_of_scalars(C::SixJCategory, K::Ring)
    b,f = is_subfield(K,base_ring(C))

    !b && error("Restriction not possible")

    D = six_j_category(K, multiplication_table(C), simples_names(C))

    D.ass = [matrix(K, size(m)..., [preimage(f, a) for a ∈ m]) for m ∈ C.ass]

    isdefined(C, :spherical) && set_pivotal!(D, [preimage(f, a) for a ∈ C.spherical])

    try 
        D.braiding = [matrix(K, size(m)..., [preimage(f, a) for a ∈ m]) for m ∈ C.braiding]
    catch 
    end

    D.one = C.one 
    D.name = C.name 

    D 
end

function simplify(C::SixJCategory)
    K = base_ring(C)
end

#=----------------------------------------------------------
    Endofunctors    
----------------------------------------------------------=#

function autoequivalences(C::SixJCategory)
    if is_tambara_yamagami(C)
        return tambara_yamagami_tensor_autoequivalences(C)
    end
    error("Not implemented")
end

#=----------------------------------------------------------
    Reverse braided  
----------------------------------------------------------=#

function reverse_braiding(C::SixJCategory)
    @assert is_braided(C)

    D = six_j_category(
        base_ring(C),
        multiplication_table(C),
        simples_names(C)
    )

    set_associator!(D, associator(C))

    isdefined(C, :pivotal) && set_pivotal!(D, C.pivotal)
    isdefined(C, :one) && set_one!(D, C.one)

    n = length(simples(D))
    set_braiding!(D, Array{MatElem,3}(undef, n, n, n))
    for i ∈ 1:n, j ∈ 1:n
        D.braiding[i,j,:] = inv.(C.braiding[j,i,:])
    end

    set_name!(D, "$(C.name) with reversed braiding")
    D
end

function trivial_fusion_category(K::Field)
    C = six_j_category(K, [1 for _ ∈ 1:1, _ ∈ 1:1, _ ∈ 1:1], ["1"])

    C.spherical = [1]
    C.one = [1]
    C.name = "Trivial fusion category over $K"

    return C
end

#=----------------------------------------------------------
    Reversed monoidal structure 
----------------------------------------------------------=#

function reversed_monoidal_category(C::SixJCategory)

    D = six_j_category(base_ring(C), simples_names(C))
    n = length(simples(D))

    set_tensor_product!(D, [C.tensor_product[j,i,k] for i in 1:n, j in 1:n, k ∈ 1:n])

    set_associator!(D, Array{MatElem,4}(undef, n, n, n, n))
    set_attribute!(D, :six_j_symbol, (i,j,k,l) -> inv(six_j_symbol(C, k,j,i,l)))

    if isdefined(C, :pivotal)
        set_pivotal!(D, inv.(C.pivotal))
    end
    if isdefined(C, :one)
        set_one!(D, C.one)
    end
    set_name!(D, "$(C.name) with reversed monoidal structure")
    D
end