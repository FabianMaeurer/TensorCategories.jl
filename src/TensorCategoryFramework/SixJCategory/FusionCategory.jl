

#=----------------------------------------------------------
    Category defined by 6j-Symbols with finitely many
    simple objects 
----------------------------------------------------------=#
@attributes mutable struct SixJCategory <: Category
    base_ring::Ring
    rank::Int64
    simples_names::Vector{String}
    ass::Array{MatElem,4}
    braiding::Array{MatElem,3}
    tensor_product::Array{Int,3}
    pivotal::Vector
    twist::Vector
    one::Vector{Int}
    name::String
    embedding::AbsSimpleNumFieldEmbedding # optional: With this you can fix the category as a subcategory of a complex category

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

# Component vectors are coordinates relative to a particular skeletal category.
==(X::SixJObject, Y::SixJObject) =
    parent(X) === parent(Y) && X.components == Y.components

Base.hash(X::SixJObject, h::UInt) =
    hash((objectid(parent(X)), X.components), h)

function _check_sixj_parents(X::SixJObject...)
    isempty(X) && throw(ArgumentError(
        "at least one object is needed to determine the parent"))
    all(Y -> parent(Y) === parent(X[1]), X) || throw(ArgumentError(
        "skeletal objects must share one parent instance; transport objects explicitly between categories"))
    return nothing
end

# Categories carry mutable structural data. Equality means the same parent,
# not coincident fusion rules/F-symbols (which can have different braidings).
==(C::SixJCategory,D::SixJCategory) = C === D
Base.hash(C::SixJCategory,h::UInt) = hash(objectid(C),h)

object_type(::SixJCategory) = SixJObject
#-------------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------------

#six_j_category(x...) = six_j_category(x...)

function morphism(X::SixJObject, Y::SixJObject, m)
    _check_sixj_parents(X, Y)
    C = parent(X)
    n = rank(C)
    all(Z -> length(Z.components) == n && all(>=(0),Z.components), (X,Y)) ||
        throw(ArgumentError("invalid skeletal object multiplicities"))
    m isa AbstractVector && length(m) == n ||
        throw(ArgumentError("one matrix block per simple object is required"))
    # In the row-vector convention, block i has one row per copy in the
    # source and one column per copy in the target.
    for i in 1:n
        m[i] isa MatElem ||
            throw(ArgumentError("morphism blocks must be matrices"))
        size(m[i]) == (X.components[i],Y.components[i]) ||
            throw(ArgumentError("block $i has the wrong dimensions for its endpoints"))
        base_ring(m[i]) === base_ring(C) ||
            throw(ArgumentError("morphism blocks must use the category's coefficient field"))
    end
    SixJMorphism(X, Y, m)
end

@doc raw""" 

    six_j_category(F::Ring, mult::Array{Int,3}, [names::Vector{String}])
    six_j_category(F::Ring, names::Vector{String})

Initialize a fusion category. Associativity constraints are all set to 1, i.e. are most likely false. 
"""
function six_j_category(F::Ring, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1,1,:])])
    C = SixJCategory()
    C.base_ring = F
    C.rank = length(mult[1,1,:])
    C.simples_names = names
    set_tensor_product!(C,mult)
    set_pivotal!(C, [F(1) for _ ∈ names])
    #C.ass = [id(⊗(X,Y,Z)) for X ∈ simples(C), Y ∈ simples(C), Z ∈ simples(C)]
    return C
end

function six_j_category(F::Ring, names::Vector{String})
    C = SixJCategory()
    C.base_ring = F
    C.rank = length(names)
    C.simples_names = names
    set_pivotal!(C, [F(1) for _ ∈ names])
    
    return C
end
#-------------------------------------------------------------------------------
#   Setters/Getters
#-------------------------------------------------------------------------------

# Invalidate derived values, but preserve lazy F/R-symbol providers.
function _invalidate_sixj_structure!(C::SixJCategory;
                                     pivotal::Bool=true,
                                     fusion::Bool=true)
    if isdefined(C,:__attrs)
        keys = Symbol[:smatrix,:smatrix_by_objects,:twists_current,
                      :modular,:is_unitary]
        pivotal && append!(keys,[:pivotal,:spherical,:is_pivotal,:is_spherical])
        fusion && push!(keys,:multiplication_table)
        for key in keys
            delete!(getfield(C,:__attrs),key)
        end
    end
    nothing
end

@doc raw""" 

    set_tensor_product!(F::SixJCategory, mult::Array{Int,4})

Set the fusion rules of ``F``.
"""
function set_tensor_product!(F::SixJCategory, tensor)
    _invalidate_sixj_structure!(F)
    F.tensor_product = tensor
    n = size(tensor,1)

    ass = Array{MatElem,4}(undef,n,n,n,n)
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        ass[i,j,k,:] = matrices(id(F[i]⊗F[j]⊗F[k]))
    end
    F.ass = ass

end

function set_braiding!(F::SixJCategory, braiding)
    # Braiding changes twists, S, and modularity, but not the pivotal or
    # spherical structure of the underlying monoidal category.
    _invalidate_sixj_structure!(F;pivotal=false,fusion=false)
    F.braiding = braiding
end

function rank(C::SixJCategory)
    C.rank 
end

@doc raw""" 
    set_associator!(F::SixJCategory, ass::Array{MatElem,4})
    set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, ass::Vector{<:MatElem})
    set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int, ass::MatElem) 
    set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int, m::Int, n::Int, v::RingElem) 

Set the ``F``-symbols of ``F``. Pass `check=true` to validate normalized unit
associators; supplied data are trusted by default.
"""
# SixJCategory uses strictly normalised unit associators (see associator).
# Reject conflicting input instead of silently discarding a supplied F=2
# in the rank-one category. The pentagon itself would force F³=F².
function _check_normalised_unit_associator(C::SixJCategory,i,j,k,M)
    if isdefined(C,:one) && sum(C.one) == 1
        u = findfirst(!iszero,C.one)
        if u in (i,j,k)
            M == identity_matrix(base_ring(C),number_of_rows(M)) ||
                throw(ArgumentError("skeletal unit associators are normalised to identities"))
        end
    end
    nothing
end

function set_associator!(F::SixJCategory,ass; check::Bool=false)
    if check
        for i in 1:F.rank,j in 1:F.rank,k in 1:F.rank,l in 1:F.rank
            isassigned(ass,i,j,k,l) || continue
            _check_normalised_unit_associator(F,i,j,k,ass[i,j,k,l])
        end
    end
    _invalidate_sixj_structure!(F;fusion=false)
    F.ass = ass
end

function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int,
                         ass::Vector{<:MatElem}; check::Bool=false)
    _invalidate_sixj_structure!(F;fusion=false)
    check && foreach(M -> _check_normalised_unit_associator(F,i,j,k,M),ass)
    F.ass[i,j,k,:] = ass
end

function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int,
                         ass::MatElem; check::Bool=false)
    _invalidate_sixj_structure!(F;fusion=false)
    check && _check_normalised_unit_associator(F,i,j,k,ass)
    F.ass[i,j,k,l] = ass
end

function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int,
                         ass::Array{T,N}; check::Bool=false) where {T,N}
    set_associator!(F,i,j,k,l,
        matrix(base_ring(F),(N > 1 ? size(ass) : (1,1))...,ass);check)
end

function set_associator!(F::SixJCategory, i::Int, j::Int, k::Int, l::Int,
                         m::Int, n::Int, v::RingElem; check::Bool=false)
    if check
        M = deepcopy(F.ass[i,j,k,l])
        M[m,n] = v
        set_associator!(F,i,j,k,l,M;check=true)
    else
        _invalidate_sixj_structure!(F;fusion=false)
        F.ass[i,j,k,l][m,n] = v
    end
end

@doc raw""" 

    set_pivotal!(F::SixJCategory, p::Vector{<:RingElem})

Set the pivotal structure of ``F``. Warning: No checks are performed.
"""
function set_pivotal!(F::SixJCategory, sp)
    _invalidate_sixj_structure!(F;fusion=false)
    F.pivotal = sp
    # Structural data are trusted by default, just like supplied F-symbols.
    # `is_pivotal(F; check=true)` remains the explicit coherence certificate.
    set_attribute!(F,:pivotal,true)
end

function set_spherical!(C::SixJCategory,sp; check::Bool=false)
    if check
        base_ring(C) isa Union{ArbField,AcbField,ComplexField} &&
            throw(ArgumentError("a checked spherical setter requires exact coefficients"))
        is_spherical(C,sp) ||
            throw(ArgumentError("not a spherical pivotal structure"))
    end
    set_pivotal!(C,base_ring(C).(sp))
    set_attribute!(C,:spherical,true)
end

function is_spherical(C::SixJCategory,sp)
    length(sp) == rank(C) || return false
    sp = base_ring(C).(sp)
    any(iszero,sp) && return false
    old = C.pivotal
    try
        set_pivotal!(C,copy(sp))
        return is_spherical(C;check=true)
    finally
        set_pivotal!(C,old)
    end
end

function is_spherical(C::SixJCategory; check::Bool=false)
    check && return _sixj_is_spherical(C)
    get_attribute!(() -> _sixj_is_spherical(C),C,:spherical)
end

function _sixj_is_spherical(C::SixJCategory)
    isdefined(C,:pivotal) && length(C.pivotal) == rank(C) || return false
    any(iszero,C.pivotal) && return false
    # Ising [1,1,2] has equal dimensions on duals but is NOT pivotal:
    # the σ⊗σ summands force the σ component to square to one.
    # EGNO, Definitions 4.7.7 and 4.7.14.
    is_pivotal(C;check=true) || return false
    if base_ring(C) isa Union{ArbField,AcbField,ComplexField}
        return all(overlaps(dim(X),dim(dual(X))) for X in simples(C))
    end
    all(dim(X) == dim(dual(X)) for X in simples(C))
end

"""
    set_canonical_spherical!(C::SixJCategory; embedding=nothing)

Set the canonical spherical structure of a pseudounitary category whose
dimensions are the FP dimensions. Pass `check=true` to verify the resulting
pivotal and spherical identities before storing it.
This requires exact characteristic-zero coefficients. For a number field,
specify the complex embedding that determines positivity. Failure preserves
the previous pivotal data. Existence is not automatic; see EGNO Proposition
9.5.1.
"""
function set_canonical_spherical!(C::SixJCategory; embedding=nothing,
                                  check::Bool=false)
    K = base_ring(C)
    K isa Union{ArbField,AcbField,ComplexField} &&
        throw(ArgumentError("canonical spherical certification requires exact coefficients"))
    characteristic(K) == 0 && is_fusion(C) ||
        throw(ArgumentError("canonical spherical normalization requires a characteristic-zero fusion category"))
    if K isa NumField
        embedding === nothing && isdefined(C,:embedding) &&
            (embedding = C.embedding)
        embedding isa AbsSimpleNumFieldEmbedding &&
            number_field(embedding) === K ||
            throw(ArgumentError("specify a complex embedding of the coefficient number field"))
        embedding = _qqbar_embedding(embedding)
    end
    target = [_fpdim_in_base_field(K,fpdim(s),embedding) for s in simples(C)]
    old = C.pivotal
    candidate = try
        set_pivotal!(C,fill(K(1),rank(C)))
        target ./ dim.(simples(C))
    finally
        set_pivotal!(C,old)
    end
    set_spherical!(C,candidate;check)
    copy(C.pivotal)
end

function _fpdim_in_base_field(K,d,embedding)
    K isa NumField || return K(d)
    candidates = roots(change_base_ring(K,minpoly(d)))
    isempty(candidates) &&
        throw(ArgumentError("the FP dimension is not in the coefficient field"))
    matches = filter(c -> embedding(c) == d,candidates)
    length(matches) == 1 ||
        throw(ArgumentError("the FP dimension is not in the specified embedded field"))
    only(matches)
end


@doc raw""" 

    set_one!(F::SixJCategory, v::Vector{Int})
    set_one!(F::SixJCategory, i::Int)   

Set the unit of ``F``. Pass `check=true` to scan stored associators for the
normalized unit constraint.
"""
function set_one!(F::SixJCategory,v::Vector; check::Bool=false)
    length(v) == rank(F) && all(c -> c isa Integer && c >= 0,v) ||
        throw(ArgumentError("invalid unit multiplicities"))
    # Unit-normalised associators must be valid even when the unit is set
    # AFTER an entire array of F-symbols (the usual constructor order).
    if check && sum(v) == 1 && isdefined(F,:ass)
        u = findfirst(!iszero,v)
        for i in 1:rank(F),j in 1:rank(F),k in 1:rank(F),l in 1:rank(F)
            if u in (i,j,k)
                isassigned(F.ass,i,j,k,l) || continue
                M = F.ass[i,j,k,l]
                M == identity_matrix(base_ring(F),number_of_rows(M)) ||
                    throw(ArgumentError("skeletal unit associators are normalised to identities"))
            end
        end
    end
    _invalidate_sixj_structure!(F;fusion=false)
    F.one = copy(v)
end

function set_one!(F::SixJCategory,i::Int; check::Bool=false)
    1 <= i <= rank(F) || throw(ArgumentError("unit index out of range"))
    set_one!(F,[Int(k == i) for k in 1:rank(F)];check)
end

function set_ribbon!(F::SixJCategory, r)
    F.ribbon = r
end

function set_twist!(F::SixJCategory, t)
    _invalidate_sixj_structure!(F;pivotal=false,fusion=false)
    F.twist = t
    set_attribute!(F,:twists_current,true)
end

@doc raw""" 

    set_name!(F::SixJCategory, name::String)    

Set the display name of ``F``.
"""
function set_name!(F::SixJCategory, name)
    F.name = name
end

function set_simples_names!(F::SixJCategory, names::Vector{String})
    F.simples_names = names
end

simples_names(C::SixJCategory) = C.simples_names
indecomposables_names(C::SixJCategory) = C.simples_names
multiplication_table(C::SixJCategory) = C.tensor_product
#(::Type{Int})(x::QQFieldElem) = Int(numerator(x))

function braiding(X::SixJObject, Y::SixJObject) 
    C = parent(X)
    if is_simple(X) && is_simple(Y)
        i = findfirst(e -> e != 0, X.components)
        j = findfirst(e -> e != 0, Y.components)

        if ! all(isassigned(C.braiding, i,j,k) for k ∈ 1:C.rank)
            r_symbol = get_attribute(C, :r_symbol)
            C.braiding[i,j,:] = [r_symbol(i,j,k) for k ∈ 1:C.rank]
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

    if one(C) ∈ [X,Y,Z]
        return id(X ⊗ Y ⊗ Z)
    end

    F = base_ring(C)
    n = C.rank
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
    #@assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"
    C = parent(X)

    if zero(C) == X ⊗ Y ⊗ Z
        return zero_morphism(zero(C),zero(C))
    end

    if one(C) ∈ [X,Y,Z]
        return id(X ⊗ Y ⊗ Z)
    end

    F = base_ring(C)
    n = C.rank
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

    _,ix,px = direct_sum(X_summands)
    _,iy,py = direct_sum(Y_summands)
    _,iz,pz = direct_sum(Z_summands)

    distr_before= horizontal_direct_sum([(f⊗g)⊗h for f ∈ ix, g ∈ iy, h ∈ iz][:])
    distr_after = vertical_direct_sum([f⊗(g⊗h) for f ∈ px, g ∈ py, h ∈ pz][:])

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

# function vector_permutation(A::Vector,B::Vector)
#     perm = Int[]
#     for a ∈ A
#         i = findall(e -> e == a, B)
#         j = filter(e -> !(e ∈ perm), i)[1]
#         perm = [perm; j]
#     end
#     return perm
# end

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
        C.braiding[i,j,k] = r_symbol(i,j,k)
    end
    return C.braiding[i,j,k]
end

# Capture the old array and provider rather than the category being mutated.
# Assigned blocks are reused, while missing blocks remain deferred.
function _sixj_symbol_reader(data, provider)
    function read(indices...)
        if !isassigned(data,indices...)
            data[indices...] = provider(indices...)
        end
        data[indices...]
    end
end

# Some consumers, such as coefficient reduction through an unspecified
# embedding, genuinely require every structural block at once.
function _materialize_sixj_symbols!(C::SixJCategory)
    n = rank(C)
    for i in 1:n,j in 1:n,k in 1:n,l in 1:n
        six_j_symbol(C,i,j,k,l)
    end
    if is_braided(C)
        for i in 1:n,j in 1:n,k in 1:n
            r_symbol(C,i,j,k)
        end
    end
    C
end

function (K::AcbField)(f::SixJMorphism)
    domain(f) == codomain(f) || throw(ArgumentError("a scalar morphism must be an endomorphism"))
    M = matrix(f)
    n = number_of_rows(M)
    n > 0 || throw(ArgumentError("the scalar of an endomorphism of zero is not unique"))
    x = sum(M[i,i] for i in 1:n) / K(n)
    all(_acb_approximately_equal(M[i,j],i == j ? x : zero(K))
        for i in 1:n,j in 1:n) ||
        throw(ArgumentError("not a numerically represented scalar matrix"))
    K(x)
end

function _acb_approximately_equal(x::AcbFieldElem,y::AcbFieldElem)
    overlaps(x,y) && return true
    p = precision(parent(x))
    A = parent(abs(x))
    guard = min(32,div(p,2))
    abs(x-y) < (one(A) + abs(x) + abs(y)) * A(2)^(-p+guard)
end

#-------------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------------
is_semisimple(::SixJCategory) = true
# Objects have finite multiplicity vectors over a fixed finite simple list.
is_finite(::SixJCategory) = true
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
function ==(f::SixJMorphism,g::SixJMorphism)
    domain(f) == domain(g) && codomain(f) == codomain(g) || return false
    if base_ring(f) isa Union{ArbField,AcbField,ComplexField}
        # Equality of represented enclosures is transitive; overlap is not.
        return Base.isequal(f.m,g.m)
    end
    f.m == g.m
end
function overlaps(f::SixJMorphism,g::SixJMorphism)
    domain(f) == domain(g) && codomain(f) == codomain(g) && overlaps(matrix(f),matrix(g))
end


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

    mats = [hcat([g.m[i] for g ∈ f]...) for i ∈ 1:C.rank]

    morphism(domain(f[1]), cod, mats)
end

function horizontal_direct_sum(f::Vector{SixJMorphism})
    @assert all(codomain(g) == codomain(f[1]) for g ∈ f[2:end])

    C = parent(f[1])
    dom = SixJObject(C, reduce(.+, [X.components for X ∈ domain.(f)]))

    mats = [vcat([g.m[i] for g ∈ f]...) for i ∈ 1:C.rank]

    morphism(dom, codomain(f[1]), mats)
end

function +(f::SixJMorphism, g::SixJMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    SixJMorphism(domain(f), codomain(f), [m + n for (m,n) ∈ zip(f.m,g.m)])
end

function tr(f::SixJMorphism)
    domain(f) == codomain(f) ||
        throw(ArgumentError("trace requires an endomorphism"))
    C = parent(f)
    result = zero_morphism(one(C),one(C))
    # On the isotypic block S_i^n, the pivotal trace is the ordinary matrix
    # trace times dim(S_i); see EGNO, Definition 4.7.1 and Proposition 4.7.3.
    for (S,M) in zip(simples(C),matrices(f))
        t = tr(M)
        iszero(t) || (result += t*left_trace(id(S)))
    end
    result
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
        for k ∈ 1:C.rank 
            if C.one[k] == 1
                j = [j; findall(e -> C.tensor_product[i,e,k] >= 1, 1:C.rank)]
            end
        end
        if length(j) != 1
            throw(ErrorException("Object not rigid."))
        end
        return SixJObject(C,[i == j[1] ? 1 : 0 for i ∈ 1:C.rank])
    end

    # Build dual from simple objects
    return direct_sum([dual(Y)^(X.components[i]) for (Y,i) ∈ zip(simples(C), 1:C.rank)])[1]
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

    unscaled_ev = basis(Hom(DX ⊗ X, one(C)))[1]
    unscaled_coev = basis(Hom(one(C), cod))[1]

    factor = F((id(X)⊗unscaled_ev)∘associator(X,DX,X)∘(unscaled_coev⊗id(X)))

    # Delete close to zero part if inexact
    if base_ring(C) isa Union{AcbField,ComplexField,ArbField}
        if overlaps(F(real(factor)), zero(base_ring(C)))
            factor = F(imag(factor))
        elseif overlaps(F(imag(factor)), zero(base_ring(C)))
            factor = F(real(factor))
        end
    end

    if base_ring(C) isa Union{QQBarField, ComplexField, AcbField, ArbField}
        return inv(sqrt(factor)) * unscaled_coev
    end

    return unscaled_coev
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

function spherical(X::SixJObject; check::Bool=false)
    check && !is_spherical(parent(X);check=true) && throw(ArgumentError("Not spherical"))
    pivotal(X)
end


function pivotal(X::SixJObject)
    C = parent(X)
    F = base_ring(C)
    sp = C.pivotal
    mats = [diagonal_matrix(θ, k) for (θ,k) ∈ zip(sp, X.components)]
    return morphism(X,X,mats)
end

"""
    twists(C::SixJCategory; check=false)

Return the twist scalars of the supplied braided pivotal structure. The
structure is assumed valid; `check=true` verifies pivotal coherence. Values
are cached until structural data are changed through a setter.
"""
function twists(C::SixJCategory; check::Bool=false)
    check && !is_pivotal(C;check=true) &&
        throw(ArgumentError("a pivotal structure is required"))
    get_attribute(C,:twists_current,false) && return C.twist
    K = base_ring(C)
    # theta_X = u_X^{-1}j_X (EGNO, Section 8.10). On a split simple both
    # maps are scalars, so divide their coefficients directly.
    C.twist = [C.pivotal[i] / K(drinfeld_morphism(X))
               for (i,X) in enumerate(simples(C))]
    set_attribute!(C,:twists_current,true)
    C.twist
end

# The cache key records the ordered object coordinates. Setters invalidate the
# cache, and a copy prevents callers from corrupting the stored matrix.
function smatrix(C::SixJCategory, objects=simples(C))
    all(X -> parent(X) === C,objects) ||
        throw(ArgumentError("objects outside the category"))
    key = Tuple(Tuple(X.components) for X in objects)
    cache = get_attribute!(() -> Dict{Tuple,MatElem}(),C,:smatrix_by_objects)
    S = get!(cache,key) do
        invoke(smatrix,Tuple{Category,typeof(objects)},C,objects)
    end
    deepcopy(S)
end

function twist(X::SixJObject)
    C = parent(X)
    F = base_ring(C)
    sp = twists(C)
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
    _check_sixj_parents(X, Y)
    C = parent(X)
    n = C.rank
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

    for i ∈ 1:C.rank, j ∈ 1:C.rank
       # @show parent(f.m[i]), parent(g.m[j])
        A = kronecker_product(f.m[i],g.m[j])
        d1,d2 = size(A)
        #if d1*d2 == 0 continue end
        for k ∈ 1:C.rank
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

#     for i ∈ 1:parent(X).rank
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
    _check_sixj_parents(X...)
    if length(X) == 1
        return X...,[id(X...)], [id(X...)]
    end

    S = SixJObject(parent(X[1]), vec(sum(hcat([x.components for x in X]...), dims = 2)))

    inc = [matrices(zero_morphism(x,S)) for x ∈ X]
    proj = [matrices(zero_morphism(S,x)) for x ∈ X]

    for i ∈ 1:parent(X[1]).rank
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
    _check_sixj_parents(X...)
    SixJObject(parent(X[1]), vec(sum(hcat([x.components for x in X]...), dims = 2)))
end

function ^(X::SixJObject, k::Int)
    SixJObject(parent(X), k.*(X.components))
end

function direct_sum(f::SixJMorphism...)
    dom = ⊕(domain.(f)...)
    cod = ⊕(codomain.(f)...)
    F = base_ring(dom)

    mats = [diagonal_matrix([g.m[i] for g ∈ f]) for i ∈ 1:parent(dom).rank]

    return morphism(dom,cod, mats)
end

# function vertical_direct_sum(f::Vector{SixJMorphism})
#     if length(f) == 1
#         return f[1]
#     end
    
#     #@assert length(unique!([domain.(f)...])) == 1 "Not compatible"

#     ms = matrices.(f)
#     m = [hcat([n[i] for n ∈ ms]...) for i ∈ 1:parent(f[1]).rank]
#     return morphism(domain(f[1]), ⊕(codomain.(f)...), m)
# end

# function horizontal_direct_sum(f::Vector{SixJMorphism})
#     if length(f) == 1
#         return f[1]
#     end
#     # @assert length(unique!([codomain.(f)...])) == 1 "Not compatible"

#     ms = matrices.(f)
#     m = [vcat([n[i] for n ∈ ms]...) for i ∈ 1:parent(f[1]).rank]
#     return morphism(⊕(domain.(f)...), codomain(f[1]), m)
# end


zero(C::SixJCategory) = SixJObject(C,[0 for i ∈ 1:C.rank])

function zero_morphism(X::SixJObject, Y::SixJObject)
    _check_sixj_parents(X, Y)
    return SixJMorphism(X,Y,[zero(matrix_space(base_ring(X), cX, cY)) for (cX,cY) ∈ zip(X.components, Y.components)])
end

function is_isomorphic(X::SixJObject, Y::SixJObject)
    if X != Y
        return false, nothing
    else
        return true, morphism(X, Y, matrices(id(X)))
    end
end
#-------------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------------

function simples(C::SixJCategory)
    n = C.rank
    [SixJObject(C, [i == j ? 1 : 0 for j ∈ 1:n]) for i ∈ 1:n]
end


"""
    sort_simples!(C::SixJCategory, order::Vector{Int})

Relabel simple `order[i]` as `i`, transporting the intermediate fusion-channel
bases as well as the outer F- and R-symbol indices. Deferred symbols remain
deferred. Existing objects and functors are not transported and should not be
reused after this mutation.
"""
function sort_simples!(C::SixJCategory, order::Vector{Int})
    n = C.rank
    structure_flags = Dict(key => get_attribute(C,key)
        for key in (:pivotal,:spherical) if has_attribute(C,key))
    sort(order) == collect(1:n) ||
        throw(ArgumentError("order must be a permutation of the simple labels"))
    order = copy(order)
    N = C.tensor_product
    old_ass = C.ass
    readF = _sixj_symbol_reader(old_ass,
                                get_attribute(C,:six_j_symbol,nothing))
    transportedF = function(i,j,k,l)
        a,b,c,d = order[i],order[j],order[k],order[l]
        rows = _fusion_channel_permutation(
            [N[a,b,e]*N[e,c,d] for e in 1:n],order)
        cols = _fusion_channel_permutation(
            [N[b,c,f]*N[a,f,d] for f in 1:n],order)
        readF(a,b,c,d)[rows,cols]
    end
    ass = Array{MatElem,4}(undef,n,n,n,n)
    for I in CartesianIndices(ass)
        i,j,k,l = Tuple(I)
        isassigned(old_ass,order[i],order[j],order[k],order[l]) || continue
        ass[I] = transportedF(i,j,k,l)
    end
    if isdefined(C,:braiding)
        old_R = C.braiding
        readR = _sixj_symbol_reader(old_R,get_attribute(C,:r_symbol,nothing))
        transportedR = (i,j,k) -> readR(order[i],order[j],order[k])
        R = Array{MatElem,3}(undef,n,n,n)
        for I in CartesianIndices(R)
            i,j,k = Tuple(I)
            isassigned(old_R,order[i],order[j],order[k]) || continue
            R[I] = transportedR(i,j,k)
        end
        C.braiding = R
        set_attribute!(C,:r_symbol,transportedR)
    end
    C.tensor_product = N[order,order,order]
    C.ass = ass
    set_attribute!(C,:six_j_symbol,transportedF)
    C.simples_names = C.simples_names[order]

    isdefined(C, :one) && (C.one = C.one[order])
    isdefined(C, :pivotal) && (C.pivotal = C.pivotal[order])
    isdefined(C, :twist) && (C.twist = C.twist[order])
    _invalidate_sixj_structure!(C)
    for (key,value) in structure_flags
        set_attribute!(C,key,value)
    end
    return C
end

function _fusion_channel_permutation(sizes,order)
    offsets = cumsum([0; sizes])
    [r for e in order for r in offsets[e]+1:offsets[e+1]]
end

function sort_simples_by_dimension!(C::SixJCategory)
    sort_simples!(C, sortperm(fpdim.(simples(C))))
end

#-------------------------------------------------------------------------------
#   Kernel and Cokernel
#-------------------------------------------------------------------------------

function kernel(f::SixJMorphism)
    C = parent(domain(f))
    kernels = [kernel(m, side = :left) for m ∈ f.m]
    
    ker = SixJObject(C,[number_of_rows(k) for k ∈ kernels])
    k = morphism(ker, domain(f), [m for m ∈ kernels])
    if is_unitary(C)
        k = orthonormalization(k)
    end
    return ker, k
end

function cokernel(f::SixJMorphism)
    C = parent(domain(f))

    cokernels = [kernel(m, side = :right) for m ∈ f.m]
    
    coker = SixJObject(C,[number_of_columns(k) for k ∈ cokernels])
    c = morphism(codomain(f),coker, [m for m ∈ cokernels])

    if is_unitary(C)
        c = dagger(orthonormalization(dagger(c)))
    end
    return coker, c
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
    _check_sixj_parents(X, Y)
    Xi, Yi = X.components, Y.components
    F = base_ring(X)

    d = sum([x*y for (x,y) ∈ zip(Xi,Yi)])

    if d == 0 return SixJHomSpace(X,Y,SixJMorphism[]) end

    basis = [zero_morphism(X,Y).m for i ∈ 1:d]
    next = 1
    for k ∈ 1:parent(X).rank

        for i ∈ 1:Xi[k], j ∈ 1:Yi[k]
            basis[next][k][i,j] = 1
            next = next + 1
        end
    end
    basis_mors = [SixJMorphism(X,Y,m) for m ∈ basis]
    return SixJHomSpace(X,Y,basis_mors)
end

function express_in_basis(f::SixJMorphism, H::SixJHomSpace)
    domain(f) == domain(H) && codomain(f) == codomain(H) ||
        throw(ArgumentError("morphism and Hom basis must have the same endpoints"))
    # Hom constructs matrix units in block, row, column order. Julia's
    # column-major vec would permute the coefficients within rectangular blocks.
    [M[i,j] for M in matrices(f)
            for i in 1:number_of_rows(M) for j in 1:number_of_columns(M)]
end


#-------------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------------

function show(io::IO, C::SixJCategory)
    if isdefined(C,:name)
        print(io, "$(C.name)")
    else
      print(io, "Fusion Category with $(C.rank) simple objects")
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
function extension_of_scalars(C::SixJCategory, L::Ring;
                              embedding = _scalar_extension_embedding(base_ring(C),L))
    # Avoid placeholder associators and do not force deferred structural data.
    D = SixJCategory()
    D.base_ring = L
    D.rank = rank(C)
    D.tensor_product = copy(C.tensor_product)
    D.simples_names = copy(simples_names(C))
    isdefined(C,:name) && set_name!(D,C.name)
    if isdefined(C,:ass)
        D.ass = _transport_sixj_array!(D,C,:ass,:six_j_symbol,L,embedding)
    end
    if isdefined(C,:one)
        D.one = copy(C.one)
    end
    if isdefined(C,:pivotal)
        D.pivotal = embedding.(C.pivotal)
        L isa Union{ArbField,AcbField} && (D.pivotal = L.(D.pivotal))
        for key in (:pivotal,:spherical)
            has_attribute(C,key) && set_attribute!(D,key,get_attribute(C,key))
        end
    end
    if isdefined(C,:braiding)
        D.braiding = _transport_sixj_array!(D,C,:braiding,:r_symbol,L,embedding)
    end
    if get_attribute(C,:twists_current,false)
        set_twist!(D,embedding.(C.twist))
    end
    D
end

function _transport_sixj_array!(D,C,field,key,L,embedding)
    data = getfield(C,field)
    read = _sixj_symbol_reader(data,get_attribute(C,key,nothing))
    convert_matrix = a -> matrix(L,size(a)...,embedding.(collect(a)))
    out = similar(data)
    for I in CartesianIndices(data)
        isassigned(data,Tuple(I)...) || continue
        out[I] = convert_matrix(data[I])
    end
    set_attribute!(D,key,(indices...) -> convert_matrix(read(indices...)))
    out
end

complex_embedding_of_base_ring(C::SixJCategory) = C.embedding

function extension_of_scalars(C::SixJCategory, K::FqField; embedding=nothing)
    if embedding !== nothing || is_finite(base_ring(C))
        e = embedding === nothing ?
            _scalar_extension_embedding(base_ring(C),K) : embedding
        return invoke(extension_of_scalars,Tuple{SixJCategory,Ring},C,K;
                      embedding=e)
    end
    _materialize_sixj_symbols!(C)
    denom = if base_ring(C) == QQ 
        lcm([isempty(m) ? ZZ(1) : lcm(denominator.(collect(m))[:]) for m ∈ C.ass][:])
    else 
        lcm([isempty(m) ? ZZ(1) : lcm(denominator.(hcat(coefficients.(collect(m))[:]...))) for m ∈ C.ass][:])
    end

    gcd(denom, characteristic(K)) > 1 && error("Not able to define over $K. $denom is not coprime to $(characteristic(K))")

    rs = base_ring(C) == QQ ? [nothing] : roots(change_base_ring(K, minpoly(gen(base_ring(C)))))

    length(rs) == 0 && error("Not able to define over $K. \n Minpoly has no roots: $(minpoly(gen(base_ring(C)))) \n Denomitator is $(denom)")
    
    conv = if base_ring(C) == QQ 
        [x -> K(numerator(x))//K(denominator(x))]
    else 
        [x -> sum([K(numerator(c))//K(denominator(c)) * r^(i-1) for (i,c) in enumerate(coefficients(x)) ]) for r ∈ rs]
    end

    for c ∈ conv 
        F = extension_of_scalars(C, K, embedding = c)
        pentagon_axiom(F) && return F
    end

    error("Not able to define over $K.")
end

@doc raw""" 

    complex_embedding(C::SixJCategory)
    complex_embedding(C::SixJCategory, e::AbsSimpleNumFieldEmbedding)

Return the complex embedding of C if an embedding of the ground field is specified or given.
"""
function complex_embedding(C::SixJCategory)
    !isdefined(C, :embedding) && error("No embedding has been specified")
    complex_embedding(C, getfield(C, :embedding))
end

function complex_embedding(C::SixJCategory, e::AbsSimpleNumFieldEmbedding)
    extension_of_scalars(C, QQBarField(), e)
end

@doc raw""" 

    complex_embeddings(C::SixJCategory)

Return all complex_embeddings of C.
"""
function complex_embeddings(C::SixJCategory)
    [extension_of_scalars(C, QQBarField(), e) for e ∈ complex_embeddings(base_ring(C))]
end


function extension_of_scalars(C::SixJCategory, K::QQBarField; embedding=nothing)
    if embedding isa AbsSimpleNumFieldEmbedding
        return extension_of_scalars(C,K,embedding)
    end
    if embedding !== nothing || base_ring(C) isa Union{QQField,QQBarField}
        e = embedding === nothing ? K : embedding
        return invoke(extension_of_scalars,Tuple{SixJCategory,Ring},C,K;
                      embedding=e)
    end
    e = isdefined(C,:embedding) ? C.embedding :
        first(complex_embeddings(base_ring(C)))
    extension_of_scalars(C,K,e)
end

function extension_of_scalars(C::SixJCategory, K::QQBarField,
                              e::AbsSimpleNumFieldEmbedding)
    number_field(e) === base_ring(C) ||
        throw(ArgumentError("embedding has the wrong source field"))
    invoke(extension_of_scalars,Tuple{SixJCategory,Ring},C,K;
           embedding=_qqbar_embedding(e))
end

# Choose one algebraic image of the primitive element using the certified
# complex embedding, then evaluate every coefficient through that same root.
# This gives a field homomorphism rather than independent algebraic guesses.
function _qqbar_embedding(e::AbsSimpleNumFieldEmbedding)
    L = number_field(e)
    K = QQBarField()
    a = gen(L)
    candidates = roots(K,minpoly(a))
    for prec in (64,128,256,512,1024,2048,4096)
        image = e(a,prec)
        A = AcbField(prec)
        matches = filter(r -> overlaps(image,A(r)),candidates)
        if length(matches) == 1
            r = only(matches)
            return x -> foldr((c,y) -> K(c)+r*y,coefficients(L(x));
                              init=K(0))
        end
    end
    throw(ArgumentError("could not isolate the specified number-field embedding"))
end

""" 

    extension_of_scalars(X::SixJObject, K::Field)

Return the object ``X`` as an object of the category ``C⊗K``.
"""
function extension_of_scalars(X::SixJObject, L::Ring, CL::SixJCategory; embedding = nothing)
    if CL === nothing
            CL = extension_of_scalars(parent(m), L, embedding = embedding)
    end
    SixJObject(CL, X.components)
end

""" 

    extension_of_scalars(f::SixJMorphism, K::Field)

Return the category ``C⊗K``.
"""
function extension_of_scalars(m::SixJMorphism, L::Ring, CL::SixJCategory;
                              embedding = _scalar_extension_embedding(base_ring(m),L))
    try 
        if CL === nothing
            CL = extension_of_scalars(parent(m), L, embedding = embedding)
        end

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
    if is_tambara_yamagami(C) && is_modular(C)
        return tambara_yamagami_tensor_autoequivalences(C)
    end

    fusion_ring_autos = automorphisms(split_grothendieck_ring(C))

    equivs = MonoidalFunctor[]

    for f ∈ fusion_ring_autos 
        images = [findfirst(==(a), f.images) for a in basis(domain(f))]
        F = functor(C,C, simples(C)[images])
        append!(equivs, monoidal_structures(F))
    end
    equivs
end

"""
    autoequivalence_candidates(C; check=false)

Return tensor autoequivalence candidates; completeness is not claimed.  The
solved equations impose coherence, and `check=true` additionally re-evaluates
the monoidal functor axioms on every returned candidate.
"""
function autoequivalence_candidates(C::SixJCategory; check::Bool=false)
    fusion_ring_autos = automorphisms(split_grothendieck_ring(C))
    result = MonoidalFunctor[]
    for f in fusion_ring_autos
        images = [findfirst(==(a),f.images) for a in basis(domain(f))]
        F = functor(C,C,simples(C)[images])
        append!(result,monoidal_structure_candidates(F;check))
    end
    result
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
    Unitary
----------------------------------------------------------=#    

"""
    is_unitary(C::SixJCategory)

Test whether the supplied structure is unitary in its stored bases. For ball
coefficients, the equations are tested at the precision of the base field.
"""
function is_unitary(C::SixJCategory)
    get_attribute!(C, :is_unitary) do
        K = base_ring(C)
        K isa Union{ArbField,AcbField,ComplexField} &&
            return is_unitary_numeric(C)
        K isa Union{QQField,NumField,FqField} && return false

        !is_spherical(C;check=true) && return false
        !all(fpdim(s) == dim(s) for s in simples(C)) && return false

        for x ∈ simples(C), y ∈ simples(C), z ∈ simples(C) 
            !is_unitary(associator(x,y,z)) && return false 
        end 
        true
    end
end

function is_unitary(f::SixJMorphism)
    base_ring(f) isa Union{ArbField,AcbField,ComplexField} &&
        return is_unitary_numeric(f)
    !is_invertible(f) && return false
    f ∘ dagger(f) == id(codomain(f))    
end

"""
    is_unitary_numeric(f::SixJMorphism)
    is_unitary_numeric(C::SixJCategory)

Test compatibility of ball enclosures with the unitarity equations in the
stored bases. This is numerical evidence, not an exact certificate or a test
for the existence of another unitary gauge.
"""
function is_unitary_numeric(f::SixJMorphism)
    base_ring(f) isa Union{ArbField,AcbField,ComplexField} ||
        throw(ArgumentError("ball coefficients required"))
    domain(f).components == codomain(f).components || return false
    overlaps(f ∘ dagger(f),id(codomain(f))) &&
        overlaps(dagger(f) ∘ f,id(domain(f)))
end

function is_unitary_numeric(C::SixJCategory)
    base_ring(C) isa Union{ArbField,AcbField,ComplexField} ||
        throw(ArgumentError("ball coefficients required"))
    is_spherical(C;check=true) || return false
    S = simples(C)
    all(s -> overlaps(fpdim(s),dim(s)),S) || return false
    all(is_unitary_numeric(associator(x,y,z)) for x in S,y in S,z in S)
end

function dagger(f::SixJMorphism)
    mats = [transpose(conj.(m)) for m ∈ matrices(f)]
    morphism(codomain(f), domain(f), mats)
end

function inner_product(f::SixJMorphism, g::SixJMorphism)
   base_ring(f)(tr(dagger(f) ∘ g))//dim(domain(g))
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

#=----------------------------------------------------------
    subcategories
----------------------------------------------------------=#    

function fusion_subcategory(X::SixJObject)
    C = parent(X)
    R = split_grothendieck_ring(C)
    basis = _topologize(R(ZZ.(X.components)))
    S = six_j_category(base_ring(X),
        multiplication_table(C)[basis, basis, basis],
        simples_names(C)[basis])

    set_associator!(S, associator(C)[basis, basis, basis, basis])
    set_one!(S, C.one[basis])
    isdefined(C, :pivotal) && set_pivotal!(S, C.pivotal[basis])
    isdefined(C, :one) && set_one!(S, C.one[basis])
    isdefined(C, :braiding) && set_braiding!(S, C.braiding[basis, basis, basis])
    isdefined(C, :embedding) && (S.embedding = C.embedding)
    set_name!(S, "Fusion subcategory of $(C.name) generated by $X")
    return S
end

#=----------------------------------------------------------
    Export F-symbols as Dict 
----------------------------------------------------------=#

# Dictionary conventions change coordinates at the API boundary, never the
# associator/braiding matrices, fusion bases, or the category itself.
function _check_symbol_convention(convention::Symbol)
    convention in (:column_major_packing, :bonderson) ||
        throw(ArgumentError("unknown symbol convention $convention; use :column_major_packing or :bonderson"))
    convention
end

@doc raw"""
    F_symbols(C::SixJCategory; convention=:column_major_packing)

Return a dictionary containing all F-symbol entries, including zeros.

* `:bonderson`: keys are `[a,b,c,d,e,f]` in a multiplicity-free category,
  otherwise `[a,b,c,d,e,μ,ν,f,ρ,σ]`. The left fusion path is
  ``a b \xrightarrow{\mu} e,\ e c \xrightarrow{\nu} d``;
  the right path is ``b c \xrightarrow{\rho} f,\ a f \xrightarrow{\sigma} d``.
  Writing ``L_u`` and ``R_v`` for these projection trees and
  ``\alpha:(a\otimes b)\otimes c\to a\otimes(b\otimes c)``, the entry is
  ``A_{u v}`` in ``R_v\circ\alpha=\sum_u A_{u v}L_u``.
  Equivalently, in composition-dual splitting bases,
  ``\alpha\circ L^u=\sum_v A_{u v}R^v``. These are the coefficients of
  Bonderson, *Non-Abelian Anyons and Interferometry* (2007), Eq. (2.6), p. 14,
  https://thesis.caltech.edu/2447/02/thesis.pdf.
  Rows of the stored matrix are ordered `(e,μ,ν)`, columns `(f,ρ,σ)`,
  with the last index varying fastest.
* `:column_major_packing` (default): preserves the historical dictionary
  layout. In each `(a,b,c,d)` block, consume the stored matrix down columns
  while traversing keys `[a,b,c,d,f,e]` in `(e,f)` order, or keys
  `[a,b,c,d,f,σ,ρ,e,μ,ν]` in `(e,f,ν,μ,ρ,σ)` order. These key positions
  do not, in general, label the matrix entry's fusion paths. Conversion
  requires unpacking the block; a fixed permutation of key positions is
  insufficient.

The result is a plain `Dict`: it does not record its convention. Keep track of
it when passing dictionaries to other functions. Neither choice changes gauge,
inverts the associator, nor changes the internal matrices. These coordinates
use the split skeletal bases of `SixJCategory`.
"""
function F_symbols(C::SixJCategory; convention::Symbol=:column_major_packing)
    _check_symbol_convention(convention)
    convention == :bonderson && return _bonderson_F_symbols(C)

    S = simples(C)

    N = length(S)
    C_morphism_type = morphism_type(C)
    K = base_ring(C) 

    F_dict = Dict{Vector{Int}, elem_type(K)}()

    one_indices = findall(s -> int_dim(Hom(s,one(C))) > 0 , S)
    one_components = simple_subobjects(one(C))

    # Set unitors to identity
    S[one_indices] = one_components

    prods = [X ⊗ Y for X ∈ S, Y ∈ S]

    homs = [basis(Hom(prods[i,j],S[k])) for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N]

    m = multiplicity(C) 
    mult = multiplication_table(C)

    for (a,b,c,d) ∈ collect(Base.product(1:N, 1:N, 1:N, 1:N))

        sym = collect(six_j_symbol(C,a,b,c,d))[:]
  
        for e in 1:N, f in 1:N 
            if mult[a,b,e] * mult[e,c,d] * mult[b,c,f] * mult[a,f,d] == 0 
                continue 
            end
            if m == 1
                # F = K(compose(
                #     Hom(C[d], C[e] ⊗ C[c])[1],
                #     Hom(C[e], C[a] ⊗ C[b])[1] ⊗ id(C[c]),
                #     associator(C[[a,b,c]]...),
                #     id(C[a]) ⊗ Hom(C[b] ⊗ C[c], C[f])[1],
                #     Hom(C[a] ⊗ C[f], C[d])[1]
                # ))
                F_dict[[a,b,c,d,f,e]] = popfirst!(sym)
            else
                for (i,g) in pairs(basis(Hom(C[d], C[e] ⊗ C[c])))
                    for (j,g2) in pairs(basis(Hom(C[e], C[a] ⊗ C[b])))
                        for (i2,h) in pairs(basis(Hom(C[b] ⊗ C[c], C[f])))
                            for (j2,h2) in pairs(basis(Hom(C[a] ⊗ C[f], C[d])))
                                # F = K(compose(
                                #     g,
                                #     g2 ⊗ id(C[c]),
                                #     associator(C[[a,b,c]]...),
                                #     id(C[a]) ⊗ h,
                                #     h2
                                # ))
                                #a == b == c == d == 2 && println((f,j2,i2,e,j,i))
                                F_dict[[a,b,c,d,f,j2,i2,e,j,i]] = popfirst!(sym)
                            end
                        end
                    end
                end
            end
        end
    end

    return F_dict           
end 

function _bonderson_F_symbols(C::SixJCategory)
    N = multiplication_table(C)
    n = rank(C)
    mf = multiplicity(C) == 1
    F = Dict{Vector{Int},elem_type(base_ring(C))}()
    for a in 1:n, b in 1:n, c in 1:n, d in 1:n
        A = six_j_symbol(C,a,b,c,d)
        row = 0
        for e in 1:n, μ in 1:N[a,b,e], ν in 1:N[e,c,d]
            row += 1
            col = 0
            for f in 1:n, ρ in 1:N[b,c,f], σ in 1:N[a,f,d]
                col += 1
                key = mf ? [a,b,c,d,e,f] : [a,b,c,d,e,μ,ν,f,ρ,σ]
                F[key] = A[row,col]
            end
        end
    end
    F
end

@doc raw""" 

    numeric_F_symbols(C::SixJCategory; precision=128, convention=:column_major_packing)
    numeric_F_symbols(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision=128, convention=:column_major_packing)

Return [`F_symbols`](@ref) evaluated under the embedding ``e``. The keyword
`convention=:column_major_packing` has the same meaning as for `F_symbols`.
"""
function numeric_F_symbols(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 128, convention::Symbol=:column_major_packing)
    F = F_symbols(C; convention)

    if base_ring(C) == QQ 
        return Dict(k => e(number_field(e)(v), precision) for (k,v) ∈ F)
    else
        return Dict(k => e(v, precision) for (k,v) ∈ F)
    end
end

function numeric_F_symbols(C::SixJCategory; precision = 128, convention::Symbol=:column_major_packing)
    !isdefined(C, :embedding) && error("No embedding has been specified")
    numeric_F_symbols(C, getfield(C, :embedding); precision, convention)
end

@doc raw"""
    R_symbols(C::SixJCategory; convention=:column_major_packing)

Return all R-symbol entries, including zeros. In a multiplicity-free category,
keys are `[a,b,c]` and both conventions agree. Otherwise keys are `[a,b,c,μ,ν]`.

For `:bonderson`, the value is ``B_{\mu\nu}`` defined by
``p^{ba}_{c,\nu}\circ c_{a,b}=\sum_\mu B_{\mu\nu}p^{ab}_{c,\mu}``,
where ``p`` are the fixed projection bases and ``c_{a,b}:a\otimes b\to b\otimes a``
is the braiding. Equivalently, on composition-dual splitting bases,
``c_{a,b}\circ s_{ab}^\mu=\sum_\nu B_{\mu\nu}s_{ba}^\nu``.
Thus the input basis index comes first, as in Bonderson,
Shtengel, and Slingerland, *Interferometry of non-Abelian Anyons*,
Annals of Physics 323 (2008), Eqs. (2.31)--(2.32),
https://doi.org/10.1016/j.aop.2008.01.012.

The default `:column_major_packing` preserves the package's dictionary packing:
the same key has value ``B_{\nu\mu}``. This is a transpose of the two
multiplicity indices, **not** the inverse braiding or a change of gauge.
As with [`F_symbols`](@ref), the plain dictionary does not record its convention.
"""
function R_symbols(C::SixJCategory; convention::Symbol=:column_major_packing)
    _check_symbol_convention(convention)
    if convention == :bonderson
        N = multiplication_table(C)
        mf = multiplicity(C) == 1
        R = Dict{Vector{Int},elem_type(base_ring(C))}()
        for a in 1:rank(C), b in 1:rank(C), c in 1:rank(C)
            B = r_symbol(C,a,b,c)
            for μ in 1:N[a,b,c], ν in 1:N[b,a,c]
                R[mf ? [a,b,c] : [a,b,c,μ,ν]] = B[μ,ν]
            end
        end
        return R
    end

    S = simples(C)

    N = length(S)
    C_morphism_type = morphism_type(C)
    K = base_ring(C) 

    R_dict = Dict{Vector{Int}, elem_type(K)}()

    one_indices = findall(s -> int_dim(Hom(s,one(C))) > 0 , S)
    one_components = simple_subobjects(one(C))

    # Set unitors to identity
    S[one_indices] = one_components

    prods = [X ⊗ Y for X ∈ S, Y ∈ S]


    homs = [basis(Hom(prods[i,j],S[k])) for i ∈ 1:N, j ∈ 1:N, k ∈ 1:N]

    m = multiplicity(C) 
    mult = multiplication_table(C)

    for (a,b,c) ∈ collect(Base.product(1:N, 1:N, 1:N))
        sym = collect(r_symbol(C,a,b,c))[:]

        if mult[a,b,c] == 0 
            continue 
        end

        if m == 1
            R_dict[[a,b,c]] = popfirst!(sym)
        else
            for i ∈ 1:mult[a,b,c]
                for j ∈ 1:mult[b,a,c]
                    R_dict[[a,b,c,i,j]] = popfirst!(sym)
                end
            end
        end
    
    end

    return R_dict           
end

@doc raw""" 

    numeric_R_symbols(C::SixJCategory; precision=2048, convention=:column_major_packing)
    numeric_R_symbols(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision=2048, convention=:column_major_packing)

Return [`R_symbols`](@ref) evaluated under the embedding ``e``. The keyword
`convention=:column_major_packing` has the same meaning as for `R_symbols`.
"""
function numeric_R_symbols(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 2048, convention::Symbol=:column_major_packing)
    R = R_symbols(C; convention)

    Dict(k => e(base_ring(C) == QQ ? number_field(e)(v) : v, precision) for (k,v) ∈ R)
end

function numeric_R_symbols(C::SixJCategory; precision = 2048, convention::Symbol=:column_major_packing)
    !isdefined(C, :embedding) && error("No embedding has been specified")
    numeric_R_symbols(C, getfield(C, :embedding); precision, convention)
end

@doc raw""" 

    P_symbols(C::SixJCategory)

Return a Dictionary of the Pivotal symbols of ``C``.
"""
function P_symbols(C::SixJCategory)
    Dict([k] => C.pivotal[k] for k in 1:rank(C)) 
end

@doc raw""" 

    numeric_P_symbols(C::SixJCategory; precision = 2048)
    numeric_P_symbols(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 2048)

Return a Dictionary of the Pivotal symbols of ``C`` evaluated under the embedding ``e``.
"""
function numeric_P_symbols(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 2048)
    P = P_symbols(C)

    Dict(k => e(v, precision) for (k,v) ∈ P)
end

function numeric_P_symbols(C::SixJCategory; precision = 2048)
    !isdefined(C, :embedding) && error("No embedding has been specified")
    numeric_P_symbols(C, getfield(C, :embedding), precision = precision)
end

@doc raw""" 

    numeric_twists(C::SixJCategory; precision = 2048)
    numeric_twists(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 2048)

Return an array containing the values for the twists of the simples of ``C`` evaluated under the embedding ``e``.
"""
function numeric_twists(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 2048)
    P = twists(C)

    [e(v, precision) for v ∈ P]
end

function numeric_twists(C::SixJCategory; precision = 2048)
    !isdefined(C, :embedding) && error("No embedding has been specified")
    numeric_twists(C, getfield(C, :embedding), precision = precision)
end

@doc raw""" 

    numeric_smatrix(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 2048)

Return the S-matrix of ``C`` evaluated under the embedding ``e``.
"""
function numeric_smatrix(C::SixJCategory, e::AbsSimpleNumFieldEmbedding; precision = 2048)
    S = smatrix(C)

    n = size(S,1)
    M = [e(S[i,j], precision) for i ∈ 1:n, j ∈ 1:n]
end

function numeric_smatrix(C::SixJCategory; precision = 2048)
    !isdefined(C, :embedding) && error("No embedding has been specified")
    numeric_smatrix(C, getfield(C, :embedding), precision = precision)
end

function numeric(C::SixJCategory, precision, max_bits)
    K = base_ring(C)
    CC = AcbField(max_bits)
    if K == QQ 
        _K, = rationals_as_number_field()
        return extension_of_scalars(C, CC, embedding = x -> C.embedding(_K(x), max_bits))
    end
    if !(typeof(K) <: Union{QQBarField, ArbField, AcbField})
        return extension_of_scalars(C, CC, embedding = x -> CC(C.embedding(x, max_bits)))
    end

    extension_of_scalars(C, CC, embedding = x -> qqbar_to_acb_with_error(x,precision,max_bits))
end

#=----------------------------------------------------------
    save and load
----------------------------------------------------------=#
#@register_serialization_type SixJCategory "SixJCategory"

# function save_object(s::SerializerState, C::SixJCategory)
#     fields = collect(fieldnames(SixJCategory)[1:end-1])
#     filter!(e -> isdefined(C, e), fields)

#     # setup serialization. Every vector has to be a Tuple
#     data= NamedTuple(Dict(e => getfield(C, e) isa Array ? Tuple(getfield(C, e)) : getfield(C, e) for e ∈ fields))
#     @show typeof(data)
#     save_object(s, data)
# end

# function load_object(s::DeserializerState, ::Type{SixJCategory})
#     data = load_object(s, @NamedTuple{simples::Int64, base_ring::AbsSimpleNumField, tensor_product::NTuple{27, Int64}, ass::NTuple{81, MatSpaceElem{AbsSimpleNumFieldElem}}, simples_names::Tuple{String, String, String}, pivotal::Tuple{AbsSimpleNumFieldElem, AbsSimpleNumFieldElem, AbsSimpleNumFieldElem}, name::String, one::Tuple{Int64, Int64, Int64}}
# )
#     C = six_j_category(data.base_ring, data.tensor_product, collect(data.simples_names))
#     for (k,v) ∈ pairs(data)
#         if k != :base_ring && k != :tensor_product && k != :simples_names
#             setfield!(C, k, v isa Tuple ? collect(v) : v)
#         end
#     end
#     return C
# end
