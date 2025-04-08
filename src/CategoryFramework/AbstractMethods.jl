#=----------------------------------------------------------
    A B S T R A C T   M E T H O D S
    
    Here are all generic methods dealing with at most the 
    structure of 
    
    - additive
    - (semisimple) abelian
    - 𝕜-linear

    categories.
----------------------------------------------------------=#

#=----------------------------------------------------------
    Endomorphism Ring
----------------------------------------------------------=#

""" 

endomorphism_ring(X::Object)

Return the endomorphism ring of ``X`` as a matrix algebra.
"""
function endomorphism_ring(X::Object, base = basis(End(X)))
    @assert is_abelian(parent(X))
    if length(base) == 0
        return matrix_algebra(base_ring(X), [zero_matrix(base_ring(X),1,1)], isbasis = true)
    end

    if hasmethod(matrix, Tuple{typeof(base[1])})
        try 
            return endomorphism_ring_by_matrices(X,base)
        catch
            @warn "matrix approach failed"
        end
    end
    
    return endomorphism_ring_by_basis(X,base)
end

function endomorphism_ring(X::Object, H::AbstractHomSpace)
    @assert domain(H) == codomain(H) == X
    endomorphism_ring(X, basis(H))
end

function extension_of_scalars(R::AbstractAssociativeAlgebra, F::Ring)
    A = StructureConstantAlgebra(F, F.(multiplication_table(R)), F.(coefficients(one(R))))
end

⊗(F::Ring, R::AbstractAssociativeAlgebra) = extension_of_scalars(R,F)
⊗(R::AbstractAssociativeAlgebra, F::Ring) = extension_of_scalars(R,F)

function endomorphism_ring_by_matrices(X::Object, base = basis(End(X)))

    try 
        mats = matrix.(base)
        matrix_algebra(base_ring(mats[1]), mats, isbasis = true)
    catch
        error("Non-matrix support comming")
    end
end

function endomorphism_ring_by_basis(X::Object, base = basis(End(X)))
    n = length(base)
    mult = Array{elem_type(base_ring(X)),3}(undef,n,n,n)

    for i ∈ 1:n, j ∈ 1:n
        mult[i,j,:] = express_in_basis(base[i] ∘ base[j], base)
    end

    A = StructureConstantAlgebra(base_ring(X), mult)
    A.one = express_in_basis(id(X), base)
    return A
end

function is_invertible(f::Morphism)
    try
        inv(f) 
        return true
    catch
        return false
    end
end
#------------------------------------------------------
#   Image & Kernel
#------------------------------------------------------

@doc raw""" 

    image(f::Morphism)

Return the image ``Im(f)`` of ``f:X → Y`` together with a monomorphism 
``Im(f) ↪ Y``.
"""
function image(f::Morphism)
    C,c = cokernel(f)
    return kernel(c)
end

#=----------------------------------------------------------
    Hom Space Functionality 
----------------------------------------------------------=#
@doc raw""" 

    express_in_basis(f::Morphism, B::Vector{Morphism})

Return a vector of coefficients expressing ``f``in the basis ``B``.
"""
function express_in_basis(f::T, B::Vector{T}) where T <: Morphism
    F = base_ring(f)    
    # b_mat = hcat([[x for x ∈ matrix(b)][:] for b ∈ B]...)
    # B_mat = matrix(F, size(b_mat,1), size(b_mat, 2), b_mat)
    # f_mat = matrix(F, *(size(matrix(f))...), 1, [x for x ∈ matrix(f)][:])
    vec_f = morphism(matrix(f))
    vec_basis = VSMorphism[morphism(matrix(b)) for b ∈ B]
    return express_in_basis(vec_f, vec_basis)
end

@doc raw""" 

    express_in_basis(f::Morphism)

Return a vector of coefficients expressing ``f: X → Y`` in the basis o f ``\mathrm{Hom}(X,Y)``.
"""
function express_in_basis(f::Morphism, H::AbstractHomSpace = Hom(domain(f), codomain(f)))
    express_in_basis(f, basis(H))
end

@doc raw""" 

    is_simple(X::Object)

Check whether ``X`` is a simple object.  
"""
function is_simple(X::Object)
    if hasfield(typeof(X), :__attrs)
        get_attribute!(X, :is_simple) do
            B = basis(End(X))
            if length(B) == 0
                return false
            elseif length(B) == 1
                return true
            end
            
            try 
                is_simple(endomorphism_ring(X, B))
            catch
                false
            end
        end
    end

    B = basis(End(X))
    if length(B) == 0
        return false
    elseif length(B) == 1
        return true
    end
    
    try 
        is_simple(endomorphism_ring(X, B))
    catch
        false
    end
end

function is_indecomposable(X::Object)
    dec = decompose(X)
    length(dec) == 1 && dec[1][2] == 1
end

@doc raw""" 

    decompose(X::Object)

Decompose an object ``X`` in an abelian category.
"""
decompose(X::Object) = decompose_by_endomorphism_ring(X)

@doc raw""" 

    decompose(X::Object, S::Vector{Object})

Decompose an object ``X`` in a semisimple category into simple objects of ``S``.
"""
decompose(X::T, S::Vector{T}) where T <: Object = decompose_by_simples(X,S)

# function decompose_by_endomorphism_ring(X::Object)
#     base = basis(End(X))
#     end_ring = endomorphism_ring(X, base)

#     idems = central_primitive_idempotents(end_ring)
#     idems_coefficients = coefficients.(idems)

#     cat_idems = [sum(c .* base) for c ∈ idems_coefficients]
#     subobjects = [image(f)[1] for f ∈ cat_idems]

#     uniques = unique_indecomposables(subobjects)
#     [(s, length(findall(r -> is_isomorphic(s,r)[1], subobjects))) for s ∈ uniques]
# end

function decompose_by_endomorphism_ring(X::Object, E = End(X))
    K = base_ring(X)
    if K == QQBarField() || typeof(K) == CalciumField
        return decompose_over_qqbar(X, E)
    end

    if int_dim(E) == 0 return Tuple{object_type(parent(X)), Int}[] end

    idems = central_primitive_idempotents(E)

    _images = [image(f)[1] for f ∈ idems]

    # Check for matrix algebras
    images = hcat([simple_subobjects(i, End(i), true) for i ∈ _images]...)
    

    tuples = Tuple{typeof(X), Int}[]

    for Y ∈ images
        i = findfirst(Z -> is_isomorphic(Z[1],Y)[1], tuples)
        k = div(int_dim(Hom(Y,X)), int_dim(End(Y)))
        if i === nothing
            push!(tuples, (Y,k))
        else
            tuples[i] = (Y, tuples[i][2] + k)
        end
    end

    return tuples
end


function decompose_over_qqbar(X::Object, E = End(X))
    #@assert is_semisimple(parent(X))
    if int_dim(E) == 1
        try 
            set_attribute!(X, :is_simple, true)
        catch 
        end
        return [(X,1)]
    end

    f = [f for f ∈ basis(E) if degree(minpoly(f)) > 1][1]
    
    K = collect(values(eigenvalues(f)))
    subs = vcat([decompose_over_qqbar(k) for k ∈ K]...)
    subs = typeof(X)[s for (s,_) ∈ subs]
    unique_simples_with_multiplicity(subs)
end

function simple_subobjects_over_qqbar(X::Object, E = End(X))
    return [s for (s,_) ∈ decompose_over_qqbar(X,E)]
end


function decompose_by_simples(X::Object, S = simples(parent(X)))
    C = parent(X)
    dimensions = [div(int_dim(Hom(s,X)), int_dim(End(s))) for s ∈ S]
    return [(s,d) for (s,d) ∈ zip(S,dimensions) if d > 0]
end


function direct_sum_decomposition(X::Object, S = simples(parent(X)))
    C = parent(X)
    @assert is_semisimple(C) "Semisimplicity required"
    
    if X == zero(C) return X, id(X), [], [] end

    components = decompose(X,S)
    Z, incl, proj = direct_sum(vcat([[s for _ ∈ 1:d] for (s,d) ∈ components]...)...)

    # temporary solution!
    iso = is_isomorphic(X,Z)[2]
    
    inv_iso = inv(iso)

    return Z, iso, [inv_iso∘i for i ∈ incl], [p∘iso for p ∈ proj]

    #----------------------------------
    f = zero_morphism(X,Z)

    for (p,i) ∈ zip(proj, incl)
        g = i∘p
        f = f + g
    end
    return Z, f, incl, proj
end

@doc raw""" 

    central_primitive_idempotents(H::AbstractHomSpace)

Compute the central primitive idempotents of an endomorphism space ``H``.
"""
function central_primitive_idempotents(H::AbstractHomSpace)
    @assert domain(H) == codomain(H) "Not an endomorphism algebra"

    base = basis(H)
    A = endomorphism_ring(H.X, base)
    one(A)
 

    if !is_semisimple(A) && characteristic(base_ring(H)) == 0 
        X = semisimplify(domain(H))
        base = basis(semisimplify(H))
        A = endomorphism_ring(X, base)
        base = morphism.(base)
    end

    A.issemisimple = true

    idems = central_primitive_idempotents(A)
    [sum(base .* coefficients(i)) for i ∈ idems]
end

function radical(H::AbstractHomSpace)
    @assert domain(H) == codomain(H)

    A = endomorphism_ring(H.X, basis(H))
    R = radical(A)
    [sum(basis(H) .* coefficients(r)) for r ∈ basis(R)]
end

@doc raw""" 

    gens(H::AbstractHomSpace)

Compute generators of an endomorphism space ``H``.
"""
function gens(H::AbstractHomSpace)
    @assert H.X == H.Y "Not an endomorphism algebra"

    A = endomorphism_ring(H.X, basis(H))
    A.issemisimple = true
    gs = gens(A)
    [sum(basis(H) .* coefficients(i)) for i ∈ gs]
end

function is_subobject(X::Object, Y::Object)
    @assert parent(X) == parent(Y)
    @assert is_semisimple(parent(X))
    
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

function is_isomorphic(X::Object, Y::Object)
    EX = End(X)
    EY = End(Y)
    H = Hom(X,Y) 

    mats = [express_in_basis]

    @error "Not implemented"
end

#-------------------------------------------------------
# Semisimple: Subobjects
#-------------------------------------------------------

@doc raw""" 

    eigenvalues(f::Morphism)

Compute the eigenvalues of ``f``. Return a dictonary with 
entries `λ => ker(f - λid)`.
"""
function eigenvalues(f::Morphism)
    @assert domain(f) == codomain(f) "Not an endomorphism"

    #@show factor(minpoly(matrix(f)))

    vals = eigenvalues(matrix(f))


    return Dict(λ => kernel(f-λ*id(domain(f)))[1] for λ ∈ vals)
end




# function indecomposable_subobjects_by_matrix_algebra(X::Object, E = End(X))
#     if length(basis(E)) == 0
#         return typeof(X)[]
#     end
#     A = endomorphism_ring(X, basis(E))
#     dec = decompose(A)
#     if length(dec) == 1
#         return [X]
#     end

#     s,f = dec[1]

#     b = sum(coefficients(image(f,basis(s)[1])) .* basis(E))

#     eig_spaces = eigenvalues(b)
#     λ,_ = collect(eig_spaces)[1]
#     K,i = kernel(b - λ*id(X))
#     C,_ = cokernel(i) 

#     return unique_simples([indecomposable_subobjects(K); indecomposable_subobjects(C)])
# end

function indecomposable_subobjects(X::Object)
    [x for (x,_) ∈ decompose(X)]
end

@doc raw""" 

    minpoly(f::Morphism)

Compute the minimal polynomial of an endomorphism ``f: X \to X``.
"""
function minpoly(f::Morphism)
    @assert domain(f) == codomain(f) "Not an edomorphism"
    
    if hasmethod(matrix, Tuple{typeof(f)})
        # if typeof(base_ring(f)) == CalciumField
        #     return change_base_ring(base_ring(f), minpoly(change_base_ring(QQBarField(), matrix(f))))
        # end
        return minpoly(matrix(f))
    else
        error("Generic minpoly coming soon")
    end
end 

function subst(f::PolyRingElem{T}, m::Morphism) where T <: FieldElem
    @assert domain(m) == codomain(m)

    R = parent(f)

    return sum([c * compose([id(domain(m)); [m for _ ∈ 1:d]]...) for (c,d) ∈ zip(coefficients(f), 0:degree(f))])
end

function _decompose_by_simple_endomorphism_ring(X::Object, E = End(X))

    K = base_ring(X)
    R = endomorphism_ring(X, E)
    CR,_ = center(R)
    dR = dim(R)

    if !is_square(div(dR,dim(CR))) || dR == 1 || is_commutative(R) 
        return (X,1)
    end
        
    G = basis(E)

    n,_ = size(matrix(G[1]))
    mats = matrix.(G)

    for f ∈ E
        eig = eigenvalues(f)

        length(eig) < 1  && continue

        Z = collect(values(eig))[1]

        length(eig) == 1 && dim(Z) == dim(X) && continue
                
        EZ = End(Z)

        k = sqrt(int_dim(E)) - sqrt(int_dim(EZ))

        (Z2, k2) = _decompose_by_simple_endomorphism_ring(Z,EZ)

        return (Z2, k*k2)
    end     

    for i ∈ 1:n
        m = matrix(K,n,length(G), hcat([M[i,:] for M ∈ mats]...))

        d,N = nullspace(m)

        if d > 0 
            Y,y = kernel(sum(collect(N[:,1]) .* G))
            Z,k = _decompose_by_simple_endomorphism_ring(Y, End(Y))
            return (Z, sqrt(dR // int_dim(End(Z))))
        end
    end

    return (X,1)
end


# function indecomposable_subobjects(X::Object, E = End(X))
#     _indecomposable_subobjects(X,E)
# end


rank(C::Category) = length(simples(C))

is_isomorphic_simples(X::Object, Y::Object) = is_isomorphic(X,Y)


function unique_simples(simples::Vector{<:Object})
    unique_simples = simples[1:1]
    for s ∈ simples[2:end]
        if sum([is_isomorphic_simples(s,u)[1] for u ∈ unique_simples]) == 0
            unique_simples = [unique_simples; s]
        end
    end

    return unique_simples
end

function unique_simples_with_multiplicity(simples::Vector{<:Object})
    unique_simples = [(simples[1], 1)]
    for s ∈ simples[2:end]
        i = findfirst(u -> is_isomorphic_simples(s,u[1])[1], unique_simples)
        if i === nothing
            unique_simples = [unique_simples; (s, 1)]
        else
            unique_simples[i] = (unique_simples[i][1], unique_simples[i][2] + 1)
        end
    end
    return unique_simples
end

function unique_indecomposables(simples::Vector{<:Object})
    if is_semisimple(parent(simples[1]))
        return unique_simples(simples)
    end
    uniques = simples[1:1]
    for s ∈ simples[2:end]
        if *([!is_isomorphic(s,u)[1] for u ∈ uniques]...)
            uniques = [uniques; s]
        end
    end
    return uniques
end

function simples_names(C::Category) 
    @assert is_semisimple(C)
    return ["X$i" for i ∈ 1:length(simples(C))]
end

function indecomposables_names(C::Category)
    try 
        return simples_names(C)
    catch
        return ["X$i" for i ∈ 1:length(indecomposables(C))]
    end
end

function indecomposables(C::Category)
    if is_semisimple(C)
        return simples(C)
    end
    error("Cannot compute indecomposables")
end

function simples(C::Category)
    simpls = indecomposables(C)
    if is_semisimple(C)
        return simpls
    end
    return [s for s ∈ simpls if is_simple(s)]
end


@doc raw""" 

    left_inverse(f::Morphism)

Compute a morphism ``g`` such that ``g ∘ f = id``. Errors if 
``f``is not mono. 
"""
function left_inverse(f::Morphism)
    X = domain(f)
    Y = codomain(f)

    if X == zero(parent(X))
        return zero_morphism(Y,X)
    end

    HomYX = basis(Hom(Y,X))
    base = basis(End(X))

    K = base_ring(f)
    Kx,x = polynomial_ring(K, length(HomYX))
    eqs = [zero(Kx) for _ ∈ length(base)]

    for (g,y) ∈ zip(HomYX,x)
        eqs = eqs .+ (y.*express_in_basis(g∘f, base))
    end

    one_coeffs = express_in_basis(id(X), base)

    M = zero_matrix(K,length(x),length(eqs))
    for (i,e) ∈ zip(1:length(eqs), eqs)
        M[:,i] = [coeff(e,y) for y ∈ x]
    end
  
    try 
         N = solve(M, matrix(K,1,length(one_coeffs), one_coeffs), side = :left)
        return sum(HomYX .* collect(N[1,:])[:])
    catch e
        error("Morphism does not have a left inverse")
    end
end

@doc raw""" 

    right_inverse(f::Morphism)

Compute a morphism ``g`` such that ``f ∘ g = id``. Errors if ``f``
is not epi.
"""
function right_inverse(f::Morphism)
    X = domain(f)
    Y = codomain(f)

    if Y == zero(parent(X))
        return zero_morphism(Y,X)
    end

    HomYX = basis(Hom(Y,X))
    base = basis(End(Y))

    K = base_ring(f)
    Kx,x = polynomial_ring(K, length(HomYX))
    eqs = [zero(Kx) for _ ∈ length(base)]

    for (g,y) ∈ zip(HomYX,x)
        eqs = eqs .+ (y.*express_in_basis(f∘g, base))
    end
    
    one_coeffs = express_in_basis(id(Y), base)

    M = zero_matrix(K,length(x), length(eqs))
    for (i,e) ∈ zip(1:length(eqs), eqs)
        M[:,i] = [coeff(e,y) for y ∈ x]
    end

    try 
        N = solve(M, matrix(K,1,length(one_coeffs), one_coeffs), side = :left)
        return sum(HomYX .* collect(N[1,:])[:])
    catch e
        error("Morphism does not have a right inverse")
    end
end

function has_right_inverse(f::Morphism)
    try 
        right_inverse(f)
        true
    catch
        false
    end
end

@doc raw""" 

    zero_morphism(X::Object, Y::Object)

Compute the zero morphism between ``X``and ``Y``.
"""
function zero_morphism(X::Object, Y::Object) 
    C = parent(X)
    if is_additive(C)
        return basis(Hom(zero(parent(X)), Y))[1] ∘ basis(Hom(X, zero(parent(X))))
    elseif is_linear(C) 
        return zero(base_ring(X)) * basis(Hom(X,Y))[1]
    end
    @error "There might be no zero morphism"
end

@doc raw""" 

    zero_morphism(X::Object)

Compute the zero morphism on ``X``.
"""
zero_morphism(X::Object) = zero_morphism(X,X)

zero_morphism(C::Category) = zero_morphism(zero(C))

is_zero(f::Morphism) = f == zero_morphism(domain(f), codomain(f))

is_zero(X::Object) = int_dim(End(X)) == 0

function ==(f::Morphism, x::T) where T <: Union{RingElem, Int} 
    if x == 0 
        return is_zero(f)
    else
        if domain(f) != codomain(f)
            return false
        end
        try
            express_in_basis(f, [id(domain(f))]) == [x]
        catch end
    end
    false
end

@doc raw""" 

    is_monomoprhism(f::Morphism)

Check whether ``f`` mono.
"""
function is_monomorphism(f::Morphism)
    try 
        left_inverse(f)
        return true
    catch
        false
    end
end

@doc raw""" 

    is_epimorphism(f::Morphism)

Check wether ``f``is epi.
"""
function is_epimorphism(f::Morphism)
    try 
        right_inverse(f)
        return true
    catch
        false
    end
end


#---------------------------------------------------------
#   Horizontal and Vertical direct sums
#---------------------------------------------------------

"""
    function horizontal_direct_sum(f::Morphism, g::Morphism)

Return the sum of ``f:X → Z``, ``g:Y → Z`` as ``f+g:X⊕Y → Z.
"""
function horizontal_direct_sum(f::Morphism, g::Morphism)
    #@assert codomain(f) == codomain(g) "Codomains do not coincide"
    sum = f ⊕ g
    _,_,(p1,p2) = direct_sum(codomain(f),codomain(g))
    return p1∘sum + p2∘sum
end

function horizontal_direct_sum(f::Vector{M}) where M <: Morphism
    #@assert codomain(f) == codomain(g) "Codomains do not coincide"
    f_sum = direct_sum(f...)
    _,_,p = direct_sum([codomain(fi) for fi ∈ f]...)
    return sum([p1∘f_sum for p1 ∈ p])
end

"""
    function vertical_direct_sum(f::Morphism, g::Morphism)

Return the sum of ``f:X → Y``, ``g:X → Z`` as ``f+g: X → Y⊕Z.
"""
function vertical_direct_sum(f::Morphism, g::Morphism)
    #@assert domain(f) == domain(g) "Domains do not coincide"

    sum = f ⊕ g
    _,(i1,i2),_ = direct_sum(domain(f), domain(g))
    return sum∘i1 + sum∘i2
end

function vertical_direct_sum(f::Vector{M}) where M <: Morphism
    f_sum = direct_sum(f...)
    _,i,_ = direct_sum([domain(fi) for fi ∈ f]...)
    return sum([f_sum∘ix for ix ∈ i])

end

