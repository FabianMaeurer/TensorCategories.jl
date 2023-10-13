function zero_morphism(X::Object, Y::Object) 
    C = parent(X)
    if is_additive(C)
        return basis(Hom(zero(parent(X)), Y))[1] ∘ basis(Hom(X, zero(parent(X))))
    elseif is_linear(C) 
        return zero(base_ring(X)) * basis(Hom(X,Y))[1]
    end
    @error "There might be no zero morphism"
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
    try 
        mats = matrix.(base)
        matrix_algebra(base_ring(X), mats, isbasis = true)
    catch
        error("Non-matrix support comming")
    end
end

#------------------------------------------------------
#   Image & Kernel
#------------------------------------------------------


function image(f::Morphism)
    C,c = cokernel(f)
    return kernel(c)
end

#=----------------------------------------------------------
    Hom Space Functionality 
----------------------------------------------------------=#
function express_in_basis(f::T, B::Vector{T}) where T <: Morphism
    F = base_ring(f)    
    b_mat = hcat([[x for x ∈ matrix(b)][:] for b ∈ B]...)
    B_mat = matrix(F, size(b_mat,1), size(b_mat, 2), b_mat)
    f_mat = matrix(F, 1, *(size(matrix(f))...), [x for x ∈ matrix(f)][:])

    return [x for x ∈ solve_left(transpose(B_mat),f_mat)][:]
end

function is_simple(X::Object)
    B = basis(End(X))
    if length(B) == 0
        return false
    end
    is_simple(endomorphism_ring(X, B))
end

decompose(X::Object) = decompose_by_endomorphism_ring(X)

decompose(X::T, S::Vector{T}) where T <: Object = decompose_by_simples(X,S)

function decompose_by_endomorphism_ring(X::Object)
    base = basis(End(X))
    end_ring = endomorphism_ring(X, base)

    idems = central_primitive_idempotents(end_ring)
    idems_coefficients = coefficients.(idems)

    cat_idems = [sum(c .* base) for c ∈ idems_coefficients]
    subobjects = [image(f)[1] for f ∈ cat_idems]

    uniques = unique_indecomposables(subobjects)
    [(s, length(findall(r -> is_isomorphic(s,r)[1], subobjects))) for s ∈ uniques]
end

function decompose_by_simples(X::Object, S = simples(parent(X)))
    C = parent(X)
    dimensions = [div(int_dim(Hom(s,X)), int_dim(End(s))) for s ∈ S]
    return [(s,d) for (s,d) ∈ zip(S,dimensions) if d > 0]
end

function direct_sum_decomposition(X::Object, S = simples(parent(X)))
    C = parent(X)
    @assert is_semisimple(C) "Semisimplicity required"
    
    if X == zero(C) return id(X), [], [] end

    components = decompose(X,S)
    Z, incl, proj = direct_sum(vcat([[s for _ ∈ 1:d] for (s,d) ∈ components]...)...)

    # temporary solution!
    iso = is_isomorphic(X,Z)[2]
    return Z, iso, [inv(iso)∘i for i ∈ incl], [p∘iso for p ∈ proj]

    #----------------------------------
    f = zero_morphism(X,Z)

    for (p,i) ∈ zip(proj, incl)
        g = i∘p
        f = f + g
    end
    return Z, f, incl, proj
end

function central_primitive_idempotents(H::AbstractHomSpace)
    @assert H.X == H.Y "Not an endomorphism algebra"

    A = endomorphism_ring(H.X, basis(H))
    A.issemisimple = true
    idems = central_primitive_idempotents(A)
    [sum(basis(H) .* coefficients(i)) for i ∈ idems]
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

#-------------------------------------------------------
# Semisimple: Subobjects
#-------------------------------------------------------

function eigenvalues(f::Morphism)
    @assert domain(f) == codomain(f) "Not an endomorphism"

    #@show factor(minpoly(matrix(f)))
    if base_ring(f) == QQBar
        vals = eigenvalues(matrix(f))
    else
        vals = keys(spectrum(matrix(f)))
    end

    return Dict(λ => kernel(f-λ*id(domain(f)))[1] for λ ∈ vals)
end

function indecomposable_subobjects_by_matrix_algebra(X::Object, E = End(X))
    if length(basis(E)) == 0
        return typeof(X)[]
    end
    A = endomorphism_ring(X, basis(E))
    dec = decompose(A)
    if length(dec) == 1
        return [X]
    end

    s,f = dec[1]

    b = sum(coefficients(image(f,basis(s)[1])) .* basis(E))

    eig_spaces = eigenvalues(b)
    λ,_ = collect(eig_spaces)[1]
    K,i = kernel(b - λ*id(X))
    C,_ = cokernel(i) 

    return unique_simples([indecomposable_subobjects(K); indecomposable_subobjects(C)])
end

function indecomposable_subobjects(X::Object, E = End(X))
    @assert is_semisimple(parent(X)) "Non semisimple categories are not yet supported"
    B = basis(E)
   
    if length(B) == 1 return [X] end

    for f ∈ B
        eig_spaces = eigenvalues(f)
        if length(eig_spaces) == 0 
            continue
            
        elseif length(eig_spaces) == 1 && dim(collect(values(eig_spaces))[1]) == dim(X)
            continue
        end

        λ = collect(keys(eig_spaces))[1]
        K,i = kernel(f - λ*id(X))
        C,_ = cokernel(i)

        return unique_indecomposables([indecomposable_subobjects(K); indecomposable_subobjects(C)])
    end

    if is_simple(X) 
        return [X]
    end
    indecomposable_subobjects_by_matrix_algebra(X,E)    
end


# function indecomposable_subobjects(X::Object, E = End(X))
#     _indecomposable_subobjects(X,E)
# end

function simple_subobjects(X::Object, E = End(X))
    indecomposables = indecomposable_subobjects(X, E)
    if is_semisimple(parent(X))
        return indecomposables
    else
        return indecomposables[is_simple.(indecomposables)]
    end
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

function is_simple(X::Object, S = simples(parent(X)))
    for s ∈ S
        if int_dim(Hom(s,X)) != 0 
            return is_isomorphic(s,X)[1]
        end
    end
    error("You might miss some simples")
end


function left_inverse(f::Morphism)
    X = domain(f)
    Y = codomain(f)

    if X == zero(parent(X))
        return zero_morphism(Y,X)
    end

    HomYX = basis(Hom(Y,X))
    base = basis(End(X))

    K = base_ring(f)
    Kx,x = PolynomialRing(K, length(HomYX))
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
         N = solve_left(M, matrix(K,1,length(one_coeffs), one_coeffs))
        return sum(HomYX .* collect(N[1,:])[:])
    catch e
        error("Morphism does not have a left inverse")
    end
end


function right_inverse(f::Morphism)
    X = domain(f)
    Y = codomain(f)

    if Y == zero(parent(X))
        return zero_morphism(Y,X)
    end

    HomYX = basis(Hom(Y,X))
    base = basis(End(Y))

    K = base_ring(f)
    Kx,x = PolynomialRing(K, length(HomYX))
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
        N = solve_left(M, matrix(K,1,length(one_coeffs), one_coeffs))
        return sum(HomYX .* collect(N[1,:])[:])
    catch e
        error("Morphism does not have a right inverse")
    end

end