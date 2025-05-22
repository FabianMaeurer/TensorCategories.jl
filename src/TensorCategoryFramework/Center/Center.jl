@attributes mutable struct CenterCategory <: Category
    base_ring::Ring
    category::Category
    simples::Vector{O} where O <: Object
    inductions::Dict{<:Object,<:Object}
    induction_gens::Vector{Object}

    function CenterCategory(F::Ring, C::Category)
        Z = new()
        Z.base_ring = F
        Z.category = C
        return Z
    end

    function CenterCategory()
        new()
    end
end

@attributes mutable struct CenterObject <: Object
    parent::CenterCategory
    object::Object
    γ::Vector{M} where M <: Morphism

    function CenterObject(parent::CenterCategory, object::Object,
        γ::Vector{<: Morphism})
        new(parent, object, γ)
    end
end

struct CenterMorphism <: Morphism
    domain::CenterObject
    codomain::CenterObject
    m::Morphism
end

function ==(C::CenterCategory, D::CenterCategory)
    if !isdefined(C, :simples) || !isdefined(D, :simples)
        if !isdefined(C, :simples) ⊻ !isdefined(D, :simples)
            return false
        else
            return base_ring(C) == base_ring(D) && C.category == D.category
        end
    elseif length(C.simples) != length(D.simples)
        return false
    end
    return base_ring(C) == base_ring(D) && C.category == D.category && *([isequal_without_parent(s,t) for (s,t) ∈ zip(C.simples, D.simples)]...)
end

function isequal_without_parent(X::CenterObject, Y::CenterObject)
    return object(X) == object(Y) && half_braiding(X) == half_braiding(Y)
end

is_multifusion(C::CenterCategory) = is_multifusion(category(C))

is_modular(C::CenterCategory) = true 
is_braided(C::CenterCategory) = true
is_rigid(C::CenterCategory) = true

function induction_generators(C::CenterCategory) 
    if isdefined(C, :induction_gens)
        return C.induction_gens
    end
    simpls = simples(category(C))
    invertibls = invertibles(category(C))

    ind_gens = object_type(category(C))[]

    indices = collect(eachindex(simpls))
    while !isempty(indices)
        i = popfirst!(indices)
        x = simpls[i]
        push!(ind_gens, x)
        indices = [s for s ∈ indices if all([!is_isomorphic(simpls[s], i ⊗ x ⊗ dual(i))[1] for i ∈ invertibls])]
    end
    C.induction_gens = ind_gens
end
#-------------------------------------------------------------------------------
#   Center Constructor
#-------------------------------------------------------------------------------
"""
    center(C::Category)

Return the Drinfeld center of ```C```.
"""
function center(C::Category; equivalence = false)
    #@assert is_semisimple(C) "Semisimplicity required"
    return CenterCategory(base_ring(C),C)
end

function morphism(dom::CenterObject, cod::CenterObject, m::Morphism)
    return CenterMorphism(dom,cod,m)
end

"""
    half_braiding(Z::CenterObject)

Return  a vector with half braiding morphisms ```Z⊗S → S⊗Z``` for all simple
objects ```S```.
"""
half_braiding(Z::CenterObject) = Z.γ


"""
    object(X::CenterObject)

Return the image under the forgetful functor.
"""
object(X::CenterObject) = X.object

@doc raw""" 

    morphism(f::CenterMorphism)

Return the image under the forgetful functor.
"""
morphism(f::CenterMorphism) = f.m

is_weakly_fusion(C::CenterCategory) = dim(category(C)) != 0
is_fusion(C::CenterCategory) = all([int_dim(End(s)) == 1 for s ∈ simples(C)])
is_abelian(C::CenterCategory) = true
is_linear(C::CenterCategory) = true
is_monoidal(C::CenterCategory) = true
is_spherical(C::CenterCategory) = is_spherical(category(C))

"""
    add_simple!(C::CenterCategory, S::CenterObject)

Add the simple object ```S``` to the vector of simple objects.
"""
function add_simple!(C::CenterCategory, S::CenterObject)
    @assert dim(End(S)) == 1 "Not simple"
    if isdefined(C, :simples)
        C.simples = filter(e -> e != zero(C), unique_simples([simples(C); S]))
    else
        C.simples = filter(e -> e != zero(C), unique_simples([S]))
    end
end


function add_simple!(C::CenterCategory, S::Array{CenterObject})
    @assert prod(dim(End(s)) for s ∈ S) == 1 "Not simple"
    if isdefined(C, :simples)
        C.simples = unique_simples([simples(C); S])
    else
        C.simples = unique_simples(S)
    end
end

"""
    spherical(X::CenterObject)

Return the spherical structure ```X → X∗∗``` of ```X```.
"""
spherical(X::CenterObject) = morphism(X,dual(dual(X)), spherical(X.object))

(F::Field)(f::CenterMorphism) = F(f.m)

#=-------------------------------------------------
    MISC 
-------------------------------------------------=#

==(f::CenterMorphism, g::CenterMorphism) = f.m == g.m

#-------------------------------------------------------------------------------
#   Direct Sum & Tensor Product
#-------------------------------------------------------------------------------

"""
    direct_sum(X::CenterObject, Y::CenterObject)

Return the direct sum object of ```X``` and ```Y```.
"""
function direct_sum(X::CenterObject, Y::CenterObject)
    S = simples(parent(X.object))
    Z,(ix,iy),(px,py) = direct_sum(X.object, Y.object)

    γZ = [(id(S[i])⊗ix)∘(X.γ[i])∘(px⊗id(S[i])) + (id(S[i])⊗iy)∘(Y.γ[i])∘(py⊗id(S[i])) for i ∈ 1:length(S)]

    CZ = CenterObject(parent(X), Z, γZ)
    ix,iy = CenterMorphism(X,CZ,ix), CenterMorphism(Y,CZ, iy)
    px,py = CenterMorphism(CZ,X,px), CenterMorphism(CZ,Y,py)
    return CZ,[ix,iy],[px,py]
end



"""
    direct_sum(f::CenterMorphism, g::CenterMorphism)

Return the direct sum of ```f``` and ```g```.
"""
function direct_sum(f::CenterMorphism, g::CenterMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    m = f.m ⊕ g.m
    return morphism(dom,cod, m)
end

"""
    tensor_product(X::CenterObject, Y::CenterObject)

Return the tensor product of ```X``` and ```Y```.
"""
function tensor_product(X::CenterObject, Y::CenterObject)
    Z = X.object ⊗ Y.object
    γ = Morphism[]
    simple_objects = simples(parent(X.object))

    x,y = X.object, Y.object

    for (S, yX, yY) ∈ zip(simple_objects, half_braiding(X), half_braiding(Y))

        half_braiding_with_S = associator(S,x,y) ∘ 
                                (yX⊗id(y)) ∘
                                inv_associator(x,S,y) ∘ 
                                (id(x)⊗yY) ∘ 
                                associator(x,y,S)
                                
        push!(γ, half_braiding_with_S)
    end
    return CenterObject(parent(X), Z, γ)
end


"""
    tensor_product(f::CenterMorphism,g::CenterMorphism)

Return the tensor product of ```f``` and ```g```.
"""
function tensor_product(f::CenterMorphism,g::CenterMorphism)
    dom = domain(f)⊗domain(g)
    cod = codomain(f)⊗codomain(g)
    return morphism(dom,cod,f.m⊗g.m)
end

"""
    zero(C::CenterCategory)

Return the zero object of ```C```.
"""
function zero(C::CenterCategory)
    Z = zero(C.category)
    CenterObject(C,Z,[zero_morphism(Z,Z) for _ ∈ simples(C.category)])
end

"""
    one(C::CenterCategory)

Return the one object of ```C```.
"""
function one(C::CenterCategory)
    Z = one(C.category)
    CenterObject(C,Z,[id(s) for s ∈ simples(C.category)])
end


#-------------------------------------------------------------------------------
#   Is central?
#-------------------------------------------------------------------------------

"""
    is_central(Z::Object)

Return true if ```Z``` is in the categorical center, i.e. there exists a half-braiding on ```Z```.
"""
function is_central(Z::Object, simples::Vector{<:Object} = simples(parent(Z)))
    if prod([is_isomorphic(Z⊗s,s⊗Z)[1] for s ∈ simples]) == 0
        return false
    end
    return dim(build_center_ideal(Z,simples)) >= 0
end

function build_natural_center_ideal(Z::Object, indecs = indecomposables(parent(Z)))
    @assert is_additive(parent(Z))

    # Compute a basis for the natural transformations
    nat_trans = additive_natural_transformations(Z⊗-, (-)⊗Z, indecs)

    K = base_ring(Z)

    Kx,x = polynomial_ring(K, length(nat_trans))

    eqs = []

    i_O = findfirst(e -> is_isomorphic(e,one(parent(Z)))[1], indecs)
    O = indecs[i_O]
    indecs_without_one = filter(e -> (O != e), indecs)

    for X ∈ indecs_without_one, Y ∈ indecs_without_one
        base_ZXY = basis(Hom((Z⊗X)⊗Y, X⊗(Y⊗Z)))
        
        length(base_ZXY) == 0 && continue

        tops = [compose(
            eᵢ(X)⊗id(Y),
            associator(X,Z,Y),
            id(X) ⊗ eⱼ(Y)
        ) for eᵢ ∈ nat_trans, eⱼ ∈ nat_trans]
        
        coeffs = [express_in_basis(t, base_ZXY) for t ∈ tops]
        ab = [a*b for a in x, b in x]
        e =  [a .* c for ((a), c) ∈ zip(ab, coeffs)]

        e = reduce(.+, e)

        bottoms = [compose(
            associator(Z,X,Y),
            eᵢ(X⊗Y),
            associator(X,Y,Z)
        ) for eᵢ ∈ nat_trans]

        coeffs = [express_in_basis(b, base_ZXY) for b ∈ bottoms]

        e2 = [a .* c for (a,c) ∈ zip(x,coeffs)]

        e2 = reduce(.+, e2)

        eqs = [eqs; e .- e2]
    end

    end_Z = basis(End(Z))
    one_coeffs = [express_in_basis(eᵢ(O), end_Z) for eᵢ ∈ nat_trans]
    id_coeffs = express_in_basis(id(Z), end_Z)

    one_eqs = reduce(.+, [a .* c for (a,c) ∈ zip(x,one_coeffs)]) .- id_coeffs
    ideal(unique([eqs; one_eqs]))
end


function build_center_ideal(Z::Object, simples::Vector = simples(parent(Z)))
    #@assert is_semisimple(parent(Z)) "Not semisimple"

    Homs = [Hom(Z⊗Xi, Xi⊗Z) for Xi ∈ simples]
    n = length(simples)
    ks = [int_dim(Homs[i]) for i ∈ 1:n]

    var_count = sum([int_dim(H) for H ∈ Homs])

    K = base_ring(Z)
    R,x = polynomial_ring(K, var_count, internal_ordering = :lex)

    # For convinience: build arrays with the variables xi
    vars = []
    q = 1
    for i ∈ 1:n
        m = int_dim(Homs[i])
        vars = [vars; [x[q:q+m-1]]]
        q = q + m
    end

    eqs = []

    one_index = findfirst(e -> is_isomorphic(one(parent(Z)), e)[1], simples)

    for k ∈ 1:n, i ∈ 1:n, j ∈ 1:n
        if i == one_index || j == one_index continue end

        base = basis(Hom(Z⊗simples[k], simples[i]⊗(simples[j]⊗Z)))

        for t ∈ basis(Hom(simples[k], simples[i]⊗simples[j]))

            l1 = [zero(R) for i ∈ base]
            l2 = [zero(R) for i ∈ base]

            for ai ∈ 1:int_dim(Homs[k])
                a = basis(Homs[k])[ai]
                l1 = l1 .+ (vars[k][ai] .* K.(express_in_basis(associator(simples[i],simples[j],Z)∘(t⊗id(Z))∘a, base)))
            end
            for bi ∈ 1:int_dim(Homs[j]), ci ∈ 1:int_dim(Homs[i])
                b,c = basis(Homs[j])[bi], basis(Homs[i])[ci]
                l2 = l2 .+ ((vars[j][bi]*vars[i][ci]) .* K.(express_in_basis((id(simples[i])⊗b)∘associator(simples[i],Z,simples[j]) ∘ (c⊗id(simples[j])) ∘ inv_associator(Z,simples[i],simples[j]) ∘ (id(Z) ⊗ t), base)))
            end
            push!(eqs, l1 .-l2)
        end
    end
    ideal_eqs = []
    for p ∈ eqs
        push!(ideal_eqs, p...)
    end

    I = ideal([f for f ∈ unique(ideal_eqs) if f != 0])

    #Require e_Z(1) = id(Z)
    
    one_c = K.(express_in_basis(id(Z), basis(End(Z))))
    push!(ideal_eqs, (vars[one_index] .- one_c)...)

    I = ideal([f for f ∈ unique(ideal_eqs) if f != 0])
end

function braidings_from_ideal(Z::Object, I::Ideal, simples::Vector{<:Object}, C)
    Homs = [Hom(Z⊗Xi, Xi⊗Z) for Xi ∈ simples]
    I = rational_lift(I)
    coeffs = recover_solutions(real_solutions(I),base_ring(Z))
    ks = [int_dim(H) for H ∈ Homs]
    centrals = CenterObject[]

    for c ∈ coeffs
        k = 1
        ex = Morphism[]
        c = [k for k ∈ c]
        for i ∈ 1:length(simples)
            if ks[i] == 0 continue end
   
            e = sum(c[k:k + ks[i] - 1] .* basis(Homs[i]))
            ex = [ex ; e]
            k = k + ks[i]
        end
        centrals = [centrals; CenterObject(C, Z, (ex))]
    end
    return centrals
end

"""
    half_braidings(Z::Object)

Return all objects in the center lying over ```Z```.
"""
function half_braidings(Z::Object; simples = simples(parent(Z)), parent = center(parent(Z)))

    I = build_center_ideal(Z,simples)

    d = dim(I)

    if d < 0 return CenterObject[] end

    if d == 0 return braidings_from_ideal(Z,I,simples, parent) end

    solutions = guess_solutions(Z,I,simples,CenterObject[],gens(base_ring(I)),d, parent)

    if length(solutions) == 0
        return CenterObject[]
    end
    unique_sols = solutions[1:1]

    for s ∈ solutions[2:end]
        if sum([dim(Hom(s,u)) for u ∈ unique_sols]) == 0
            unique_sols = [unique_sols; s]
        end
    end
    return unique_sols
end

function guess_solutions(Z::Object, I::Ideal, simples::Vector{<:Object}, solutions::Vector{CenterObject}, vars, d = dim(I), C = center(parent(Z)))
    for y in vars
        J = I + ideal([y*(y^2-1)])
        d2 = dim(J)
        if d2 == 0
            return [solutions; braidings_from_ideal(Z,J,simples,C)]
        elseif d2 < 0
            return solutions
        else
            vars_new = filter(e -> e != y, vars)
            return [solutions; guess_solutions(Z,J,simples,solutions,vars_new,d2,C)]
        end
    end
end

function center_simples(C::CenterCategory, simples = simples(C.category))
    d = dim(C.category)^2

    simples_indices = []
    c_simples = CenterObject[]
    d_max = Int(QQ(ceil(fpdim(C.category))))
    d_rem = d
    k = length(simples)

    coeffs = [i for i ∈ Base.product([0:d_max for i ∈ 1:k]...)][:][2:end]

    for c ∈ sort(coeffs, by = t -> (sum(t),length(t) - length([i for i ∈ t if i != 0])))
        if sum((c .* dim.(simples)).^2) > d_rem continue end

        if simples_covered(c,simples_indices) continue end

        X = direct_sum([simples[j]^c[j] for j ∈ 1:k])[1]

        ic = is_central(X)

        if ic
            so = half_braidings(X, simples = simples, parent = C)
            c_simples = [c_simples; so]
            d_rem = d_rem - sum([dim(x)^2 for x in so])
            if d_rem == 0 return c_simples end
            push!(simples_indices, c)
        end
    end
    if d_rem > 0
        @warn "Not all halfbraidings found"
    end
    return c_simples
end

# function monoidal_completion(simples::Vector{CenterObject})
#     complete_simples = simples
#     for i ∈ 1:length(simples)
#         for j ∈ i:length(simples)
#             X,Y = simples[[i,j]]
#             complete_simples = [complete_simples; [x for (x,m) ∈ simple_subobjects(X⊗Y)]]
#             @show complete_simples
#             complete_simples = unique_simples(complete_simples)
#         end
#     end
#     if length(complete_simples) > length(simples)
#         return monoidal_completion(complete_simples)
#     end
#     return complete_simples
# end

function simples_covered(c::Tuple, v::Vector)
    for w ∈ v
        if *((w .<= c)...)
            return true
        end
    end
    false
end

function is_independent(c::Vector,v::Vector...)
    if length(v) == 0 return true end
    m = matrix(ZZ, [vi[j] for vi ∈ v, j ∈ 1:length(v[1])])

    try
        x = solve(m,matrix(ZZ,c))
    catch
        return true
    end

    return !(*((x .>=0)...))
end

function find_centrals(simples::Vector{<:Object})
    c_simples = typeof(simples[1])[]
    non_central = typeof(simples[1])[]
    for s ∈ simples
        ic, so = is_central(s)
        if ic
            c_simples = [c_simples; so]
        else
            non_central = [non_central; s]
        end
    end
    return c_simples, non_central
end

function partitions(d::Int64,k::Int64)
    parts = []
    for c ∈ Base.product([0:d for i ∈ 1:k]...)
        if sum([x for x ∈ c]) == d
            parts = [parts; [[x for x ∈ c]]]
        end
    end
    return parts
end

"""
    braiding(X::CenterObject, Y::CenterObject)

Return the braiding isomorphism ```γ_X(Y): X⊗Y → Y⊗X```.
"""
function braiding(X::CenterObject, Y::CenterObject)
    # dom = X.object⊗Y.object
    # cod = Y.object⊗X.object
    # braid = zero_morphism(dom, cod)
    # for (s,ys) ∈ zip(simples(parent(X).category), X.γ)
    #     proj = basis(Hom(Y.object,s))
    #     if length(proj) == 0 continue end
    #     incl = basis(Hom(s,Y.object))
    #     braid = braid + sum([(i⊗id(X.object))∘ys∘(id(X.object)⊗p) for i ∈ incl, p ∈ proj][:])
    # end
    braid = half_braiding(X,object(Y))
    return morphism(X⊗Y,Y⊗X,braid)
end

@doc raw""" 

    half_braiding(X::CenterObject, Y::Object)

Return the half braiding isomorphism ```γ_X(Y): X⊗Y → Y⊗X```.
"""
function half_braiding(X::CenterObject, Y::Object)
    simpls = simples(parent(Y))

    if is_simple(Y) 
        if !(Y ∈ simpls)
            k = findfirst(e -> is_isomorphic(e, Y)[1], simpls)
            iso = is_isomorphic(Y,simpls[k])[2]
            return (inv(iso)⊗id(X.object)) ∘ X.γ[k] ∘ (id(X.object)⊗iso)
        else
            k = findfirst(e -> e == Y, simpls)
            return X.γ[k]
        end
    end
    dom = X.object⊗Y
    cod = Y⊗X.object
    braid = zero_morphism(dom, cod)
   
  
    _,iso, incl, proj = direct_sum_decomposition(Y)

    for (p,i) ∈ zip(proj, incl)
        k = findfirst(e -> is_isomorphic(e, domain(i))[1], simpls)
        incliso = is_isomorphic(domain(i), domain(i))[2]
        #projiso = is_isomorphic(codomain(p), )[2]

        i = i ∘ incliso
        p = inv(incliso) ∘ p 

        braid = braid + (i⊗id(X.object))∘X.γ[k]∘(id(X.object)⊗p)
    end

    return braid
end


#-------------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------------

"""
    dim(X::CenterObject)

Return the categorical dimension of ```X```.
"""
dim(X::CenterObject) = dim(X.object)

"""
    fpdim(X::CenterObject)

Return the Frobenius-Perron dimension of ```X```.
"""
#fpdim(X::CenterObject) = fpdim(object(X))

"""
    simples(C::CenterCategory)

Return a vector containing the simple objects of ```C```. 
"""
function simples(C::CenterCategory; sort = false, show_progress = false)
    if isdefined(C, :simples) 
        return C.simples 
    end
    if is_modular(category(C))
        C.simples = center_simples_by_braiding(category(C), C)
        return C.simples
    end
    simples_by_induction!(C, show_progress)
    if sort 
        sort_simples_by_dimension!(C)
    end
    return C.simples
end


function decompose(X::CenterObject)
    C = parent(X)
    K = base_ring(X)

    if K == QQBarField() || typeof(K) == CalciumField
        return decompose_over_qqbar(X)
    end

    if isdefined(C, :simples) && is_semisimple(C)
        return decompose_by_simples(X,simples(C))
    else
        try
            return decompose_by_endomorphism_ring(X)
        catch
            @assert is_semisimple(C)
            error("cannot decompose")
            indecs = indecomposable_subobjects(X)
            return [(x, div(int_dim(Hom(x,X)), int_dim(End(x)))) for x ∈ indecs]
        end
    end
end

function decompose(X::CenterObject, S::Vector{CenterObject})
    decompose_by_simples(X,S)
end

# function indecomposable_subobjects(X::CenterObject)

#     !is_split_semisimple(category(parent(X))) && return _indecomposable_subobjects(X)


#     B = basis(End(object(X)))

#     if length(B) == 1
#         return [X]
#     end

#     S = simples(parent(object(X)))

#     if length(B) ≤ length(S)^2
#         return _indecomposable_subobjects(X)
#     end

#     eig_spaces = []


#     while length(eig_spaces) ≤ 1 && length(B) > 0
#         f = popat!(B, rand(eachindex(B)))
#         proj_f = central_projection(X,X,f,S)
#         eig_spaces = collect(values(eigenvalues(proj_f)))
#     end

#     if length(eig_spaces) ≤ 1 
#         is_simple(X) && return [X]
#         return [x for (x,k) ∈ decompose(X)]
#     end

#     return unique_simples(vcat([indecomposable_subobjects(Y) for Y ∈ eig_spaces]...))
# end

# function indecomposable_subobjects_of_induction(X::Object, IX::CenterObject = induction(X))
#     @assert object(IX) == X
#     B = basis(Hom(X, object(IX)))

#     while length(B) > 0
#         f = popat!(B, rand(eachindex(B)))
#         f = induction_adjunction(f,IX,IX)
        
#         eig = collect(values(eigenspaces(f)))

#         length(B) ≥ 2 && break
#     end
# end


"""
    associator(X::CenterObject, Y::CenterObject, Z::CenterObject)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
function associator(X::CenterObject, Y::CenterObject, Z::CenterObject)
    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)
    return morphism(dom,cod, associator(X.object, Y.object, Z.object))
end

matrices(f::CenterMorphism) = matrices(f.m)
matrix(f::CenterMorphism) = matrix(f.m)

"""
    compose(f::CenterMorphism, g::CenterMorphism)

Return the composition ```g∘f```.
"""
compose(f::CenterMorphism, g::CenterMorphism) = morphism(domain(f), codomain(g), g.m ∘ f.m) 

"""
    dual(X::CenterObject)

Return the (left) dual object of ```X```.
"""
function dual(X::CenterObject)
    a = associator
    inv_a = inv_associator
    e = ev(X.object)
    c = coev(X.object)
    γ = Morphism[]
    dX = dual(X.object)
    for (Xi,yXi) ∈ zip(simples(parent(X).category), X.γ)
        f = (e ⊗ id(Xi ⊗ dX)) ∘ 
            inv_a(dX, X.object, Xi ⊗ dX) ∘ 
            (id(dX) ⊗ a(X.object, Xi, dX)) ∘ 
            (id(dX) ⊗ (inv(yXi) ⊗ id(dX))) ∘ 
            (id(dX) ⊗ inv_a(Xi, X.object, dX)) ∘ 
            a(dX, Xi, X.object ⊗ dX) ∘ 
            (id(dX ⊗ Xi) ⊗ c)
        γ = [γ; f]
    end
    return CenterObject(parent(X),dX,γ)
end

"""
    ev(X::CenterObject)

Return the evaluation morphism ``` X⊗X → 1```.
"""
function ev(X::CenterObject)
    morphism(dual(X)⊗X,one(parent(X)),ev(X.object))
end

"""
    coev(X::CenterObject)

Return the coevaluation morphism ```1 → X⊗X∗```.
"""
function coev(X::CenterObject)
    morphism(one(parent(X)),X⊗dual(X),coev(X.object))
end

"""
    id(X::CenterObject)

Return the identity on ```X```.
"""
id(X::CenterObject) = morphism(X,X,id(X.object))

"""
    tr(f:::CenterMorphism)

Return the categorical trace of ```f```.
"""
function tr(f::CenterMorphism)
    C = parent(domain(f))
    return CenterMorphism(one(C),one(C),tr(f.m))
end

"""
    inv(f::CenterMorphism)

Return the inverse of ```f```if possible.
"""
function inv(f::CenterMorphism)
    return morphism(codomain(f),domain(f), inv(f.m))
end


"""
    is_isomorphic(X::CenterObject, Y::CenterObject)

Check if ```X≃Y```. Return ```(true, m)``` where ```m```is an isomorphism if true,
else return ```(false,nothing)```.
"""
function is_isomorphic(X::CenterObject, Y::CenterObject)
    # TODO: Fix This. How to compute a central isomorphism?

    if ! is_isomorphic(object(X),object(Y))[1]
        return false, nothing
    end
    # if is_simple(X) && is_simple(Y)
    #     H = hom_by_linear_equations(X,Y)
    #     if int_dim(H) > 0
    #         return true, basis(H)[1]
    #     else
    #         return false, nothing
    #     end
    # end

    EX = endomorphism_ring(X) 
    EY = endomorphism_ring(Y) 

    if dim(EX) != dim(EY) 
        return false, nothing 
    end

    S = simples(parent(X))

    if [dim(Hom(X,s)) for s ∈ S] == [dim(Hom(Y,s)) for s ∈ S]
        _, iso = is_isomorphic(X.object, Y.object)
        return true, central_projection(X,Y,[iso])[1]
    else
        return false, nothing
    end
end

function is_isomorphic_simples(X::CenterObject, Y::CenterObject)
    !is_isomorphic(object(X), object(Y))[1] && return (false, nothing)

    H = hom_by_adjunction(X,Y)
    int_dim(H) > 0 ? (true , H[1]) : (false, nothing)
end

function +(f::CenterMorphism, g::CenterMorphism)
    return morphism(domain(f), codomain(f), g.m +f.m)
end

function *(x, f::CenterMorphism)
    return morphism(domain(f),codomain(f),x*f.m)
end
#-------------------------------------------------------------------------------
#   Functionality: Image
#-------------------------------------------------------------------------------

"""
    kernel(f::CenterMoprhism)

Return a tuple ```(K,k)``` where ```K```is the kernel object and ```k```is the inclusion.
"""
function kernel(f::CenterMorphism)
    ker, incl = kernel(f.m)
    #f_inv = left_inverse(incl)

    if ker == zero(parent(f.m))
        return zero(parent(f)), zero_morphism(zero(parent(f)), domain(f))
    end

    inv_incl = left_inverse(incl)

    braiding = [id(s)⊗inv_incl∘γ∘(incl⊗id(s)) for (s,γ) ∈ zip(simples(parent(domain(f.m))), domain(f).γ)]

    Z = CenterObject(parent(domain(f)), ker, braiding)
    return Z, morphism(Z,domain(f), incl)
end

"""
    cokernel(f::CenterMorphism)

Return a tuple ```(C,c)``` where ```C```is the cokernel object and ```c```is the projection.
"""
function cokernel(f::CenterMorphism)
    coker, proj = cokernel(f.m)
    #f_inv = right_inverse(proj)

    if coker == zero(parent(f.m))
        return zero(parent(f)), zero_morphism(codomain(f), zero(parent(f)))
    end

    inv_proj = right_inverse(proj)
    braiding = [(id(s)⊗proj)∘γ∘(inv_proj⊗id(s)) for (s,γ) ∈ zip(simples(parent(domain(f.m))), codomain(f).γ)]

    Z = CenterObject(parent(domain(f)), coker, braiding)
    return Z, morphism(codomain(f),Z, proj)
end

function image(f::CenterMorphism)
    I, incl = image(f.m)

    if I == zero(parent(f.m))
        return zero(parent(f)), zero_morphism(zero(parent(f)), domain(f))
    end

    braiding = [id(s)⊗left_inverse(incl)∘γ∘(incl⊗id(s)) for (s,γ) ∈ zip(simples(parent(I)), codomain(f).γ)]

    Z = CenterObject(parent(domain(f)), I, braiding)
    return Z, morphism(Z,domain(f), incl)
end


# function left_inverse(f::CenterMorphism)
#     X = domain(f)
#     Y = codomain(f)
#     l_inv = central_projection(Y,X,left_inverse(morphism(f)))
#     return morphism(Y,X,l_inv)
# end

function quotient(Y::CenterObject, X::Object)
    # TODO: Compute quotient
    @assert parent(X) == parent(Y).Category
end

#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct CenterHomSpace <: AbstractHomSpace
    X::CenterObject
    Y::CenterObject
    basis::Vector{CenterMorphism}
    parent::VectorSpaces
end


function Hom(X::CenterObject, Y::CenterObject) 

    alg_closed = (base_ring(X) == QQBarField() || typeof(base_ring(X)) == CalciumField)

    if alg_closed && has_attribute(X, :is_simple) && 
        get_attribute(X, :is_simple) &&
        has_attribute(Y, :is_simple) &&
        get_attribute(Y, :is_simple) &&
        return hom_by_linear_equations(X,Y)
    end

    if is_zero(dim(category(parent(X))))
        return hom_by_linear_equations(X,Y)
    else 
        return hom_by_adjunction(X,Y)
    end
end

@doc raw""" 

    central_projection(X::CenterObject, Y::CenterObject, f::Morphism)

Compute the image under the projection ```Hom(F(X),F(Y)) → Hom(X,Y)```.
"""
function central_projection(dom::CenterObject, cod::CenterObject, f::Vector{<:Morphism}, simpls = simples(category(parent(dom))))

    if length(f) == 0 
        return CenterMorphism[]
    end

    X = domain(f[1])
    Y = codomain(f[1])
    C = parent(X)
    D = dim(C)
    proj = [zero_morphism(X, Y) for _ ∈ f]
    a = associator
    inv_a = inv_associator

    for (Xi, yX) ∈ zip(simpls, dom.γ)
        dXi = dual(Xi)

        yY = half_braiding(cod, dXi)
        
        ϕ_before = (yX⊗id(dXi))∘inv_a(X,Xi,dXi)∘(id(X)⊗coev(Xi))

        ϕ_after = (ev(dXi)⊗id(Y))∘inv_a(dual(dXi),dXi,Y)∘(spherical(Xi)⊗yY)∘a(Xi,Y,dXi)

        proj = [p + dim(Xi)*(ϕ_after ∘ ((id(Xi)⊗g)⊗id(dXi)) ∘ ϕ_before) for (p,g) ∈ zip(proj,f)]
    end
    return [morphism(dom, cod, inv(D*base_ring(dom)(1))*p) for p ∈ proj]
end

"""
    zero_morphism(X::CenterObject, Y::CenterObject)

Return the zero morphism ```0:X → Y```.
"""
zero_morphism(X::CenterObject, Y::CenterObject) = morphism(X,Y,zero_morphism(X.object,Y.object))

#-------------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------------

function show(io::IO, X::CenterObject)
    print(io, "Central object: $(X.object)")
end

function show(io::IO, C::CenterCategory)
    print(io, "Drinfeld center of $(C.category)")
end

function show(io::IO, f::CenterMorphism)
    print(io, "Morphism in $(parent(domain(f)))")
end


#=------------------------------------------------
    Center by Induction
------------------------------------------------=#
function add_induction!(C::CenterCategory, X::Object, IX::CenterObject)
    if isdefined(C, :inductions)
        if !(X ∈ keys(C.inductions))
            push!(C.inductions, X => IX)
        end
    else
        C.inductions = Dict(X => IX)
    end
end

function simples_by_induction!(C::CenterCategory, log = true)
    S = CenterObject[]

    if rank(category(C)) == 1
        C.simples = S
        return S 
    end

    d = dim(C.category)^2
    C.induction_gens = induction_generators(C)
    simpls = simples(C.category)
    K = base_ring(C)

    # FI_simples = [(s, )]

    # ind_res = [induction_restriction(s) for s ∈ simpls]

    # # Group the simples by isomorphic inductions
    # is_iso = [s == t ? true : is_isomorphic(s,t)[1] for s ∈ ind_res, t ∈ ind_res]
    # groups = connected_components(graph_from_adjacency_matrix(Undirected, is_iso))
    
    # for gr ∈ groups 
    #     Is = induction(simpls[gr[1]], simpls, parent_category = C)
    #     push!(FI_simples, (simpls[gr[1]], ind_res[gr[1]], Is))
    #     push!(C.induction_gens, simpls[gr[1]])
    # end
    
    log && println("Simples:")

    #center_dim = 0
    used_gens = []
    remove_gens = []
    for s ∈ induction_generators(C)

        Is = induction_restriction(s, simpls)
        S_in_Z = []

        if dim(category(C)) != 0
            #Test which simples in S are included
            multiplicities = K == QQBarField() || typeof(K) == CalciumField ? [int_dim(Hom(object(t),s)) for t ∈ S] : [div(int_dim(Hom(object(t),s)), int_dim(End(t))) for t ∈ S]

            S_in_Z = [(t, k) for (t,k) ∈ zip(S, multiplicities) if k != 0]

            if !isempty(S_in_Z) && sum([fpdim(t)*k for (t,k) ∈ S_in_Z]) == fpdim(Is)
                push!(remove_gens, s)
                continue
            end
        end

        Z = induction(s, simpls, parent_category = C)

        if ! isempty(S_in_Z)
            #factor out all already known simples

            incl_basis = vcat([basis(induction_right_adjunction(Hom(object(t),s), t, Z)) for (t,_) ∈ S_in_Z]...)

            incl = horizontal_direct_sum(incl_basis)
            
            Q, proj = cokernel(incl)

            H = induction_right_adjunction(Hom(object(Q), s), Q, Z)
            
            H = HomSpace(Q,Q, [proj ∘ f for f ∈ basis(H)])

            Z = Q
            # for (t,k) ∈ S_in_Z
            #     incl = horizontal_direct_sum(basis(Hom(t,Z))...)
            #     @show Z = cokernel(incl)[1]
            # end
            # H = End(Z)
        else
            H = end_of_induction(s, Z)
        end

        #H = end_of_induction(s, Z)

        #contained_simples = filter(x -> int_dim(Hom(object(x),s)) != 0, S)
        # if length(contained_simples) > 0
        #     if is_isomorphic(Is, direct_sum(object.(contained_simples))[1])[1]
        #         continue
        #     end
        # end

        #Z = induction(s, simpls, parent_category = C)

        # for x ∈ contained_simples
        #     f = horizontal_direct_sum(basis(Hom(x,Z)))
        #     Z = cokernel(f)[1]
        # end

        # Compute subobjects by computing central primitive central_primitive_idempotents of End(I(X))
        
        # idems = central_primitive_idempotents(H)
        new_simples = simple_subobjects(Z,H)

        # Every simple such that Hom(s, Zᵢ) ≠ 0 for an already dealt with s is not new
        filter!(Zi -> sum(Int[int_dim(Hom(s,object(Zi))) for s ∈ used_gens]) == 0, new_simples)

        # if length(new_simples) == 0
        #     continue
        # end
        
        log && println.(["    " * "$s" for s ∈ new_simples])

        S = [S; new_simples]
        used_gens = [used_gens; s]
        #center_dim += sum(dim.(new_simples).^2)
        # if d == center_dim
        #     break
        # end
    end
    filter!(e -> e ∉ remove_gens, C.induction_gens)

    C.simples = S
end

function sort_simples_by_dimension!(C::CenterCategory)  
    
    one_ind = findfirst(==(one(C)), C.simples)

    C.simples[[1,one_ind]] = C.simples[[one_ind, 1]]
    fp_dims = [fpdim(s) for s ∈ simples(C)]

    σ = sortperm(fp_dims, by = abs)

    has_attribute(C, :multiplication_table) && delete!(C.__attrs, :multiplication_table)
    has_attribute(C, :smatrix) && delete!(C.__attrs, :smatrix)


    C.simples = C.simples[σ]
end

function split(X::CenterObject, E = End(X), 
    C = _extension_of_scalars(parent(X), QQBarField(), 
    extension_of_scalars(category(parent(X)), QQBarField())), 
    e = complex_embeddings(base_ring(X))[1])
    #Assume X simple
    # int_dim(E) ≤ 1 && return [X]

    # g = gens(E)
    # f = g[findmax(degree ∘ minpoly, g)[2]]

    # L = splitting_field(minpoly(f))
    # collect(values(eigenvalues(f ⊗ L)))
    if base_ring(X) == QQ 
        to_qqbar = QQBarField()
    else
        to_qqbar = x -> guess(QQBarField(), e(x,512), maximum([1,degree(x)]))
    end
    # C = _extension_of_scalars(parent(X), QQBarField(), extension_of_scalars(category(parent(X)), QQBarField(), embedding = to_qqbar))

    ext_X = extension_of_scalars(X, QQBarField(), C, embedding = to_qqbar)
    ext_E = HomSpace(ext_X,ext_X, [morphism(ext_X,ext_X, extension_of_scalars(morphism(f), QQBarField(), category(C), embedding = to_qqbar)) for f ∈ basis(E)])

    simple_subobjects(ext_X, ext_E)
end

function split_cyclotomic(C::CenterCategory; absolute = false)
    Ends = End.(simples(C))
    K = base_ring(C)

    if all([int_dim(e) == 1 for e ∈ Ends])
        return C
    end

    n = []

    # find smallest cyclotomic extension such that Z(C) splits
    for e ∈ Ends
        if int_dim(e) > 1
            i = findfirst(f -> degree(minpoly(f)) == int_dim(e), basis(e))
           
            f = if i === nothing 
                r = endomorphism_ring(domain(e),e)
                 r_gens = gens(r)
                 _,i = findmax(degree.(minpoly.(r_gens)))
                 minpoly(r_gens[i])
            else
                minpoly(e[i])
            end

            m = max(int_dim(e),4)
            while true 
                L = if K == QQ 
                        cyclotomic_field(m)[1]
                    else
                        cyclotomic_extension(K,m).Kr
                    end
               
                rs = roots(L, f)
                if length(rs) > 0 
                    push!(n,m) 
                    break
                end
                m += 1
            end
            
        end
    end
    m = lcm(n...)
    L = if K == QQ 
            cyclotomic_field(m)[1]
        else
            L2 = cyclotomic_extension(K,m)
            absolute ? L2.Ka : L2.Kr
        end
    

    if absolute
        extension_of_scalars(C,L)
    else
        extension_of_scalars(C, L, embedding = x -> L(x))
    end
end


function split_algebraic(C::CenterCategory, e = complex_embeddings(base_ring(C))[1])
    ends = End.(simples(C))

    non_split = filter(e -> int_dim(e) > 1, ends)

    length(non_split) == 0 && return C 

    # minpolys = [minpoly.(basis(e)) for e ∈ non_split]

    # K = base_ring(C) 
    # _,x = K[:x]
    # minpolys = [B[findfirst(f -> degree(f) == length(B), B)](x) for B ∈ minpolys]

    # L = simplify(absolute_simple_field(splitting_field(minpolys))[1])[1]

    # incl = K == QQ ? L : is_subfield(K,L)[2]
    if base_ring(C) == QQ 
        to_qqbar = QQBarField()
    else
        to_qqbar = x -> guess(QQBarField(), e(x,512), maximum([degree(x),1]))
    end
    CL = _extension_of_scalars(C, QQBarField(), extension_of_scalars(category(C), QQBarField(), embedding = to_qqbar))

    simpls = vcat([split(s, E, CL, e) for (s,E) ∈ zip(simples(C), ends)]...)

    CL.simples = simpls

    return CL
end


#=----------------------------------------------------------
    Hom Spaces 2.0 
----------------------------------------------------------=#


function hom_by_adjunction(X::CenterObject, Y::CenterObject)
    Z = parent(X)
    C = category(Z)
    S = induction_generators(Z)

    X_Homs = [Hom(object(X),s) for s ∈ S]
    Y_Homs = [Hom(s,object(Y)) for s ∈ S]

    candidates = [int_dim(H)*int_dim(H2) > 0 for (H,H2) ∈ zip(X_Homs,Y_Homs)]

    !any(candidates) && return HomSpace(X,Y, CenterMorphism[]) 

    
    # X_Homs = X_Homs[candidates]
    # Y_Homs = Y_Homs[candidates]
    

    M = zero_matrix(base_ring(C),0,*(size(matrix(zero_morphism(X,Y)))...))

    mors = []

    @threads for i ∈ findall(==(true), candidates)
        s, X_s, s_Y = S[i], X_Homs[i], Y_Homs[i]
        Is = induction(s, parent_category = Z)
        B = induction_right_adjunction(X_s, X, Is)
        B2 = induction_adjunction(s_Y, Y, Is)
        
        # Take all combinations
        B3 = [h ∘ b for b ∈ B, h in B2][:]

        mors = [mors; B3]
        # Build basis
    end
    mats = matrix.(mors)
    M = transpose(matrix(base_ring(C), hcat(hcat([collect(m)[:] for m in mats]...))))

    Mrref = hnf(M)
    base = CenterMorphism[]
    mats_morphisms = morphism.(mats)

    for k ∈ 1:rank(Mrref)
        coeffs = express_in_basis(morphism(transpose(matrix(base_ring(C), size(mats[1])..., Mrref[k,:]))), mats_morphisms)
        f = sum([m*bi for (m,bi) ∈ zip(coeffs, mors)])
        push!(base, f)
    end

    return HomSpace(X,Y, base)
end


function hom_by_linear_equations(X::CenterObject, Y::CenterObject)
    #@assert parent(X) == parent(Y)

    H = Hom(object(X), object(Y))
    B = basis(H)
    F = base_ring(X)
    n = length(basis(H))

    if n == 0 
        return HomSpace(X,Y, CenterMorphism[])
    end 

    Fx,poly_basis = polynomial_ring(F,n)
    
    eqs = []

    S = simples(parent(object(X)))

    for (s,γₛ,λₛ) ∈ zip(S,half_braiding(X), half_braiding(Y))
        Hs = Hom(object(X)⊗s, s⊗object(Y))
        base = basis(Hs)
        if length(base) == 0
            continue
        end
        eq_i = [zero(Fx) for _ ∈ 1:length(base)]
        for (f,a) ∈ zip(B,poly_basis)
            coeffs = express_in_basis((id(s)⊗f)∘γₛ - λₛ∘(f⊗id(s)), base)
            eq_i = eq_i .+ (a .* coeffs)
        end
        
        eqs = [eqs; eq_i]

    end

    M = zero(matrix_space(F,length(eqs),n))

    for (i,e) ∈ zip(1:length(eqs),eqs)
        M[i,:] = [coeff(e, a) for a ∈ poly_basis]
    end

    N = nullspace(M)[2]

    _,cols = size(N)

    basis_coeffs = [collect(N[:,i]) for i ∈ 1:cols]

    center_basis = [CenterMorphism(X,Y,sum(b .* B)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,center_basis)
end

function hom_by_projection(X::CenterObject, Y::CenterObject)
    b = basis(Hom(X.object, Y.object))

    projs = central_projection(X,Y,b)

    proj_exprs = [express_in_basis(morphism(p),b) for p ∈ projs]

    M = zero(matrix_space(base_ring(X), length(b),length(b)))
    for i ∈ 1:length(proj_exprs)
        M[i,:] = proj_exprs[i]
    end
    r, M = rref(M)
    H_basis = CenterMorphism[]
    for i ∈ 1:r
        f = morphism(X,Y,sum([m*bi for (m,bi) ∈ zip(M[i,:], b)]))
        H_basis = [H_basis; f]
    end
    return HomSpace(X,Y,H_basis)
end


#=----------------------------------------------------------
    Modular Stuff 
----------------------------------------------------------=#    

function smatrix(C::CenterCategory)
    get_attribute!(C, :smatrix) do
        simpls = simples(C)
        n = length(simpls)
        K = base_ring(C)
        S = [zero_morphism(category(C)) for _ ∈ 1:n, _ ∈ 1:n]
        @threads for i ∈ 1:n
            for j ∈ i:n
                S[i,j] = S[j,i] = tr(half_braiding(simpls[i], object(simpls[j])) ∘ half_braiding(simpls[j], object(simpls[i])))
            end
        end

        try
            return matrix(K, n, n, [K(s) for s ∈ S])
        catch
            return S
        end
    end
end

#=----------------------------------------------------------
    extension_of_scalars 
----------------------------------------------------------=#    

function extension_of_scalars(C::CenterCategory, L::Field; embedding = is_subfield(base_ring(C), L)[2])

    CL = _extension_of_scalars(C,L, extension_of_scalars(category(C), L, embedding = embedding))

    S = simples(C)
    Ends = End.(S)

    new_simples = CenterObject[]

    for i ∈ eachindex(S)
        if int_dim(Ends[i]) == 1
            # All simples with End(s) == k can be just written over
            push!(new_simples, extension_of_scalars(S[i], L, CL, embedding = embedding))
        else
            s_L = extension_of_scalars(S[i], L, CL, embedding = embedding)
            h = HomSpace(s_L, s_L, [morphism(s_L,s_L, extension_of_scalars(morphism(f), L, category(CL), embedding = embedding)) for f ∈ Ends[i]])
            # Split all other simples by decomposing the extension
            append!(new_simples, simple_subobjects(extension_of_scalars(S[i], L, CL, embedding = embedding), h))
        end
    end

    CL.simples = new_simples
    # if isdefined(C, :inductions)
    #     CL.inductions = Dict(extension_of_scalars(x, L, category(CL), embedding = embedding) =>
    #                     extension_of_scalars(Ix, L, CL, embedding = embedding) for (x,Ix) ∈ C.inductions)
    # end

    if isdefined(C, :induction_gens)
        CL.induction_gens = [extension_of_scalars(is, L, category(CL),  embedding = embedding) for is ∈ C.induction_gens]
    end

    return CL
end

function _extension_of_scalars(C::CenterCategory, L::Field, cL = category(C)⊗L)
    CenterCategory(L,cL)
end

function extension_of_scalars(X::CenterObject, L::Field, CL::CenterCategory = _extension_of_scalars(parent(X),L);  embedding = embedding(base_ring(X), L))
    CenterObject(CL, extension_of_scalars(object(X), L, category(CL),  embedding = embedding), [extension_of_scalars(f, L, category(CL),  embedding = embedding) for f ∈ half_braiding(X)])
end

function extension_of_scalars(f::CenterMorphism, L::Field, CL::CenterCategory = _extension_of_scalars(parent(f),L),  embedding = embedding(base_ring(f), L))
    dom = extension_of_scalars(domain(f), L, CL,  embedding = embedding)
    cod = extension_of_scalars(codomain(f), L, CL,  embedding = embedding)
    m = extension_of_scalars(morphism(f), L, category(CL),  embedding = embedding)
    CenterMorphism(dom, cod, m)
end



function karoubian_envelope(C::CenterCategory)
    KC = CenterCategory(base_ring(C), category(C))
    simpls = unique_simples(vcat([simple_subobjects(s) for s ∈ simples(C)]...))
    KC.simples = [CenterObject(KC, object(s), half_braiding(s)) for s ∈ simpls]
    return KC
end


#=----------------------------------------------------------
    Center for non-degenerate braided fusion categories
    by C ⊠ C^rev ≃ 𝒵(C) 
----------------------------------------------------------=#

function center_simples_by_braiding(C::Category, Z::CenterCategory = center(C))
    S = simples(C)

    S_braided = [CenterObject(Z, s, [braiding(s,t) for t ∈ S]) for s ∈ S]
    S_rev_braided = [CenterObject(Z, s, [inv(braiding(t,s)) for t ∈ S]) for s ∈ S]

    [t⊗s for s ∈ S_braided, t ∈ S_rev_braided][:]
end

#=----------------------------------------------------------
    Drinfeld Morphism 
----------------------------------------------------------=#

function drinfeld_morphism(X::CenterObject) 
    morphism(X,dual(dual(X)), _drinfeld_morphism(X))
end

function _drinfeld_morphism(X::CenterObject)
    x = object(X)
    u = (ev(x)⊗id(dual(dual(x)))) ∘ 
        (half_braiding(X,dual(x))⊗id(dual(dual(x)))) ∘ 
        inv_associator(x, dual(x), dual(dual(x))) ∘ 
        (id(x)⊗coev(dual(x)))
end

function twist(X::CenterObject)
    u = _drinfeld_morphism(X)
    
    B,k = is_scalar_multiple(matrix(spherical(object(X))), matrix(u))

    !B && error("Something went wrong")

    return k
end


#=----------------------------------------------------------
    Skeletalization 
----------------------------------------------------------=#

function multiplication_table(C::CenterCategory)
    get_attribute!(C, :multiplication_table) do

        if characteristic(base_ring(C)) > 0 || !is_split_semisimple(C) || !is_spherical(C)
            return multiplication_table(simples(C))
        end

        S = simples(C) 
        dims = dim.(S)
        d = dim(C)



        S_matrix = smatrix(C)

        n = length(S) 

        multiplicities = Array{Int,3}(undef,n,n,n)

        #duals = [findfirst(e -> is_isomorphic(e,dual(s))[1], S) for s ∈ S]

        duals = [findfirst(!=(0), e) for e ∈ eachrow(S_matrix^2)]

        for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n 
            verlinde_formula = sum([*(S_matrix[l,[i,j,duals[k]]]...)//dims[l] for l ∈ 1:n])

            if typeof(base_ring(C)) <: Union{PadicField, QadicField}
                multiplicities[i,j,k] = Int(lift(coordinates(verlinde_formula//d)[1]))
            else
                multiplicities[i,j,k] = Int(QQ(verlinde_formula//d))
            end
        end
        return multiplicities
    end
end

# function multiplicity_spaces(C::CenterCategory, S = simples(C))
#     get_attribute!(C, :multiplicity_spaces) do
#         n = length(S)
#         d = dim(C)

#         if ! is_split_semisimple(C) 
#             Sij = S[i]⊗S[j]
#             return [Hom(Sij, S[k]) for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n]
#         end

#         multiplicities = multiplication_table(C)

#         spaces = [HomSpace(S[i],S[j], [zero_morphism(S[i],S[j]) for _ ∈ 1:multiplicities[i,j,k]]) for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n]

#         return spaces
#     end
# end


function six_j_symbols(C::CenterCategory, S = simples(C), mult = nothing)
    @assert is_semisimple(C)

    six_j_symbols_of_construction(C, S, multiplication_table(C))
end


function simples_names(C::CenterCategory)
    S = simples(C) 
    n = length(S)
    names = ["" for i ∈ S]
    
    not_taken = [i for i ∈ 1:n]
    while ! isempty(not_taken)
        i = popfirst!(not_taken) 
        s = object(S[i])
        
        iso = findall(e -> is_isomorphic(object(S[i]), object(e))[1], S)

        if length(iso) == 1 
            names[i] = "($(s), γ)"
            continue
        end

        l = 1
        for j ∈ iso
            names[j] = "($s, γ$l)"
            l += 1
        end
    end

    return names
end