mutable struct CenterCategory <: Category
    base_ring::Field
    category::Category
    simples::Vector{O} where O <: CategoryObject

    function CenterCategory(F::Field, C::Category)
        Z = new()
        Z.base_ring = F
        Z.category = C
        return Z
    end

    function CenterCategory()
        new()
    end
end

struct CenterCategoryObject <: CategoryObject
    parent::CenterCategory
    object::CategoryObject
    γ::Vector{M} where M <: CategoryMorphism
end

struct CenterCategoryMorphism <: CategoryMorphism
    domain::CenterCategoryObject
    codomain::CenterCategoryObject
    m::CategoryMorphism
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

function isequal_without_parent(X::CenterCategoryObject, Y::CenterCategoryObject)
    return object(X) == object(Y) && half_braiding(X) == half_braiding(Y)
end
#-------------------------------------------------------------------------------
#   Center Constructor
#-------------------------------------------------------------------------------
"""
    Center(C::Category)

Return the Drinfeld center of ```C```.
"""
function Center(C::Category; equivalence = false)
    @assert is_semisimple(C) "Semisimplicity required"
    return CenterCategory(base_ring(C),C)
end

function Morphism(dom::CenterCategoryObject, cod::CenterCategoryObject, m::CategoryMorphism)
    return CenterCategoryMorphism(dom,cod,m)
end

"""
    half_braiding(Z::CenterCategoryObject)

Return  a vector with half braiding morphisms ```Z⊗S → S⊗Z``` for all simple
objects ```S```.
"""
half_braiding(Z::CenterCategoryObject) = Z.γ


"""
    object(X::CenterCategoryObject)

Return the onderlying object in ```𝒞```.
"""
object(X::CenterCategoryObject) = X.object

morphism(f::CenterCategoryMorphism) = f.m

is_fusion(C::CenterCategory) = true

"""
    add_simple!(C::CenterCategory, S::CenterCategoryObject)

Add the simple object ```S``` to the vector of simple objects.
"""
function add_simple!(C::CenterCategory, S::CenterCategoryObject)
    @assert dim(End(S)) == 1 "Not simple"
    if isdefined(C, :simples)
        C.simples = unique_simples([simples(C); S])
    else
        C.simples = unique_simples([S])
    end
end

function add_simple!(C::CenterCategory, S::Array{CenterCategoryObject})
    @assert prod(dim(End(s)) for s ∈ S) == 1 "Not simple"
    if isdefined(C, :simples)
        C.simples = unique_simples([simples(C); S])
    else
        C.simples = unique_simples(S)
    end
end
"""
    spherical(X::CenterCategoryObject)

Return the spherical structure ```X → X∗∗``` of ```X```.
"""
spherical(X::CenterCategoryObject) = Morphism(X,dual(dual(X)), spherical(X.object))

(F::Field)(f::CenterCategoryMorphism) = F(f.m)

#=-------------------------------------------------
    MISC 
-------------------------------------------------=#

==(f::CenterCategoryMorphism, g::CenterCategoryMorphism) = f.m == g.m

#-------------------------------------------------------------------------------
#   Direct Sum & Tensor Product
#-------------------------------------------------------------------------------

"""
    direct_sum(X::CenterCategoryObject, Y::CenterCategoryObject)

Return the direct sum object of ```X``` and ```Y```.
"""
function direct_sum(X::CenterCategoryObject, Y::CenterCategoryObject)
    S = simples(parent(X.object))
    Z,(ix,iy),(px,py) = direct_sum(X.object, Y.object)

    γZ = [(id(S[i])⊗ix)∘(X.γ[i])∘(px⊗id(S[i])) + (id(S[i])⊗iy)∘(Y.γ[i])∘(py⊗id(S[i])) for i ∈ 1:length(S)]

    CZ = CenterCategoryObject(parent(X), Z, γZ)
    ix,iy = CenterCategoryMorphism(X,CZ,ix), CenterCategoryMorphism(Y,CZ, iy)
    px,py = CenterCategoryMorphism(CZ,X,px), CenterCategoryMorphism(CZ,Y,py)
    return CZ,[ix,iy],[px,py]
end



"""
    direct_sum(f::CenterCategoryMorphism, g::CenterCategoryMorphism)

Return the direct sum of ```f``` and ```g```.
"""
function direct_sum(f::CenterCategoryMorphism, g::CenterCategoryMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    m = f.m ⊕ g.m
    return Morphism(dom,cod, m)
end

"""
    tensor_product(X::CenterCategoryObject, Y::CenterCategoryObject)

Return the tensor product of ```X``` and ```Y```.
"""
function tensor_product(X::CenterCategoryObject, Y::CenterCategoryObject)
    Z = X.object ⊗ Y.object
    γ = CategoryMorphism[]
    a = associator
    inv_a = inv_associator
    s = simples(parent(X.object))
    x,y = X.object, Y.object
    for (S, yX, yY) ∈ zip(s, X.γ, Y.γ)
        push!(γ, a(S,x,y)∘(yX⊗id(y))∘inv_a(x,S,y)∘(id(x)⊗yY)∘a(x,y,S))
    end
    return CenterCategoryObject(parent(X), Z, γ)
end

"""
    tensor_product(f::CenterCategoryMorphism,g::CenterCategoryMorphism)

Return the tensor product of ```f``` and ```g```.
"""
function tensor_product(f::CenterCategoryMorphism,g::CenterCategoryMorphism)
    dom = domain(f)⊗domain(g)
    cod = codomain(f)⊗codomain(g)
    return Morphism(dom,cod,f.m⊗g.m)
end

"""
    zero(C::CenterCategory)

Return the zero object of ```C```.
"""
function zero(C::CenterCategory)
    Z = zero(C.category)
    CenterCategoryObject(C,Z,[zero_morphism(Z,Z) for _ ∈ simples(C.category)])
end

"""
    one(C::CenterCategory)

Return the one object of ```C```.
"""
function one(C::CenterCategory)
    Z = one(C.category)
    CenterCategoryObject(C,Z,[id(s) for s ∈ simples(C.category)])
end


#-------------------------------------------------------------------------------
#   Is central?
#-------------------------------------------------------------------------------

"""
    is_central(Z::CategoryObject)

Return true if ```Z``` is in the categorical center, i.e. there exists a half-braiding on ```Z```.
"""
function is_central(Z::CategoryObject, simples::Vector{<:CategoryObject} = simples(parent(Z)))
    if prod([is_isomorphic(Z⊗s,s⊗Z)[1] for s ∈ simples]) == 0
        return false
    end
    return dim(build_center_ideal(Z,simples)) >= 0
end



function build_center_ideal(Z::CategoryObject, simples::Vector = simples(parent(Z)))
    @assert is_semisimple(parent(Z)) "Not semisimple"

    Homs = [Hom(Z⊗Xi, Xi⊗Z) for Xi ∈ simples]
    n = length(simples)
    ks = [dim(Homs[i]) for i ∈ 1:n]

    var_count = sum([int_dim(H) for H ∈ Homs])

    K = base_ring(Z)
    R,x = PolynomialRing(K, var_count, ordering = :lex)

    # For convinience: build arrays with the variables xi
    vars = []
    q = 1
    for i ∈ 1:n
        m = int_dim(Homs[i])
        vars = [vars; [x[q:q+m-1]]]
        q = q + m
    end

    eqs = []

    for k ∈ 1:n, i ∈ 1:n, j ∈ 1:n
        base = basis(Hom(Z⊗simples[k], simples[i]⊗(simples[j]⊗Z)))

        for t ∈ basis(Hom(simples[k], simples[i]⊗simples[j]))
            e = [zero(R) for i ∈ base]

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
    one_index = findfirst(e -> is_isomorphic(one(parent(Z)), e)[1], simples)
    one_c = K.(express_in_basis(id(Z), basis(End(Z))))
    push!(ideal_eqs, (vars[one_index] .- one_c)...)

    I = ideal([f for f ∈ unique(ideal_eqs) if f != 0])
end

function braidings_from_ideal(Z::CategoryObject, I::Ideal, simples::Vector{<:CategoryObject}, C)
    Homs = [Hom(Z⊗Xi, Xi⊗Z) for Xi ∈ simples]
    I = rational_lift(I)
    coeffs = recover_solutions(real_solutions(I),base_ring(Z))
    ks = [int_dim(H) for H ∈ Homs]
    centrals = CenterCategoryObject[]

    for c ∈ coeffs
        k = 1
        ex = CategoryMorphism[]
        c = [k for k ∈ c]
        for i ∈ 1:length(simples)
            if ks[i] == 0 continue end
            e = sum(c[k:k + ks[i] - 1] .* basis(Homs[i]))
            ex = [ex ; e]
            k = k + ks[i]
        end
        centrals = [centrals; CenterCategoryObject(C, Z, inv.(ex))]
    end
    return centrals
end

"""
    half_braidings(Z::CategoryObject)

Return all objects in the center lying over ```Z```.
"""
function half_braidings(Z::CategoryObject; simples = simples(parent(Z)), parent = Center(parent(Z)))

    I = build_center_ideal(Z,simples)

    d = dim(I)

    if d < 0 return CenterCategoryObject[] end

    if d == 0 return braidings_from_ideal(Z,I,simples, parent) end

    solutions = guess_solutions(Z,I,simples,CenterCategoryObject[],gens(base_ring(I)),d, parent)

    if length(solutions) == 0
        return CenterCategoryObject[]
    end
    unique_sols = solutions[1:1]

    for s ∈ solutions[2:end]
        if sum([dim(Hom(s,u)) for u ∈ unique_sols]) == 0
            unique_sols = [unique_sols; s]
        end
    end
    return unique_sols
end

function guess_solutions(Z::CategoryObject, I::Ideal, simples::Vector{<:CategoryObject}, solutions::Vector{CenterCategoryObject}, vars, d = dim(I), C = Center(parent(Z)))
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
    c_simples = CenterCategoryObject[]
    d_max = dim(C.category)
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

# function monoidal_completion(simples::Vector{CenterCategoryObject})
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

function isindependent(c::Vector,v::Vector...)
    if length(v) == 0 return true end
    m = matrix(ZZ, [vi[j] for vi ∈ v, j ∈ 1:length(v[1])])

    try
        x = solve(m,matrix(ZZ,c))
    catch
        return true
    end

    return !(*((x .>=0)...))
end

function find_centrals(simples::Vector{<:CategoryObject})
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
    braiding(X::CenterCategoryObject, Y::CenterCategoryObject)

Return the braiding isomorphism ```X⊗Y → Y⊗X```.
"""
function braiding(X::CenterCategoryObject, Y::CenterCategoryObject)
    dom = X.object⊗Y.object
    cod = Y.object⊗X.object
    # braid = zero_morphism(dom, cod)
    # for (s,ys) ∈ zip(simples(parent(X).category), X.γ)
    #     proj = basis(Hom(Y.object,s))
    #     if length(proj) == 0 continue end
    #     incl = basis(Hom(s,Y.object))
    #     braid = braid + sum([(i⊗id(X.object))∘ys∘(id(X.object)⊗p) for i ∈ incl, p ∈ proj][:])
    # end
    braid = half_braiding(X,object(Y))
    return Morphism(X⊗Y,Y⊗X,braid)
end

function half_braiding(X::CenterCategoryObject, Y::CategoryObject)
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
   
  
    iso, incl, proj = decompose_morphism(Y)

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
    dim(X::CenterCategoryObject)

Return the categorical dimension of ```X```.
"""
dim(X::CenterCategoryObject) = dim(X.object)

"""
    simples(C::CenterCategory)

Return a vector containing the simple objects of ```C```. The list might be incomplete.
"""
function simples(C::CenterCategory; sort = false)
    if isdefined(C, :simples) 
        # if dim(RingSubcategory(C.category,1))^2 != sum((dim.(C.simples)).^2)
        #     @warn "List not complete"
        # end
        return C.simples 
    end
    simples_by_induction!(C)
    if sort 
        sort_simples_by_dimension!(C)
    end
    return C.simples
end


"""
    associator(X::CenterCategoryObject, Y::CenterCategoryObject, Z::CenterCategoryObject)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
function associator(X::CenterCategoryObject, Y::CenterCategoryObject, Z::CenterCategoryObject)
    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)
    return Morphism(dom,cod, associator(X.object, Y.object, Z.object))
end

matrices(f::CenterCategoryMorphism) = matrices(f.m)
matrix(f::CenterCategoryMorphism) = matrix(f.m)

"""
    compose(f::CenterCategoryMorphism, g::CenterCategoryMorphism)

Return the composition ```g∘f```.
"""
compose(f::CenterCategoryMorphism, g::CenterCategoryMorphism) = Morphism(domain(f), codomain(g), g.m∘f.m)

"""
    dual(X::CenterCategoryObject)

Return the (left) dual object of ```X```.
"""
function dual(X::CenterCategoryObject)
    a = associator
    inv_a = inv_associator
    e = ev(X.object)
    c = coev(X.object)
    γ = CategoryMorphism[]
    dX = dual(X.object)
    for (Xi,yXi) ∈ zip(simples(parent(X).category), X.γ)
        f = (e⊗id(Xi⊗dX))∘inv_a(dX,X.object,Xi⊗dX)∘(id(dX)⊗a(X.object,Xi,dX))∘(id(dX)⊗(inv(yXi)⊗id(dX)))∘(id(dX)⊗inv_a(Xi,X.object,dX))∘a(dX,Xi,X.object⊗dX)∘(id(dX⊗Xi)⊗c)
        γ = [γ; f]
    end
    return CenterCategoryObject(parent(X),dX,γ)
end

"""
    ev(X::CenterCategoryObject)

Return the evaluation morphism ``` X⊗X → 1```.
"""
function ev(X::CenterCategoryObject)
    Morphism(dual(X)⊗X,one(parent(X)),ev(X.object))
end

"""
    coev(X::CenterCategoryObject)

Return the coevaluation morphism ```1 → X⊗X∗```.
"""
function coev(X::CenterCategoryObject)
    Morphism(one(parent(X)),X⊗dual(X),coev(X.object))
end

"""
    id(X::CenterCategoryObject)

Return the identity on ```X```.
"""
id(X::CenterCategoryObject) = Morphism(X,X,id(X.object))

"""
    tr(f:::CenterCategoryMorphism)

Return the categorical trace of ```f```.
"""
function tr(f::CenterCategoryMorphism)
    C = parent(domain(f))
    return CenterCategoryMorphism(one(C),one(C),tr(f.m))
end

"""
    inv(f::CenterCategoryMorphism)

Return the inverse of ```f```if possible.
"""
function inv(f::CenterCategoryMorphism)
    return Morphism(codomain(f),domain(f), inv(f.m))
end

"""
    is_isomorphic(X::CenterCategoryObject, Y::CenterCategoryObject)

Check if ```X≃Y```. Return ```(true, m)``` where ```m```is an isomorphism if true,
else return ```(false,nothing)```.
"""
function is_isomorphic(X::CenterCategoryObject, Y::CenterCategoryObject)
    # TODO: Fix This. How to compute a central isomorphism?
    S = simples(parent(X))

    if [dim(Hom(X,s)) for s ∈ S] == [dim(Hom(Y,s)) for s ∈ S]
        _, iso = is_isomorphic(X.object, Y.object)
        return true, Morphism(X,Y,central_projection(X,Y,iso))
    else
        return false, nothing
    end
end

function +(f::CenterCategoryMorphism, g::CenterCategoryMorphism)
    return Morphism(domain(f), codomain(f), g.m +f.m)
end

function *(x, f::CenterCategoryMorphism)
    return Morphism(domain(f),codomain(f),x*f.m)
end
#-------------------------------------------------------------------------------
#   Functionality: Image
#-------------------------------------------------------------------------------

"""
    kernel(f::CenterMoprhism)

Return a tuple ```(K,k)``` where ```K```is the kernel object and ```k```is the inclusion.
"""
function kernel(f::CenterCategoryMorphism)
    ker, incl = kernel(f.m)
    #f_inv = left_inverse(incl)

    braiding = [left_inverse(id(s)⊗incl)∘γ∘(incl⊗id(s)) for (s,γ) ∈ zip(simples(parent(domain(f.m))), domain(f).γ)]

    Z = CenterCategoryObject(parent(domain(f)), ker, braiding)
    return Z, Morphism(Z,domain(f), incl)
end

"""
    cokernel(f::CenterCategoryMorphism)

Return a tuple ```(C,c)``` where ```C```is the cokernel object and ```c```is the projection.
"""
function cokernel(f::CenterCategoryMorphism)
    coker, proj = cokernel(f.m)
    #f_inv = right_inverse(proj)

    braiding = [(id(s)⊗proj)∘γ∘(right_inverse(proj⊗id(s))) for (s,γ) ∈ zip(simples(parent(domain(f.m))), codomain(f).γ)]

    Z = CenterCategoryObject(parent(domain(f)), coker, braiding)
    return Z, Morphism(codomain(f),Z, proj)
end


function left_inverse(f::CenterCategoryMorphism)
    X = domain(f)
    Y = codomain(f)
    l_inv = central_projection(Y,X,left_inverse(morphsm(f)))
    return Morphism(Y,X,l_inv)
end

function quotient(Y::CenterCategoryObject, X::CategoryObject)
    # TODO: Compute quotient
    @assert parent(X) == parent(Y).Category
end

#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct CenterCategoryHomSpace <: AbstractCategoryHomSpace
    X::CenterCategoryObject
    Y::CenterCategoryObject
    basis::Vector{CenterCategoryMorphism}
    parent::VectorSpaces
end


Hom(X::CenterCategoryObject, Y::CenterCategoryObject) = hom_by_linear_equations(X,Y)

function central_projection(dom::CenterCategoryObject, cod::CenterCategoryObject, f::CategoryMorphism, simpls = simples(parent(domain(f))))
    X = domain(f)
    Y = codomain(f)
    C = parent(X)
    D = dim(C)
    proj = zero_morphism(X, Y)
    a = associator
    inv_a = inv_associator

    for (Xi, yX) ∈ zip(simpls, dom.γ)
        dXi = dual(Xi)

        yY = half_braiding(cod, dXi)
        
        ϕ = (ev(dXi)⊗id(Y))∘inv_a(dual(dXi),dXi,Y)∘(spherical(Xi)⊗yY)∘a(Xi,Y,dXi)∘((id(Xi)⊗f)⊗id(dXi))∘(yX⊗id(dXi))∘inv_a(X,Xi,dXi)∘(id(X)⊗coev(Xi))

        proj = proj + dim(Xi)*ϕ
    end
    return inv(D*base_ring(dom)(1))*proj
end

"""
    zero_morphism(X::CenterCategoryObject, Y::CenterCategoryObject)

Return the zero morphism ```0:X → Y```.
"""
zero_morphism(X::CenterCategoryObject, Y::CenterCategoryObject) = Morphism(X,Y,zero_morphism(X.object,Y.object))

#-------------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------------

function show(io::IO, X::CenterCategoryObject)
    print(io, "Central object: $(X.object)")
end

function show(io::IO, C::CenterCategory)
    print(io, "Drinfeld center of $(C.category)")
end

function show(io::IO, f::CenterCategoryMorphism)
    print(io, "Morphism in $(parent(domain(f)))")
end


#=------------------------------------------------
    Center by Induction
------------------------------------------------=#

function simples_by_induction!(C::CenterCategory)
    S = CenterCategoryObject[]
    d = dim(C.category)^2

    if characteristic(base_ring(C)) == 0 
        ordered_simples = sort(simples(C.category), by = fpdim)
    else 
        ordered_simples = simples(C.category)
    end

    FI_simples = induction_restriction.(ordered_simples)
    center_dim = 0
    for (s, Is) ∈ zip(ordered_simples, FI_simples)
       contained_simples = filter(x -> int_dim(Hom(object(x),s)) != 0, S)
        if length(contained_simples) > 0
            if is_isomorphic(Is, direct_sum(object.(contained_simples))[1])[1]
                continue
            end
        end

        Z = induction(s)
        # for x ∈ contained_simples
        #     f = horizontal_direct_sum(basis(Hom(x,Z)))
        #     Z = cokernel(f)[1]
        # end
        new_simples = indecomposable_subobjects(Z)
        S = [S; new_simples]
        center_dim += sum(dim.(new_simples).^2)
        if d == center_dim
            break
        end
    end
    C.simples = unique_simples(S)
end

function sort_simples_by_dimension!(C::CenterCategory)  
    fp_dims = [fpdim(s) for s ∈ simples(C)]
    K = base_ring(C)
    f = complex_embeddings(K)[1]
    σ = sortperm(fp_dims, by = e -> abs(f(e)))
    C.simples = C.simples[σ]
end


#=----------------------------------------------------------
    Hom Spaces 2.0 
----------------------------------------------------------=#

function hom_by_linear_equations(X::CenterCategoryObject, Y::CenterCategoryObject)
    @assert parent(X) == parent(Y)

    H = Hom(object(X), object(Y))
    B = basis(H)
    F = base_ring(X)
    n = length(basis(H))

    if n == 0 
        return CategoryHomSpace(X,Y, CenterCategoryMorphism[], VectorSpaces(F))
    end 

    Fx,poly_basis = PolynomialRing(F,n)
    
    eqs = []

    S = simples(parent(object(X)))

    for (s,γₛ,λₛ) ∈ zip(S,half_braiding(X), half_braiding(Y))

        Hs = Hom(object(X)⊗s, s⊗object(Y))
        base = basis(Hs)
        eq_i = [zero(Fx) for _ ∈ 1:length(base)]
        for (f,a) ∈ zip(B,poly_basis)
            coeffs = express_in_basis((id(s)⊗f)∘γₛ - λₛ ∘(f⊗id(s)), base)
            eq_i = eq_i .+ (a .* coeffs)
        end
        
        eqs = [eqs; eq_i]

    end

    M = zero(MatrixSpace(F,length(eqs),n))

    for (i,e) ∈ zip(1:length(eqs),eqs)
        M[i,:] = [coeff(e, a) for a ∈ poly_basis]
    end

    N = nullspace(M)[2]

    _,cols = size(N)

    basis_coeffs = [N[:,i] for i ∈ 1:cols]

    center_basis = [CenterCategoryMorphism(X,Y,sum(b .* B)) for b ∈ basis_coeffs]

    return CategoryHomSpace(X,Y,center_basis, VectorSpaces(F))
end

function hom_by_projection(X::CenterCategoryObject, Y::CenterCategoryObject)
    b = basis(Hom(X.object, Y.object))

    projs = [central_projection(X,Y,f) for f in b]

    proj_exprs = [express_in_basis(p,b) for p ∈ projs]

    M = zero(MatrixSpace(base_ring(X), length(b),length(b)))
    for i ∈ 1:length(proj_exprs)
        M[i,:] = proj_exprs[i]
    end
    r, M = rref(M)
    H_basis = CenterCategoryMorphism[]
    for i ∈ 1:r
        f = Morphism(X,Y,sum([m*bi for (m,bi) ∈ zip(M[i,:], b)]))
        H_basis = [H_basis; f]
    end
    return CenterCategoryHomSpace(X,Y,H_basis, VectorSpaces(base_ring(X)))
end

