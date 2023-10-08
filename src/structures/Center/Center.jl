mutable struct CenterCategory <: Category
    base_ring::Field
    category::Category
    simples::Vector{O} where O <: Object

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

struct CenterObject <: Object
    parent::CenterCategory
    object::Object
    γ::Vector{M} where M <: Morphism
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

function Morphism(dom::CenterObject, cod::CenterObject, m::Morphism)
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

Return the onderlying object in ```𝒞```.
"""
object(X::CenterObject) = X.object

morphism(f::CenterMorphism) = f.m

is_fusion(C::CenterCategory) = true

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
spherical(X::CenterObject) = Morphism(X,dual(dual(X)), spherical(X.object))

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
    return Morphism(dom,cod, m)
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
    return Morphism(dom,cod,f.m⊗g.m)
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



function build_center_ideal(Z::Object, simples::Vector = simples(parent(Z)))
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
        centrals = [centrals; CenterObject(C, Z, inv.(ex))]
    end
    return centrals
end

"""
    half_braidings(Z::Object)

Return all objects in the center lying over ```Z```.
"""
function half_braidings(Z::Object; simples = simples(parent(Z)), parent = Center(parent(Z)))

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

function guess_solutions(Z::Object, I::Ideal, simples::Vector{<:Object}, solutions::Vector{CenterObject}, vars, d = dim(I), C = Center(parent(Z)))
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

Return the braiding isomorphism ```X⊗Y → Y⊗X```.
"""
function braiding(X::CenterObject, Y::CenterObject)
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
    simples(C::CenterCategory)

Return a vector containing the simple objects of ```C```. The list might be incomplete.
"""
function simples(C::CenterCategory; sort = false)
    if isdefined(C, :simples) 
        return C.simples 
    end
    if is_modular(category(C))
        C.simples = center_simples_by_braiding(category(C), C)
        return C.simples
    end
    simples_by_induction!(C)
    if sort 
        sort_simples_by_dimension!(C)
    end
    return C.simples
end

function decompose(X::CenterObject)
    C = parent(X)
    if isdefined(C, :simples)
        return decompose_by_simples(X,simples(C))
    else
        return decompose_by_endomorphism_ring(X)
    end
end

function decompose(X::CenterObject, S::Vector{CenterObject})
    decompose_by_simples(X,S)
end

"""
    associator(X::CenterObject, Y::CenterObject, Z::CenterObject)

Return the associator isomorphism ```(X⊗Y)⊗Z → X⊗(Y⊗Z)```.
"""
function associator(X::CenterObject, Y::CenterObject, Z::CenterObject)
    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)
    return Morphism(dom,cod, associator(X.object, Y.object, Z.object))
end

matrices(f::CenterMorphism) = matrices(f.m)
matrix(f::CenterMorphism) = matrix(f.m)

"""
    compose(f::CenterMorphism, g::CenterMorphism)

Return the composition ```g∘f```.
"""
compose(f::CenterMorphism, g::CenterMorphism) = Morphism(domain(f), codomain(g), g.m∘f.m)

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
        f = (e⊗id(Xi⊗dX))∘inv_a(dX,X.object,Xi⊗dX)∘(id(dX)⊗a(X.object,Xi,dX))∘(id(dX)⊗(inv(yXi)⊗id(dX)))∘(id(dX)⊗inv_a(Xi,X.object,dX))∘a(dX,Xi,X.object⊗dX)∘(id(dX⊗Xi)⊗c)
        γ = [γ; f]
    end
    return CenterObject(parent(X),dX,γ)
end

"""
    ev(X::CenterObject)

Return the evaluation morphism ``` X⊗X → 1```.
"""
function ev(X::CenterObject)
    Morphism(dual(X)⊗X,one(parent(X)),ev(X.object))
end

"""
    coev(X::CenterObject)

Return the coevaluation morphism ```1 → X⊗X∗```.
"""
function coev(X::CenterObject)
    Morphism(one(parent(X)),X⊗dual(X),coev(X.object))
end

"""
    id(X::CenterObject)

Return the identity on ```X```.
"""
id(X::CenterObject) = Morphism(X,X,id(X.object))

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
    return Morphism(codomain(f),domain(f), inv(f.m))
end


"""
    is_isomorphic(X::CenterObject, Y::CenterObject)

Check if ```X≃Y```. Return ```(true, m)``` where ```m```is an isomorphism if true,
else return ```(false,nothing)```.
"""
function is_isomorphic(X::CenterObject, Y::CenterObject)
    # TODO: Fix This. How to compute a central isomorphism?

    if is_simple(X) && is_simple(Y)
        H = Hom(X,Y)
        if int_dim(H) > 0
            return true, basis(H)[1]
        else
            return false, nothing
        end
    end

    S = simples(parent(X))

    if [dim(Hom(X,s)) for s ∈ S] == [dim(Hom(Y,s)) for s ∈ S]
        _, iso = is_isomorphic(X.object, Y.object)
        return true, Morphism(X,Y,central_projection(X,Y,iso))
    else
        return false, nothing
    end
end

function +(f::CenterMorphism, g::CenterMorphism)
    return Morphism(domain(f), codomain(f), g.m +f.m)
end

function *(x, f::CenterMorphism)
    return Morphism(domain(f),codomain(f),x*f.m)
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

    braiding = [left_inverse(id(s)⊗incl)∘γ∘(incl⊗id(s)) for (s,γ) ∈ zip(simples(parent(domain(f.m))), domain(f).γ)]

    Z = CenterObject(parent(domain(f)), ker, braiding)
    return Z, Morphism(Z,domain(f), incl)
end

"""
    cokernel(f::CenterMorphism)

Return a tuple ```(C,c)``` where ```C```is the cokernel object and ```c```is the projection.
"""
function cokernel(f::CenterMorphism)
    coker, proj = cokernel(f.m)
    #f_inv = right_inverse(proj)

    braiding = [(id(s)⊗proj)∘γ∘(right_inverse(proj⊗id(s))) for (s,γ) ∈ zip(simples(parent(domain(f.m))), codomain(f).γ)]

    Z = CenterObject(parent(domain(f)), coker, braiding)
    return Z, Morphism(codomain(f),Z, proj)
end


function left_inverse(f::CenterMorphism)
    X = domain(f)
    Y = codomain(f)
    l_inv = central_projection(Y,X,left_inverse(morphsm(f)))
    return Morphism(Y,X,l_inv)
end

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


Hom(X::CenterObject, Y::CenterObject) = hom_by_linear_equations(X,Y)

function central_projection(dom::CenterObject, cod::CenterObject, f::Morphism, simpls = simples(parent(domain(f))))
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
    zero_morphism(X::CenterObject, Y::CenterObject)

Return the zero morphism ```0:X → Y```.
"""
zero_morphism(X::CenterObject, Y::CenterObject) = Morphism(X,Y,zero_morphism(X.object,Y.object))

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

function simples_by_induction!(C::CenterCategory)
    S = CenterObject[]
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
        for x ∈ contained_simples
            f = horizontal_direct_sum(basis(Hom(x,Z)))
            Z = cokernel(f)[1]
        end
        new_simples = indecomposable_subobjects(Z)
        S = [S; new_simples]
        center_dim += sum(dim.(new_simples).^2)
        # if d == center_dim
        #     break
        # end
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

function hom_by_linear_equations(X::CenterObject, Y::CenterObject)
    #@assert parent(X) == parent(Y)

    H = Hom(object(X), object(Y))
    B = basis(H)
    F = base_ring(X)
    n = length(basis(H))

    if n == 0 
        return HomSpace(X,Y, CenterMorphism[], VectorSpaces(F))
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

    center_basis = [CenterMorphism(X,Y,sum(b .* B)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,center_basis, VectorSpaces(F))
end

function hom_by_projection(X::CenterObject, Y::CenterObject)
    b = basis(Hom(X.object, Y.object))

    projs = [central_projection(X,Y,f) for f in b]

    proj_exprs = [express_in_basis(p,b) for p ∈ projs]

    M = zero(MatrixSpace(base_ring(X), length(b),length(b)))
    for i ∈ 1:length(proj_exprs)
        M[i,:] = proj_exprs[i]
    end
    r, M = rref(M)
    H_basis = CenterMorphism[]
    for i ∈ 1:r
        f = Morphism(X,Y,sum([m*bi for (m,bi) ∈ zip(M[i,:], b)]))
        H_basis = [H_basis; f]
    end
    return CenterHomSpace(X,Y,H_basis, VectorSpaces(base_ring(X)))
end


#=----------------------------------------------------------
    Modular Stuff 
----------------------------------------------------------=#    

function smatrix(C::CenterCategory)
    simpls = simples(C)
    n = length(simpls)
    K = base_ring(C)
    S = [zero_morphism(category(C)) for _ ∈ 1:n, _ ∈ 1:n]
    for i ∈ 1:n
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

#=----------------------------------------------------------
    extension_of_scalars 
----------------------------------------------------------=#    

function extension_of_scalars(C::CenterCategory, L::Field)
    CL = _extension_of_scalars(C,L, category(C)⊗L)
    
    CL.simples = [extension_of_scalars(s, L, CL) for s ∈ simples(C)]

    return CL
end

function _extension_of_scalars(C::CenterCategory, L::Field, cL = category(C)⊗L)
    CenterCategory(L,cL)
end

function extension_of_scalars(X::CenterObject, L::Field, CL = _extension_of_scalars(parent(X),L))
    CenterObject(CL, extension_of_scalars(object(X), L, category(CL)), [f ⊗ L for f ∈ half_braiding(X)])
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

function center_simples_by_braiding(C::Category, Z = Center(C))
    S = simples(C)

    S_braided = [CenterObject(Z, s, [braiding(s,t) for t ∈ S]) for s ∈ S]
    S_rev_braided = [CenterObject(Z, s, [inv(braiding(t,s)) for t ∈ S]) for s ∈ S]

    [s⊗t for s ∈ S_braided, t ∈ S_rev_braided][:]
end
