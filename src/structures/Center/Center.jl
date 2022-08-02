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



#-------------------------------------------------------------------------------
#   Center Constructor
#-------------------------------------------------------------------------------
"""
    Center(C::Category)

Return the Drinfeld center of ```C```.
"""
function Center(C::Category)
    @assert issemisimple(C) "Semisimplicity required"
    return CenterCategory(base_ring(C), C)
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

isfusion(C::CenterCategory) = true

"""
    add_simple!(C::CenterCategory, S::CenterObject)

Add the simple object ```S``` to the vector of simple objects.
"""
function add_simple!(C::CenterCategory, S::CenterObject)
    @assert dim(End(S)) == 1 "Not simple"
    C.simples = unique_simples([simples(C); S])
end

"""
    spherical(X::CenterObject)

Return the spherical structure ```X → X∗∗``` of ```X```.
"""
spherical(X::CenterObject) = Morphism(X,dual(dual(X)), spherical(X.object))

(F::Field)(f::CenterMorphism) = F(f.m)
#-------------------------------------------------------------------------------
#   Direct Sum & Tensor Product
#-------------------------------------------------------------------------------

"""
    dsum(X::CenterObject, Y::CenterObject)

Return the direct sum object of ```X``` and ```Y```.
"""
function dsum(X::CenterObject, Y::CenterObject)
    S = simples(parent(X.object))
    Z,(ix,iy),(px,py) = dsum(X.object, Y.object,true)

    γZ = [(id(S[i])⊗ix)∘(X.γ[i])∘(px⊗id(S[i])) + (id(S[i])⊗iy)∘(Y.γ[i])∘(py⊗id(S[i])) for i ∈ 1:length(S)]
    return CenterObject(parent(X), Z, γZ)
end

"""
    dsum(f::CenterMorphism, g::CenterMorphism)

Return the direct sum of ```f``` and ```g```.
"""
function dsum(f::CenterMorphism, g::CenterMorphism)
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
    a = associator
    s = simples(parent(X.object))
    x,y = X.object, Y.object
    for (S, yX, yY) ∈ zip(s, X.γ, Y.γ)
        push!(γ, a(S,x,y)∘(yX⊗id(y))∘inv(a(x,S,y))∘(id(x)⊗yY)∘a(x,y,S))
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
    iscentral(Z::Object)

Return true if ```Z``` is in the categorical center, i.e. there exists a half-braiding on ```Z```.
"""
function iscentral(Z::Object, simples::Vector{<:Object} = simples(parent(Z)))
    if prod([isisomorphic(Z⊗s,s⊗Z)[1] for s ∈ simples]) == 0
        return false
    end
    return dim(build_center_ideal(Z,simples)) >= 0
end



function build_center_ideal(Z::Object, simples::Vector = simples(parent(Z)))
    @assert issemisimple(parent(Z)) "Not semisimple"

    Homs = [Hom(Z⊗Xi, Xi⊗Z) for Xi ∈ simples]
    n = length(simples)
    ks = [dim(Homs[i]) for i ∈ 1:n]

    var_count = sum([dim(H) for H ∈ Homs])

    R,x = PolynomialRing(QQ, var_count, ordering = :lex)

    # For convinience: build arrays with the variables xi
    vars = []
    q = 1
    for i ∈ 1:n
        m = dim(Homs[i])
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

            for ai ∈ 1:dim(Homs[k])
                a = basis(Homs[k])[ai]
                l1 = l1 .+ (vars[k][ai] .* QQ.(express_in_basis(associator(simples[i],simples[j],Z)∘(t⊗id(Z))∘a, base)))
            end
            for bi ∈ 1:dim(Homs[j]), ci ∈ 1:dim(Homs[i])
                b,c = basis(Homs[j])[bi], basis(Homs[i])[ci]
                l2 = l2 .+ ((vars[j][bi]*vars[i][ci]) .* QQ.(express_in_basis((id(simples[i])⊗b)∘associator(simples[i],Z,simples[j]) ∘ (c⊗id(simples[j])) ∘ inv(associator(Z,simples[i],simples[j])) ∘ (id(Z) ⊗ t), base)))
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
    one_index = findfirst(e -> isisomorphic(one(parent(Z)), e)[1], simples)
    one_c = QQ.(express_in_basis(id(Z), basis(End(Z))))
    push!(ideal_eqs, (vars[one_index] .- one_c)...)

    I = ideal([f for f ∈ unique(ideal_eqs) if f != 0])
end

function braidings_from_ideal(Z::Object, I::Ideal, simples::Vector{<:Object}, C)
    Homs = [Hom(Z⊗Xi, Xi⊗Z) for Xi ∈ simples]
    coeffs = recover_solutions(msolve(I),base_ring(Z))
    ks = [dim(H) for H ∈ Homs]
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
        centrals = [centrals; CenterObject(C, Z, ex)]
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

        X = dsum([simples[j]^c[j] for j ∈ 1:k])

        ic = iscentral(X)

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
        ic, so = iscentral(s)
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
    braid = zero_morphism(dom, cod)
    for (s,ys) ∈ zip(simples(parent(X).category), X.γ)
        proj = basis(Hom(Y.object,s))
        if length(proj) == 0 continue end
        incl = basis(Hom(s,Y.object))
        braid = braid + sum([(i⊗id(X.object))∘ys∘(id(X.object)⊗p) for i ∈ incl, p ∈ proj][:])
    end
    return Morphism(X⊗Y,Y⊗X,braid)
end

function half_braiding(X::CenterObject, Y::Object)
    dom = X.object⊗Y
    cod = Y⊗X.object
    braid = zero_morphism(dom, cod)
    simpls = simples(parent(Y))

    iso, incl, proj = decompose_morphism(Y)

    for (p,i) ∈ zip(proj, incl)
        k = findfirst(e -> isisomorphic(e, domain(i))[1], simpls)
        incliso = isisomorphic(simpls[k], domain(i))[2]
        #projiso = isisomorphic(codomain(p), )[2]

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
function simples(C::CenterCategory)
    if isdefined(C, :simples) return C.simples end
    C.simples = center_simples(C)
    return C.simples
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
    e = ev(X.object)
    c = coev(X.object)
    γ = Morphism[]
    dX = dual(X.object)
    for (Xi,yXi) ∈ zip(simples(parent(X).category), X.γ)
        f = (e⊗id(Xi⊗dX))∘inv(a(dX,X.object,Xi⊗dX))∘(id(dX)⊗a(X.object,Xi,dX))∘(id(dX)⊗(inv(yXi)⊗id(dX)))∘(id(dX)⊗inv(a(Xi,X.object,dX)))∘a(dX,Xi,X.object⊗dX)∘(id(dX⊗Xi)⊗c)
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
    isisomorphic(X::CenterObject, Y::CenterObject)

Check if ```X≃Y```. Return ```(true, m)``` where ```m```is an isomorphism if true,
else return ```(false,nothing)```.
"""
function isisomorphic(X::CenterObject, Y::CenterObject)
    S = simples(parent(X))
    if [dim(Hom(X,s)) for s ∈ S] == [dim(Hom(Y,s)) for s ∈ S]
        return true, inv(decompose_morphism(Y))∘decompose_morphism(X)
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
    f_inv = left_inverse(incl)

    braiding = [(id(s)⊗f_inv)∘γ∘(incl⊗id(s)) for (s,γ) ∈ zip(simples(parent(domain(f.m))), domain(f).γ)]

    Z = CenterObject(parent(domain(f)), ker, braiding)
    return Z, Morphism(Z,domain(f), incl)
end

"""
    cokernel(f::CenterMoprhism)

Return a tuple ```(C,c)``` where ```C```is the cokernel object and ```c```is the projection.
"""
function cokernel(f::CenterMorphism)
    coker, proj = cokernel(f.m)
    f_inv = right_inverse(proj)

    braiding = [(proj⊗id(s))∘γ∘(id(s)⊗f_inv) for (s,γ) ∈ zip(simples(parent(domain(f.m))), codomain(f).γ)]

    Z = CenterObject(parent(domain(f)), coker, braiding)
    return Z, Morphism(codomain(f),Z, proj)
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

function central_projection(dom::CenterObject, cod::CenterObject, f::Morphism, simples = simples(parent(domain(f))))
    X = domain(f)
    Y = codomain(f)
    C = parent(X)
    D = dim(C)
    proj = zero_morphism(X, Y)
    a = associator

    for (Xi, yX) ∈ zip(simples, dom.γ)
        dXi = dual(Xi)
        yY = half_braiding(cod, dXi)

        @show (ev(dXi)⊗id(Y))∘inv(a(dual(dXi),dXi,Y))∘(spherical(Xi)⊗yY)∘a(Xi,Y,dXi)∘((id(Xi)⊗f)⊗id(dXi))∘(yX⊗id(dXi))
        ϕ = (ev(dXi)⊗id(Y))∘inv(a(dual(dXi),dXi,Y))∘(spherical(Xi)⊗yY)∘a(Xi,Y,dXi)∘((id(Xi)⊗f)⊗id(dXi))∘(yX⊗id(dXi))∘inv(a(X,Xi,dXi))∘(id(X)⊗coev(Xi))

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
