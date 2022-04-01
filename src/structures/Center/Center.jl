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

function Center(C::Category)
    @assert issemisimple(C) "Semisimplicity required"
    return CenterCategory(base_ring(C), C)
end

function Morphism(dom::CenterObject, cod::CenterObject, m::Morphism)
    return CenterMorphism(dom,cod,m)
end

#-------------------------------------------------------------------------------
#   Direct Sum & Tensor Product
#-------------------------------------------------------------------------------

function dsum(X::CenterObject, Y::CenterObject)
    Z = X.object ⊕ Y.object
    γZ = [y ⊕ q for (y,q) ∈ zip(X.γ, Y.γ)]
    return CenterObject(parent(X), Z, γZ)
end

function dsum(f::CenterMorphism, g::CenterMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    m = f.m ⊕ g.m
    return Morphism(dom,cod, m)
end

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
#-------------------------------------------------------------------------------
#   Induction
#-------------------------------------------------------------------------------

function induction(X::Object, simples::Vector = simples(parent(X)))
    @assert issemisimple(parent(X)) "Requires semisimplicity"
    Z = dsum([dual(s)⊗X⊗s for s ∈ simples])

    function γ(W)
        r = Morphism[]
        for i ∈ simples, j ∈ simples
            b1 = basis(Hom(W⊗dual(i),j))
            b2 = basis(Hom(i,j⊗W))
            if length(b1)*length(b2) == 0 continue end
            push!(r,dim(i)*dsum([ϕ ⊗ id(X) ⊗ ψ for (ϕ,ψ) ∈ zip(b1,b2)]))
        end
        return dsum(r)
    end
    return CenterObject(CenterCategory(base_ring(X),parent(X)),Z,γ)
end



#-------------------------------------------------------------------------------
#   Is central?
#-------------------------------------------------------------------------------


function iscentral(Z::Object, simples::Vector = simples(parent(Z)))
    @assert issemisimple(parent(Z)) "Not semisimple"

    if dim(Z) == 0 return true end

    Homs = [Hom(Z⊗Xi, Xi⊗Z) for Xi ∈ simples]
    n = length(simples)
    ks = [length(basis(Homs[i])) for i ∈ 1:n]

    if prod(ks) == 0 return false, nothing end

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
        base = basis(Hom(Z⊗simples[k], simples[i]⊗simples[j]⊗Z))

        for t ∈ basis(Hom(simples[k], simples[i]⊗simples[j]))
            e = [zero(R) for i ∈ base]

            # for ai ∈ 1:dim(Homs[k]), bi ∈ 1:dim(Homs[j]), ci ∈ 1:dim(Homs[i])
            #     a,b,c = basis(Homs[k])[ai], basis(Homs[j])[bi], basis(Homs[i])[ci]
            #     l1 = express_in_basis((t⊗id(Z))∘a, base)
            #     l2 = express_in_basis((id(simples[i])⊗b)∘associator(simples[i],Z,simples[j]) ∘ (c⊗id(simples[j])) ∘ inv(associator(Z,simples[i],simples[j])) ∘ (id(Z) ⊗ t), base)
            #     e = e .+ (vars[k][ai] .* l1) .- (vars[j][bi]*vars[i][ci] .* l2)
            # end
            l1 = [zero(R) for i ∈ base]
            l2 = [zero(R) for i ∈ base]

            for ai ∈ 1:dim(Homs[k])
                a = basis(Homs[k])[ai]
                l1 = l1 .+ vars[k][ai] .* QQ.(express_in_basis(associator(simples[i],simples[j],Z)∘(t⊗id(Z))∘a, base))
            end
            for bi ∈ 1:dim(Homs[j]), ci ∈ 1:dim(Homs[i])
                b,c = basis(Homs[j])[bi], basis(Homs[i])[ci]
                l2 = l2 .+ (vars[j][bi]*vars[i][ci]) .* QQ.(express_in_basis((id(simples[i])⊗b)∘associator(simples[i],Z,simples[j]) ∘ (c⊗id(simples[j])) ∘ inv(associator(Z,simples[i],simples[j])) ∘ (id(Z) ⊗ t), base))
            end
            push!(eqs, l1 .-l2)
        end
    end
    ideal_eqs = []
    for p ∈ eqs
        push!(ideal_eqs, p...)
    end

    #normalize isomorphisms
    I = ideal([f for f ∈ unique(ideal_eqs) if f != 0])

    # for i ∈ 1:dim(I)+1
    #     push!(ideal_eqs, sum(vars[i] .^2) - dim(Homs[i]))
    # end
    #Require e_Z(1) = id(Z)
    one_c = QQ.(express_in_basis(id(Z), basis(End(Z))))
    push!(ideal_eqs, (vars[1] .- one_c)...)

    I = ideal([f for f ∈ unique(ideal_eqs) if f != 0])

    if dim(I) < 0 return false,nothing end

    coeffs = recover_solutions(msolve(I),base_ring(Z))

    centrals = CenterObject{elem_type(base_ring(Z))}[]

    C = Center(parent(Z))

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
    return true, centrals
end


function center_simples(C::Category, simples = simples(C))
    d = dim(C)^2

    simples_indices = []
    c_simples = typeof(simples[1])[]
    d_max = d
    d_rem = d
    k = length(simples)

    coeffs = [i for i ∈ Base.product([0:d_max for i ∈ 1:k]...) if (sum(i) != 0 && sum((i .* dim.(simples)).^2) <= d_rem)][:]

    for c ∈ sort(coeffs, by = t -> (sum(t),length(t) - length([i for i ∈ t if i != 0])))
        if sum((c .* dim.(simples)).^2) > d_rem continue end
        if simples_covered(c,simples_indices) continue end
        X = dsum([simples[j]^c[j] for j ∈ 1:k])

        ic, so = iscentral(X)

        if ic
            c_simples = [c_simples; so]
            d_rem = d_rem - sum([dim(x)^2 for x in so])
            if d_rem == 0 return c_simples end
            push!(simples_indices, c)
        end
    end

    return c_simples
end

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


function build_half_braiding(Z::Object, Y::Object, ex::Vector, simples::Vector = simples(parent(Z)))

end

#-------------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------------

dim(X::CenterObject) = dim(X.object)

function simples(C::CenterCategory)
    if isdefined(C, :simples) return C.simples end
    C.simples = center_simples(C.category)
    return C.simples
end

function associator(X::CenterObject, Y::CenterObject, Z::CenterObject)
    dom = (X⊗Y)⊗Z
    cod = X⊗(Y⊗Z)
    return Morphism(dom,cod, associator(X.object, Y.object, Z.object))
end

matrices(f::CenterMorphism) = matrices(f.m)
matrix(f::CenterMorphism) = matrix(f.m)
#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct CenterHomSpace <: HomSpace
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
        #index of dual
        dualXi = findfirst(x -> isisomorphic(dual(Xi),x)[1], simples)
        yY = cod.γ[dualXi]
        dXi = dual(Xi)
        ϕ = (ev(dXi)⊗id(Y))∘inv(a(dual(dXi),dXi,Y))∘(spherical(Xi)⊗yY)∘a(Xi,Y,dXi)∘((id(Xi)⊗f)⊗id(dXi))∘(yX⊗id(dXi))∘inv(a(X,Xi,dXi))∘(id(X)⊗coev(Xi))

        proj = proj + dim(Xi)*ϕ
    end
    return inv(D*base_ring(dom)(1))*proj
end


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
