#=----------------------------------------------------------
    An abstract structure to construct the arrow category
    Mor(𝒞) of any category 𝒞. 

    Objects are Morphisms f:X → Y ∈ 𝒞. Morphisms between 
    f:X→Y and g:A → B are pairs (ϕ:X→A,ψ:Y→B) such that 
    ψ∘f = g∘ϕ.
    
    If the category is abelian the arrow category is as 
    well. If 𝒞 is monoidal and has pushouts Mor(𝒞) is also
    monoidal with the pushout product.
----------------------------------------------------------=#

struct ArrowCategory <: Category
    category::Category
end

struct ArrowObject <: Object
    parent::ArrowCategory
    morphism::Morphism
end

struct ArrowMorphism <: Morphism
    domain::ArrowObject
    codomain::ArrowObject
    left::Morphism
    right::Morphism
end

domain(X::ArrowObject) = domain(morphism(X))
codomain(X::ArrowObject) = codomain(morphism(X))
morphism(X::ArrowObject) = X.morphism

left(f::ArrowMorphism) = f.left
right(f::ArrowMorphism) = f.right

function compose(f::ArrowMorphism...)
    Morphism(domain(f[1]), codomain(f[end]), compose(left.(f)...), compose(right.(f)...))
end

function id(X::ArrowObject) 
    Morphism(X,X, id(domain(X)), id(codomain(X)))
end

#=----------------------------------------------------------
    Abelian structure 
----------------------------------------------------------=#
is_abelian(C::ArrowCategory) = is_abelian(category(C))

matrix(f::ArrowMorphism) = diagonal_matrix(matrix(f.left), matrix(f.right))

function direct_sum(X::ArrowObject...)
    dom, dom_incl, dom_proj = direct_sum(domain.(X))
    cod, cod_incl, cod_proj = direct_sum(codomain.(X))

    S = ArrowObject(parent(X[1]), direct_sum(morphism.(X)))

    incl = [ArrowMorphism(f, S, i_d, i_c) for (f,i_d,i_c) ∈ zip(X, dom_incl, cod_incl)]

    proj = [ArrowMorphism(S, f, p_d, p_c) for (f,p_d,p_c) ∈ zip(X, dom_proj, cod_proj)]

    return S, incl, proj
end

function direct_sum(f::ArrowMorphism...)
    dom = direct_sum(domain.(f))[1]
    cod = direct_sum(codomain.(f))[1]

    left = direct_sum([g.left for g ∈ f])
    right = direct_sum([g.right for g ∈ f])

    return Morphism(dom, cod, left, right)
end

function *(λ, f::ArrowMorphism)
    Morphism(domain(f),codomain(f), λ * f.left, λ * f.right)
end

function +(f::ArrowMorphism, g::ArrowMorphism)
    Morphism(domain(f), codomain(f), f.left + g.left, f.right + g.right)
end

function kernel(f::ArrowMorphism)
    _, k_left = kernel(f.left)
    _, k_right = kernel(f.right)

    K = ArrowObject(parent(f), left_inverse(k_right) ∘ morphism(domain(f)) ∘ k_left)
    incl = Morphism(K, domain(f), k_left, k_right)

    return K, incl
end

function cokernel(f::ArrowMorphism)
    _, c_left = cokernel(f.left)
    _, c_right = cokernel(f.right)

    C = ArrowObject(parent(f), c_right ∘ morphism(codomain(f)) ∘ right_inverse(c_left))
    proj = Morphism(codomain(f), C, c_left, c_right)

    return C, proj
end

function zero(C::ArrowCategory)
    ArrowObject(C, zero_morphism(category(C)))
end

function zero_morphism(X::ArrowObject, Y::ArrowObject)
    ArrowMorphism(X,Y, zero_morphism(domain(X),domain(Y)), zero_morphism(codomain(X), codomain(Y)))
end

function simples(C::ArrowCategory)
    S = ArrowObject[]
    Z = zero(category(C))
    for s ∈ simples(category(C))
        push!(S, ArrowObject(C, zero_morphism(Z, s)), 
                    ArrowObject(C, zero_morphism(s,Z)))
    end
    return S
end

function indecomposables(C::ArrowCategory)
    S = ArrowObject[]
    Z = zero(category(C))

    indecs = indecomposables(category(C))

    for s ∈ indecs
        push!(S, ArrowObject(C, zero_morphism(Z, s)), 
                    ArrowObject(C, zero_morphism(s,Z)))

        for t ∈ indecs
            for f ∈ Hom(s,t)
                push!(S, ArrowObject(C, f))
            end
        end
    end 
    return S
end

# function is_isomorphic(X::ArrowObject, Y::ArrowObject)
#     SX, LX, RX = snf_with_transform(matrix(morphism(X)))
#     SY, LY, RY = snf_with_transform(matrix(morphism(Y)))

#     if SX != SY 
#         return false, nothing
#     end

#     return true, Morphism(X,Y, )
#=----------------------------------------------------------
    Monoidal structure 
----------------------------------------------------------=#

function tensor_product(X::ArrowObject, Y::ArrowObject)
    Z = pushout_product(morphism(X), morphism(Y))
    ArrowObject(parent(X), Z)
end

function tensor_product(f::ArrowMorphism, g::ArrowMorphism)
    cod_f = codomain(f)
    cod_g = codomain(g)
    dom_f = domain(f)
    dom_g = domain(g)

    dom = domain(f) ⊗ (domain(g))
    cod = codomain(f) ⊗ (codomain(g))
    base = basis(Hom(domain(dom), domain(cod)))

    if length(base) == 0 
        return Morphism(dom, cod, zero_morphism(domain(dom), domain(cod)), right(f)⊗right(g))
    end

    _,(cx,cy) = pushout(morphism(cod_f) ⊗ id(domain(cod_g)), id(domain(cod_f)) ⊗ morphism(cod_g))
    
    _,(dx,dy) = pushout(morphism(dom_f) ⊗ id(domain(dom_g)), id(domain(dom_f)) ⊗ morphism(dom_g)) 



    mor_1 = cx ∘ (right(f) ⊗ left(g))
    mor_2 = cy ∘ (left(f) ⊗ right(g))

    K = base_ring(f)
    
    base_1 = basis(Hom(domain(mor_1), codomain(mor_1)))
    base_2 = basis(Hom(domain(mor_2), codomain(mor_2)))
    n = length(base_1) + length(base_2)

    Rx,x = PolynomialRing(K, length(base))

    eqs = [zero(Rx) for _ ∈ 1:n]

    for (h,a) ∈ zip(base, x)
        e_1 = express_in_basis(h ∘ dx, base_1)
        e_2 = express_in_basis(h ∘ dy, base_2)

        eqs = eqs .+ (a.* [e_1; e_2])
    end

    M_arr = hcat([[coeff(e, a) for a ∈ x] for e ∈ eqs]...)
    b_arr = [express_in_basis(mor_1, base_1); express_in_basis(mor_2, base_1)]

    M = matrix(K, length(base), length(eqs), M_arr)
    b = matrix(K, 1, length(eqs), b_arr)

    s = solve_left(M,b)

    l = sum(collect(s)[:] .* base)
    
    Morphism(dom,cod, l, right(f) ⊗ right(g))
end

one(C::ArrowCategory) = ArrowObject(C, zero_morphism(zero(category(C)), one(category(C))))

function associator(X::ArrowObject, Y::ArrowObject, Z::ArrowObject)
    ass_right = associator(codomain.((X,Y,Z))...)

    
    Morphism((X⊗Y)⊗Z, X⊗(Y⊗Z), ass_left, ass_right)
end
#=----------------------------------------------------------
    Hom spaces 
----------------------------------------------------------=#

function Hom(X::ArrowObject, Y::ArrowObject)
    base = basis(Hom(domain(X), codomain(Y)))

    base_dom = basis(Hom(domain(X), domain(Y)))
    base_cod = basis(Hom(codomain(X), codomain(Y)))
    n,m = length(base_dom), length(base_cod)

    if n+m == 0 
        return HomSpace(X,Y, ArrowMorphism[], VectorSpaces(base_ring(X)))
    end

    F = base_ring(X)

    Rx,x = PolynomialRing(F, n+m)

    eqs = [zero(Rx) for _ ∈ length(base)]

    mX = morphism(X)
    mY = morphism(Y)

    for (f,a) ∈ zip(base_dom, x[1:n])
        eqs = eqs .+ (a .* express_in_basis(mY ∘ f, base))
    end

    for (g,a) ∈ zip(base_cod, x[n+1:n+m])
        eqs = eqs .- (a .* express_in_basis(g ∘ mX, base))
    end

    M_arr = [coeff(e, a) for e ∈ eqs, a ∈ x]
    M = matrix(F, length(eqs), length(x), M_arr)

    N = collect(nullspace(M)[2])

    sols = [collect(s) for s ∈ collect(eachcol(N))]

    if length(base_dom) == 0
        B = [(zero_morphism(domain(X), domain(Y)), sum(s[n+1:n+m] .* base_cod)) for s ∈ sols]
    elseif length(base_cod) == 0 
        B = [(sum(s[1:n] .* base_dom), zero_morphism(codomain(X), codomain(Y))) for s ∈ sols]
    else
        B = [(sum(s[1:n] .* base_dom), sum(s[n+1:n+m] .* base_cod)) for s ∈ sols]
    end

    B = ArrowMorphism[Morphism(X,Y, l, r) for (l,r) ∈ B]

    return HomSpace(X,Y, unique_without_hash(B), VectorSpaces(F))
end
    

#=----------------------------------------------------------
    pretty print 
----------------------------------------------------------=#

function show(io::IO, C::ArrowCategory)
    print(io, """Category of morphisms in $(category(C))""")
end

function show(io::IO, X::ArrowObject)
    print(io, """Arrow object: $(domain(X)) → $(codomain(X))""")
end

function show(io::IO, f::ArrowMorphism)
    print(io, """Morphism in arrow category defined by 
    $(f.left)
    and
    $(f.right)""")
end