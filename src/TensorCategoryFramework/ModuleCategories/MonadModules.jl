#=----------------------------------------------------------
    Modules over a Monad T : 𝒞op × 𝒞 → 𝒞 
----------------------------------------------------------=#

struct MonadModules <: Category
    category::Category 
    T::Monad
end

struct MonadModule <: Object
    parent::MonadModules
    object::Object 
    action::Morphism 
end

struct MonadModuleMorphism <: Morphism 
    domain::MonadModule 
    codomain::MonadModule 
    m::Morphism 
end

MonadModule(X::Object, action::Morphism) = MonadModule(parent(X),X,action)

morphism(X::MonadModule, Y::MonadModule, m) = MonadModuleMorphism(X,Y,m)

monad(C::MonadModules) = C.monad

action(X::MonadModule) = X.action

#=----------------------------------------------------------
    Abelian 
----------------------------------------------------------=#

function direct_sum(X::MonadModule...)
    
    C = parent(X[1])

    S,i,p = direct_sum(object.(X)...)

    T = monad(C)

    action = horizontal_direct_sum([action(Xi) ∘ T(pi) for (Xi,pi) ∈ zip(X, p)])

    S = MonadModule(C, S, action)

    i = [morphism(X[j], S, i[j]) for j ∈ eachindex(X)]
    p = [morphism(S, X[j], p[j]) for j ∈ eachindex(X)]

    return S,i,p
end

function direct_sum(f::MonadModuleMorphism...)
    dom = direct_sum(domain.(f))[1]
    cod = direct_sum(codomain.(f))[1]

    morphism(dom,cod, direct_sum(morphism.(f)))
end

function *(k, f::MonadModuleMorphism) 
    morphism(domain(f), codomain(f), k.morphism(f))
end

function +(f::MonadModuleMorphism...)
    morphism(domain(f), codomain(f), +(collect(morphism.(f))))
end

function kernel(f::MonadModuleMorphism) 
    K,k = kernel(morphism(f))
    C = parent(f)
    T = monad(C)
    action = compose(
        T(k),
        action(domain(f)),
        left_inverse(k)
    )
    
    K = MonadModule(C, K, action)
    k = morphism(K, domain(f), k)
    return K, k 
end

function cokernel(f::MonadModuleMorphism) 
    C,c = kernel(morphism(f))
    M = parent(f)
    T = monad(M)
    action = compose(
        T(right_inverse(c)),
        action(codomain(f)),
        c
    )
    
    C = MonadModule(M, C, action)
    c = morphism(codomain(f), C, c)
    return C,c
end

#=----------------------------------------------------------
    monoidal structure  
----------------------------------------------------------=#

function tensor_product(X::MonadModule...)
end




