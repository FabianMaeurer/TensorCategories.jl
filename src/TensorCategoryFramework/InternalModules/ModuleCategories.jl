#=----------------------------------------------------------
    Structures for left, right and bimodules over 
    algebra objects in multitensor categories
----------------------------------------------------------=#

abstract type ModuleCategory <: Category end
abstract type ModuleObject <: Object end

@attributes mutable struct LeftModuleObject <: ModuleObject
    parent::ModuleCategory
    object::Object
    left_action::Morphism

    function LeftModuleObject(C::ModuleCategory, X::Object, l::Morphism)
        M = new()
        M.parent = C
        M.object = X
        M.left_action = l
        return M
    end
end

@attributes mutable struct RightModuleObject <: ModuleObject
    parent::ModuleCategory
    object::Object
    right_action::Morphism 

    function RightModuleObject(C::ModuleCategory, X::Object, r::Morphism)
        M = new()
        M.parent = C
        M.object = X
        M.right_action = r
        return M
    end
end

@attributes mutable struct BiModuleObject <: ModuleObject
    parent::ModuleCategory
    object::Object
    left_action::Morphism
    right_action::Morphism

    function BiModuleObject(C::ModuleCategory, X::Object, l::Morphism, r::Morphism)
        M = new()
        M.parent = C
        M.object = X
        M.left_action = l
        M.right_action = r
        return M
    end
end

@attributes mutable struct LeftModuleCategory <: ModuleCategory
    category::Category
    left_algebra::AlgebraObject
    simples::Vector{LeftModuleObject}

    function LeftModuleCategory(C::Category, A::AlgebraObject)
        M = new()
        M.category = C
        M.left_algebra = A
        return M
    end
end

@attributes mutable struct RightModuleCategory <: ModuleCategory
    category::Category
    right_algebra::AlgebraObject
    simples::Vector{RightModuleObject}

    function RightModuleCategory(C::Category, A::AlgebraObject)
        M = new()
        M.category = C
        M.right_algebra = A
        return M
    end
end

@attributes mutable struct BiModuleCategory <: ModuleCategory
    category::Category
    left_algebra::AlgebraObject
    right_algebra::AlgebraObject
    simples::Vector{BiModuleObject}

    function BiModuleCategory(C::Category, A::AlgebraObject, B::AlgebraObject)
        M = new()
        M.category = C
        M.left_algebra = A
        M.right_algebra = B
        return M
    end
end



struct ModuleMorphism{T<:ModuleObject} <: Morphism
    domain::T
    codomain::T
    map::Morphism
end

function morphism(X::T, Y::T, m::Morphism) where T <: ModuleObject
    ModuleMorphism(X,Y,m)
end


#=----------------------------------------------------------
    Getter/Setter 
----------------------------------------------------------=#

@doc raw""" 

    left_action(M::ModuleObject)

Return the left action morphism ``A⊗M → M``
"""
function left_action(M::ModuleObject) 
    !isdefined(M, :left_action) && error("Not a left module")

    M.left_action
end

@doc raw""" 

    right_action(M::ModuleObject)

Return the right action morphism ``M⊗A → M``
"""
function right_action(M::ModuleObject)
    !isdefined(M, :right_action) && error("Not a right module")

    M.right_action
end


base_ring(C::ModuleCategory) = base_ring(category(C))

left_algebra(C::ModuleCategory) = C.left_algebra
right_algebra(C::ModuleCategory) = C.right_algebra

algebra(C::LeftModuleCategory) = left_algebra(C)
algebra(C::RightModuleCategory) = right_algebra(C)

morphism(f::ModuleMorphism) = f.map

matrix(f::ModuleMorphism) = matrix(morphism(f))

is_abelian(C::ModuleCategory) = is_abelian(category(C))

is_multiring(C::BiModuleCategory) = true
is_ring(C::BiModuleCategory) = is_multiring(C) &&  int_dim(End(one(C))) == 1

is_multiring(C::RightModuleCategory) = is_commutative(algebra(C))
is_multiring(C::LeftModuleCategory) = is_commutative(algebra(C))

is_tensor(C::BiModuleCategory) = int_dim(End(one(C))) == 1

is_weak_multifusion(C::BiModuleCategory) = is_weak_multifusion(category(C)) &&
                                    left_algebra(C) == right_algebra(C) && 
                                    is_separable(left_algebra(C))

is_weak_fusion(C::BiModuleCategory) = is_weak_multifusion(C) && 
                                    int_dim(End(one(C))) == 1

is_multifusion(C::BiModuleCategory) = is_weak_multifusion(C) &&
                                    all(int_dim(End(s)) == 1 for s ∈ simples(C))

is_fusion(C::BiModuleCategory) = is_multifusion(C) && int_dim(End(one(C))) == 1

function ==(M::ModuleCategory, N::ModuleCategory)
    typeof(M) != typeof(N) && return false

    category(M) != category(N) && return false

    isdefined(M, :left_algebra) && left_algebra(M) != left_algebra(N) && return false

    isdefined(M, :right_algebra) && right_algebra(M) != right_algebra(N) && return false

    true
end

function ==(X::ModuleObject, Y::ModuleObject)
    typeof(X) != typeof(Y) && return false

    parent(X) != parent(Y) && return false

    isdefined(X, :left_action) && left_action(X) != left_action(Y) && return false

    isdefined(X, :right_action) && right_action(X) != right_action(Y) && return false

    true
end

morphism_type(C::LeftModuleCategory) = ModuleMorphism{LeftModuleObject}
morphism_type(C::RightModuleCategory) = ModuleMorphism{RightModuleObject}
morphism_type(C::BiModuleCategory) = ModuleMorphism{BiModuleObject}

#=----------------------------------------------------------
    Constructors 
----------------------------------------------------------=#

@doc raw""" 

    left_module(A::AlgebraObject)

Return ``A`` as the trivial left module
"""
function left_module(A::AlgebraObject)
    C = LeftModuleCategory(parent(A), A)
    LeftModuleObject(C, object(A), multiplication(A))
end

@doc raw""" 

    right_module(A::AlgebraObject)

Return ``A`` as the trivial right module
"""
function right_module(A::AlgebraObject)
    C = RightModuleCategory(parent(A), A)
    RightModuleObject(C, object(A), multiplication(A))
end

@doc raw""" 

    bimodule(A::AlgebraObject)

Return ``A`` as the trivial bimodule
"""
function bimodule(A::AlgebraObject)
    C = BiModuleCategory(parent(A), A,A)
    BiModuleObject(C, object(A), multiplication(A), multiplication(A))
end

@doc raw""" 

    right_module(M::BiModuleObject)

Return ``M`` as a right module forgetting the left module structure
"""
function right_module(M::BiModuleObject)
    C = RightModuleCategory(category(parent(M)), right_algebra(parent(M)))
    RightModuleObject(C, object(M), right_action(M))
end

@doc raw""" 

    left_module(M::BiModuleObject)

Return ``M`` as a left module forgetting the right module structure
"""
function left_module(M::BiModuleObject)
    C = LeftModuleCategory(category(parent(M)), left_algebra(parent(M)))
    LeftModuleObject(C, object(M), left_action(M))
end

function bimodule(X::RightModuleObject)
    M = parent(X)
    C = category(M)
    !is_braided(C) && error("No canonnical left module structure")

    A = right_algebra(M)

    !is_commutative(A) && error("Algebra is not commutative")

    l = right_action(X) ∘ braiding(object(A), object(X))

    M2 = BiModuleCategory(C,A,A)

    BiModuleObject(M2, object(X), l, right_action(X))
end


function left_module(M::RightModuleObject)
    A = algebra(parent(M))
    @assert is_commutative(A)

    LeftModuleObject(
        category_of_left_modules(A), 
        object(M), 
        right_action(M) ∘ braiding(object(A), object(M))
    )
end
#=----------------------------------------------------------
    Functionality 
----------------------------------------------------------=#

id(X::ModuleObject) = ModuleMorphism(X,X, id(object(X)))

function kernel(f::ModuleMorphism)
    X = domain(f)
    typeof(X) == LeftModuleObject  && return left_module_kernel(f)
    typeof(X) == RightModuleObject && return right_module_kernel(f)
    typeof(X) == BiModuleObject    && return bimodule_kernel(f)
end

function left_module_kernel(f::ModuleMorphism)
    X = domain(f)
    A = object(algebra(parent(X)))
    l = left_action(X)
    K,k = kernel(morphism(f))
    inv_k = left_inverse(k)
    ker_action = inv_k ∘ l ∘ (id(A)⊗k)
    ker = LeftModuleObject(parent(f), K, ker_action)
    return ker, ModuleMorphism(ker, X, k)
end


function right_module_kernel(f::ModuleMorphism)
    X = domain(f)
    A = object(algebra(parent(X)))
    r = right_action(X)
    K,k = kernel(morphism(f))
    inv_k = left_inverse(k)
    ker_action = inv_k ∘ r ∘ (k ⊗ id(A))
    ker = RightModuleObject(parent(f), K, ker_action)
    return ker, ModuleMorphism(ker, X, k)
end

function bimodule_kernel(f::ModuleMorphism)
    X = domain(f)
    A = object(left_algebra(parent(X)))
    B = object(right_algebra(parent(X)))
   
    K,k = kernel(morphism(f))
    inv_k = left_inverse(k)

    r = right_action(X)
    right_ker_action = inv_k ∘ r ∘ (k ⊗ id(B))

    l = left_action(X)
    left_ker_action = inv_k ∘ l ∘ (id(A)⊗k)

    ker = BiModuleObject(parent(f), K, left_ker_action, right_ker_action)

    return ker, ModuleMorphism(ker, X, k)
end

function cokernel(f::ModuleMorphism)
    X = domain(f)
    typeof(X) == LeftModuleObject  && return left_module_cokernel(f)
    typeof(X) == RightModuleObject && return right_module_cokernel(f)
    typeof(X) == BiModuleObject    && return bimodule_cokernel(f)
end

function left_module_cokernel(f::ModuleMorphism)
    Y = codomain(f)
    A = object(algebra(parent(Y)))
    l = left_action(Y)
    C,c = cokernel(morphism(f))
    inv_c = right_inverse(c)
    coker_action = c ∘ l ∘ (id(A) ⊗ inv_c)
    coker = LeftModuleObject(parent(f), C, coker_action)
    return coker, ModuleMorphism(Y, coker, c)
end


function right_module_cokernel(f::ModuleMorphism)
    Y = codomain(f)
    A = object(algebra(parent(Y)))
    r = right_action(Y)
    C,c = cokernel(morphism(f))
    inv_c = right_inverse(c)
    coker_action = c ∘ r ∘ (inv_c ⊗ id(A))
    coker = RightModuleObject(parent(f), C, coker_action)
    return coker, ModuleMorphism(Y, coker, c)
end

function bimodule_cokernel(f::ModuleMorphism)
    X = domain(f)
    A = object(left_algebra(parent(X)))
    B = object(right_algebra(parent(X)))

    C,c = cokernel(morphism(f))
    inv_c = right_inverse(c)

    l = left_action(X)
    left_coker_action = c ∘ l ∘ (id(A) ⊗ inv_c)

    r = right_action(X)
    right_coker_action = c ∘ r ∘ (inv_c ⊗ id(B))

    coker = BiModuleObject(parent(f), C, left_coker_action, right_coker_action)

    return coker, ModuleMorphism(X, coker, c)
end

function inv(M::ModuleMorphism)
    morphism(codomain(M), domain(M), inv(morphism(M)))
end
#=----------------------------------------------------------
    Free modules   
----------------------------------------------------------=#

@doc raw""" 

    free_left_module(X::Object, A::AlgebraObject)

Return the free left module ``A⊗X``
"""
function free_left_module(X::Object, A::AlgebraObject, parent = LeftModuleCategory(parent(X), A))
    AA = object(A)
    M = AA ⊗ X
    l = compose(
        inv_associator(AA, AA, X),
        multiplication(A) ⊗ id(X)
    )
    LeftModuleObject(parent, M, l)
end

@doc raw""" 

    free_right_module(X::Object, A::AlgebraObject)

Return the free right module ``X⊗A``
"""
function free_right_module(X::Object, A::AlgebraObject, parent = RightModuleCategory(parent(X), A))
    AA = object(A)
    M = X ⊗ AA
    l = compose(
        associator(X, AA, AA),
        id(X) ⊗ multiplication(A)
    )
    RightModuleObject(parent, M, l)
end

@doc raw""" 

    free_bimodule(X::Object, A::AlgebraObject)

Return the free ``A-A`` bimodule ``A⊗X⊗A``
"""
function free_bimodule(X::Object, A::AlgebraObject, parent = BiModuleCategory(parent(X), A, A))
    free_bimodule(X,A,A,parent)
end

function free_bimodule(X::RightModuleObject, A::AlgebraObject, parent_cat = BiModuleCategory(category(parent(X)), A, right_algebra(parent(X))))
    LX = free_left_module(object(X), A)
    B = right_algebra(parent(X))
    right = compose(
        associator(object(A), object(X), object(B)),
        id(object(A)) ⊗ right_action(X)
    )
    BiModuleObject(parent_cat, object(LX), left_action(LX), right)
end

@doc raw""" 

    free_bimodule(X::Object, A::AlgebraObject, B::AlgebraObject)

Return the free ``A-B`` bimodule ``A⊗X⊗B``
"""
function free_bimodule(X::Object, A::AlgebraObject, B::AlgebraObject, parent = BiModuleCategory(parent(X), A, B))
    AA = object(A)
    BB = object(B)
    AX = AA ⊗ X
    M = AX ⊗ BB

    # left action on A⊗X
    ll = compose(
        inv_associator(AA, AA, X),
        multiplication(A) ⊗ id(X)
    )

    # left action on (A⊗X)⊗B
    l = compose(
        inv_associator(AA, AX, BB),
        ll ⊗ id(BB),
    )

    # right action on (A ⊗ X) ⊗ b
    r = compose(
        associator(AX, BB, BB),
        id(AX) ⊗ multiplication(B)
    ) 

    BiModuleObject(parent, M, l, r)
end

@doc raw""" 

    free_module(X::Object, M::ModuleCategory)

Return the free module of ``X`` in ``M``
"""
function free_module(X::Object, M::ModuleCategory)
    if typeof(M) == LeftModuleCategory
        return free_left_module(X,left_algebra(M), M)
    elseif typeof(M) == RightModuleCategory
        return free_right_module(X, right_algebra(M), M)
    else
        return free_bimodule(X, left_algebra(M), right_algebra(M), M)
    end
end


@doc raw""" 

    transposed_module(M::RightModuleObject)

Construct the left module ``(∗M, q)`` from ``(M,p)``. 
"""
function transposed_module(N::RightModuleObject)
    p = right_action(N)
    M = object(N)
    dM = dual(M)
    A  = object(algebra(parent(N)))
    dA = dual(A)


    q1 = compose(
        coev(M),
        compose(
            (id(M) ⊗ coev(A)),
            inv_associator(M,A,dA),
            p ⊗ id(dA)
         ) ⊗ id(dM),
            associator(M,dA,dM)
    )
    q2 = compose(
        id(dM) ⊗ q1,
        inv_associator(dM, M, (dA ⊗ dM)),
        ev(M) ⊗ (id(dA) ⊗ id(dM))
    )

    q = compose(
        id(A) ⊗ q2,
        inv_associator(A, dA, dM),
        (ev(dA) ∘ (spherical(A) ⊗ id(dA))) ⊗ id(dM)
    )

    C = LeftModuleCategory(
        parent(M),
        algebra(parent(N))
    )

    LeftModuleObject(C, dM, q)
end

function transposed_module(N::LeftModuleObject)
    q = left_action(N)
    M = object(N)
    dM = dual(M)
    A = object(algebra(parent(N)))
    dA = dual(A)

    p1 = compose(
        coev(dM),
        id(dM) ⊗ inv(spherical(M)),
        id(dM) ⊗ compose(
            coev(dA) ⊗ id(M) ,
            (id(dA) ⊗ inv(spherical(A))) ⊗ id(M),
            associator(dA,A,M),
            id(dA) ⊗ q
        ),
        inv_associator(dM, dA, M)
    )

    p2 = compose(
        p1 ⊗ id(dM),
        associator(dM ⊗ dA, M, dM),
        (id(dM) ⊗ id(dA)) ⊗ (ev(dM) ∘ (spherical(M) ⊗ id(dM)))
    )

    p = compose(
        p2 ⊗ id(A),
        associator(dM, dA, A),
        id(dM) ⊗ ev(A)
    ) 

    C = RightModuleCategory(
        parent(N),
        algebra(parent(N))
    )

    RightModuleObject(C, dM, p)
end



#=----------------------------------------------------------
    Monoidal Structe of Bimod
----------------------------------------------------------=#

function tensor_product(M::RightModuleObject, N::LeftModuleObject)
    _tensor_product(M,N)[1]
end

function _tensor_product(M::RightModuleObject, N::LeftModuleObject)
    @assert algebra(parent(M)) == algebra(parent(N))

    A = algebra(parent(M))
    p = right_action(M)
    q = left_action(N)

    C,c = coequilizer(
        p ⊗ id(object(N)),
        (id(object(M)) ⊗ q) ∘ associator(object(M), object(A), object(N))
    )

    return C,c
end

function tensor_product(M::RightModuleObject, N::RightModuleObject)
    right_module_tensor_product(M,N)[1]
end

@memoize Dict function right_module_tensor_product(M::RightModuleObject, N::RightModuleObject)
    A = algebra(parent(M))
    @assert is_commutative(A)
    
    one_C = free_right_module(one(category(parent(M))), A)

    if N == one_C
        return M, right_action(M)
    end 

    N2 = left_module(N)

    if M == one_C
        return N, left_action(N2)
    end

    MN,c = _tensor_product(M, N2)
    inv_c = right_inverse(c)
    RightModuleObject(
        parent(M),
        MN,
        compose(
            inv_c ⊗ id(object(A)),
            associator(object(M), object(N), object(A)),
            id(object(M)) ⊗ right_action(N),
            c 
        )
    ), c
end
         

function tensor_product(M::BiModuleObject, N::BiModuleObject)
    bimodule_tensor_product(M,N)[1]
end

function bimodule_tensor_product(M::BiModuleObject, N::BiModuleObject)
    A = left_algebra(parent(M))
    B = right_algebra(parent(M))
    C = left_algebra(parent(N))
    D = right_algebra(parent(N))

    @assert B == C "Mismatching algebras"

    if A == B == C == D
        C = parent(M)

        # Make the unitors strict
        if M == bimodule(A) 
            return N, left_action(N)
        elseif N == bimodule(A)
            return M, right_action(M)
        end
    else
        C = BiModuleCategory(category(parent(M)), A, D)
    end

    p = right_action(M)
    q = left_action(N)

    P,c = coequilizer(
        p ⊗ id(object(N)),
        (id(object(M)) ⊗ q) ∘ associator(object(M), object(A), object(N))
    )

    inv_c = right_inverse(c)

    left = compose(
        id(object(A)) ⊗ inv_c,
        inv_associator(object(A), object(M), object(N)),
        left_action(M) ⊗ id(object(N)),
        c
    )

    right = compose(
        inv_c ⊗ id(object(D)),
        associator(object(M), object(N), object(B)),
        id(object(M)) ⊗ right_action(N),
        c
    )

    BiModuleObject(C, P, left, right), c
end

# function tensor_product(X::RightModuleObject, Y::RightModuleObject)
#     right_module_tensor_product(X,Y)[1]
# end

# function right_module_tensor_product(X::RightModuleObject, Y::RightModuleObject)
#     !(parent(X) == parent(Y)) && error("Mismatching parents")

#     M = parent(X) 

#     !is_braided(category(M)) && error("Category is not braided")

#     A = object(right_algebra(M))

#     ass = associator(object(X), object(Y), A)

#     top = compose(
#         ass,
#         id(object(X)) ⊗ right_action(Y)
#     )

#     bottom = compose(
#         ass,
#         id(object(X)) ⊗ braiding(object(Y), A),
#         inv_associator(object(X), A, object(Y)),
#         right_action(X) ⊗ id(object(Y))
#     )

#     Z, c = coequilizer(
#         top, bottom
#     )

#     inv_c = right_inverse(c)

#     right = compose(
#         inv_c ⊗ id(A),
#         top,
#         c
#     )
    
#     RightModuleObject(M, Z, right), c
# end

# function right_module_tensor_product(X::RightModuleObject, Y::RightModuleObject)
#     S = object(X) ⊗ Y
#     S, morphism(id(S))
# end

function tensor_product(f::ModuleMorphism{T}, g::ModuleMorphism{T}) where T <: ModuleObject
    T == BiModuleObject && (mod_tensor = bimodule_tensor_product)
    T == RightModuleObject && (mod_tensor = right_module_tensor_product)
    T == LeftModuleObject && (mod_tensor = left_module_tensor_product)

    dom, p_dom = mod_tensor(domain(f), domain(g))
    cod, p_cod = mod_tensor(codomain(f), codomain(g))

    inv_p_dom = right_inverse(p_dom)

    ModuleMorphism(dom, cod, p_cod ∘ (morphism(f) ⊗ morphism(g)) ∘ inv_p_dom)
end


function associator(X::T, Y::T, Z::T) where T <: ModuleObject

    if one(parent(X)) ∈ [X,Y,Z]
        return id(X ⊗ Y ⊗ Z)
    end
    
    T == BiModuleObject && (mod_tensor = bimodule_tensor_product)
    T == RightModuleObject && (mod_tensor = right_module_tensor_product)
    T == LeftModuleObject && (mod_tensor = left_module_tensor_product)

    XY, p_XY = mod_tensor(X,Y)
    XY_Z, p_XY_Z = mod_tensor(XY, Z)

    YZ, p_YZ = mod_tensor(Y, Z)
    X_YZ, p_X_YZ = mod_tensor(X, YZ)

    before = compose(
        right_inverse(p_XY_Z),
        right_inverse(p_XY) ⊗ id(object(Z))
        )
    after = compose(
        id(object(X)) ⊗ p_YZ,
        p_X_YZ
    )

    a = compose(
        before,
        associator(object(X), object(Y), object(Z)),
        after
    )
    morphism(XY_Z, X_YZ, a)
end

function one(C::BiModuleCategory)
    @assert left_algebra(C) == right_algebra(C)
    A = left_algebra(C)
    bimodule(A)
end

function one(C::Union{RightModuleCategory, LeftModuleCategory})
    @assert is_commutative(algebra(C))
    return free_module(one(category(C)), C)
end

function dual(M::BiModuleObject)
    #@assert left_algebra(parent(M)) == right_algebra(parent(M))

    lM = transposed_module(right_module(M))
    rM = transposed_module(left_module(M))
    A = left_algebra(parent(M))
    B = right_algebra(parent(M))

    if A == B
        C = parent(M)
    else
        C = BiModuleCategory(
            category(parent(M)),
            B,
            A
        )
    end

    BiModuleObject(
        C,
        object(lM),
        left_action(lM),
        right_action(rM)
    )
end

# function ev(M::BiModuleObject)
#     dM = dual(M)
#     dMM,p = bimodule_tensor_product(dM,M)
#     inv_p = right_inverse(p)
#     u = unit(right_algebra(parent(M)))
#     morphism(dMM, one(parent(M)), u ∘ ev(object(M)) ∘ inv_p)
# end

# function coev(M::BiModuleObject)
#     dM = dual(M)
#     MdM, p = bimodule_tensor_product(M,dM)

#     base = basis(Hom(one(parent(M)), MdM))
#     M_basis = basis(End(M))

#     m = (id(M)⊗ev(M))∘associator(M,dM,M)

#     M = hcat([express_in_basis(m∘(f ⊗ id(M)), M_basis) for f ∈ base]...)
#     one_coeffs = express_in_basis(id(one(parent(dM))), M_basis)

#     s = solve(transpose(matrix(base_ring(dM), size(M,1), size(M,2), M)), one_coeffs)

#     sum(s .* base)
# end

function spherical(M::BiModuleObject)
    ddM = dual(dual(M))
    morphism(M, ddM, spherical(object(M)))
end

function tr(f::ModuleMorphism{BiModuleObject}) 
    C = parent(f)
    A,B = left_algebra(C), right_algebra(C)

    d = sqrt(dim(object(A))*dim(object(B)))

    t = base_ring(f)(tr(morphism(f)))*inv(d) 
    t * id(one(parent(f)))
end

#=----------------------------------------------------------
    Direct sum 
----------------------------------------------------------=#

function direct_sum(M::T...) where T <: ModuleObject
    if length(M) == 1
        return M[1], [id(M[1])], [id(M[1])]
    end

    Z, incl, proj = direct_sum(object.(M)...)
    C = parent(M[1])
    
    if isdefined(M[1], :right_action)
        B = object(right_algebra(parent(M[1])))
        r = compose(
            distribute_left([object(m) for m ∈ M], B),
            direct_sum(right_action.(M)...)
        )
    end
    
    if isdefined(M[1], :left_action)
        A = object(left_algebra(parent(M[1])))
        l = compose(
            distribute_right(A, [object(m) for m ∈ M]),
            direct_sum(left_action.(M)...)
        )
    end

    if T == RightModuleObject 
        Z =  RightModuleObject(C,Z,r)
        incl = [morphism(m,Z,i) for (m,i) ∈ zip(M,incl)]
        proj = [morphism(m,Z,p) for (m,p) ∈ zip(M,proj)]
        return Z, incl, proj
    elseif T == LeftModuleObject
        Z = LeftModuleObject(C,Z,l)
        incl = [morphism(m,Z,i) for (m,i) ∈ zip(M,incl)]
        proj = [morphism(m,Z,p) for (m,p) ∈ zip(M,proj)]
        return Z, incl, proj
    else
        Z = BiModuleObject(C,Z,l,r)
        incl = [morphism(m,Z,i) for (m,i) ∈ zip(M,incl)]
        proj = [morphism(m,Z,p) for (m,p) ∈ zip(M,proj)]
        return Z, incl, proj
    end
end

function direct_sum(f::ModuleMorphism...)
    dom = direct_sum(domain.(f)...)[1]
    cod = direct_sum(codomain.(f)...)[1]
    morphism(dom, cod, direct_sum(morphism.(f)...))
end

function decompose(M::Union{LeftModuleObject,RightModuleObject})
    A = algebra(parent(M))
    # if is_separable(A)
    #     return minimal_subquotients_with_multiplicity(M)
    # else
        return decompose_by_endomorphism_ring(M)
   # end
end
#=----------------------------------------------------------
    Internal Hom 
----------------------------------------------------------=#

function internal_hom(M::RightModuleObject, N::RightModuleObject)
    dual(M ⊗ transposed_module(N))
end

function internal_hom_adjunction(X::Object, M::RightModuleObject, N::RightModuleObject, f::Morphism)

    A = algebra(parent(M))
    p = right_action(M)
    dN = transposed_module(N)
    q = left_action(dN)
    
    # Build the internal hom with the projection 
    # M ⊗ ∗N → M ⊗ₐ ∗N 
    H,h = _tensor_product(M, transposed_module(N))

    # Internal Hom is given by (M ⊗ₐ ∗N)∗
    H = dual(H)

    # unpack (M ⊗ₐ ∗N)∗ = N ⊗ M∗
    dual_monoidal = dual_monoidal_structure(object(M), object(dN))

    f = dual_monoidal ∘ dual(h) ∘ f

    g = compose(
        f ⊗ id(object(M)),
        associator(object(N), dual(object(M)), object(M)),
        id(object(N)) ⊗ ev(object(M))
    )

    return ModuleMorphism(X ⊗ M, N, g)
end

function inverse_internal_hom_adjunction(X::Object, M::RightModuleObject, N::RightModuleObject, f::ModuleMorphism)

    A = algebra(parent(M))
    p = right_action(M)
    dN = transposed_module(N)
    q = left_action(dN)
    
    # Build the internal hom manually to obtain the projection 
    # M ⊗ ∗N → M ⊗ₐ ∗N 
  
    H,h = coequilizer(
        p ⊗ id(object(dN)),
        (id(object(M)) ⊗ q) ∘ associator(object(M), object(A), object(dN))
    ) 

    incl = right_inverse(h)

    # Internal Hom is given by (M ⊗ₐ ∗N)∗
    H = dual(H)

    # Hom(XM, N) → Hom((∗M∗X)∗, N)
    dual_monoidal = dual_monoidal_structure(X, dual(object(M)))
    g = morphism(f) ∘ dual_monoidal

    # Hom((∗M∗X)∗, N) → Hom(∗N, ∗M∗X)
    g = dual(g)

    # Hom(∗N, ∗M∗X) → Hom(M∗N, ∗X)
    g = compose(
        id(object(M)) ⊗ g,
        inv_associator(object(M), dual(object(M)), dual(X)),
        (ev(dual(object(M))) ∘ (spherical(object(M)) ⊗ id(dual(object(M))))) ⊗ id(dual(X))
    )

    # Hom(M∗N, ∗X) → Hom(X,(M∗N)∗)
    g = dual(g ∘ incl)

    return g
end

function ev(M::RightModuleObject, N::RightModuleObject)
    H = internal_hom(M,N)
    internal_hom_adjunction(H, M, N, id(H))
end

function unit(M::ModuleObject)
    H = internal_hom(M,M)
    internal_hom_adjunction(one(category(parent(M))), M, M, id(H))
end

function zero(C::RightModuleCategory)
    RightModuleObject(C,zero(category(C)), zero_morphism(category(C)))
end

function zero(C::LeftModuleCategory)
    LeftModuleObject(C,zero(category(C)), zero_morphism(category(C)))
end

function zero(C::BiModuleCategory)
    BiModuleObject(C,zero(category(C)), zero_morphism(category(C)), zero_morphism(category(C)))
end


#=----------------------------------------------------------
    Module Category Structure 
----------------------------------------------------------=#

function left_module_action(X::Object, M::RightModuleObject) 
    p = compose( 
        associator(X, object(M), object(algebra(parent(M)))),
        id(X) ⊗ right_action(M)
    )
    RightModuleObject(parent(M), X ⊗ object(M), p)
end

function right_module_action(X::Object, M::LeftModuleObject)
    p = compose( 
        inv_associator(object(algebra(parent(M))), object(M), X),
        left_action(M) ⊗ id(X)
    )
    LeftModuleObject(parent(M), object(M) ⊗ X, p)
end

function module_action(X::Object, M::T)  where T <: Union{LeftModuleObject, RightModuleObject}
    if typeof(X) == LeftModuleObject
        return right_module_action(X, M)
    else
        return left_module_action(X, M)
    end
end

function ⊗(X::Object, M::ModuleObject) 
    if parent(X) == parent(M) 
        return tensor_product(X,M) 
    end
    left_module_action(X, M)
end

⊗(M::ModuleObject, N::ModuleObject) = tensor_product(M,N)
⊗(M::ModuleObject, X::Object) = right_module_action(X, M)


#=----------------------------------------------------------
    Equivalence of Module categories 
----------------------------------------------------------=#

@doc raw""" 

    is_equivalent(M::ModuleCategory, N::ModuleCategory)     

Check if ``M`` and ``N`` are equivalent as module categories.
The equivalence is provided by a suitable bimodule object.
"""
function is_equivalent(M::T, N::T) where T <: Union{LeftModuleCategory, RightModuleCategory}
    Func = category_of_bimodules(algebra(M), algebra(N))

    S = simples(Func)

    for s ∈ S
        
        # If End(s) > 1 s cannot be an isomorphism
        int_dim(End(s)) > 1 && continue
        
        dual_s = dual(s)
        sds = s ⊗ dual(s)
        dss = dual(s) ⊗ s
        NN = parent(sds)
        MM = parent(dss)
        
        # If id(M) ≠ s∗ ∘ s continue
        !is_isomorphic(dss, one(NN))[1] && continue
        # If also s ∘ s∗ = id(N) we found an isomorphism
        is_isomorphic(sds, one(MM))[1] && return true, s
    end

    return false, nothing
end
    

#=----------------------------------------------------------
    Checks 
----------------------------------------------------------=#

function is_bimodule_morphism(f::Morphism, M::BiModuleObject, N::BiModuleObject)
    @assert domain(f) == object(M)
    @assert codomain(f) == object(N)

    A = object(left_algebra(parent(M)))
    B = object(right_algebra(parent(M)))

    left_action(N) ∘ (id(A) ⊗ f) == f ∘ left_action(M) &&
    right_action(N) ∘ (f ⊗ id(A)) == f ∘ right_action(M)
end

function is_bimodule_morphism(f::ModuleMorphism{BiModuleObject})
    is_bimodule_morphism(morphism(f), domain(f), codomain(f))
end

function is_simple(M::ModuleObject)
    has_attribute(M, :is_simple) && return get_attribute(M, :is_simple)
    bool = length(minimal_subquotients(M)) == 1
    set_attribute!(M, :is_simple, bool)
    return bool
end

#=----------------------------------------------------------
    Hom Spaces
----------------------------------------------------------=#

function Hom(X::LeftModuleObject, Y::LeftModuleObject)
    H = Hom(object(X), object(Y))
    B = basis(H)
    F = base_ring(X)
    n = length(basis(H))

    A = algebra(parent(X))

    if n == 0 
        return HomSpace(X,Y, ModuleMorphism[])
    end 

    Fx,poly_basis = polynomial_ring(F,n)

    base = basis(Hom(object(A)⊗object(X), object(Y)))

    eqs = [zero(Fx) for _ ∈ 1:length(base)]

    l_X = left_action(X)
    l_Y = left_action(Y)

    for (f,a) ∈ zip(B,poly_basis)
        coeffs = express_in_basis((f ∘ l_X) - (l_Y ∘ (id(object(A))⊗f)), base)
        eqs = eqs .+ (a .* coeffs) 
    end
    

    M = zero(matrix_space(F,length(eqs),n))

    for (i,e) ∈ zip(1:length(eqs),eqs)
        M[i,:] = [coeff(e, a) for a ∈ poly_basis]
    end

    N = nullspace(M)[2]

    _,cols = size(N)

    basis_coeffs = [N[:,i] for i ∈ 1:cols]

    T = ModuleMorphism{LeftModuleObject}
    module_hom_basis = T[morphism(X,Y,sum(b .* B)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,module_hom_basis)
end

function zero_morphism(M::ModuleObject, N::ModuleObject)
    morphism(M,N, zero_morphism(object(M), object(N)))
end

function Hom(X::RightModuleObject, Y::RightModuleObject)
    H = Hom(object(X), object(Y))
    B = basis(H)
    F = base_ring(X)
    n = length(basis(H))

    A = algebra(parent(X))

    if n == 0 
        return HomSpace(X,Y, ModuleMorphism[])
    end 

    Fx,poly_basis = polynomial_ring(F,n)

    base = basis(Hom(object(X)⊗object(A), object(Y)))

    eqs = [zero(Fx) for _ ∈ 1:length(base)]

    r_X = right_action(X)
    r_Y = right_action(Y)

    for (f,a) ∈ zip(B,poly_basis)
        coeffs = express_in_basis((f ∘ r_X) - (r_Y ∘ (f⊗id(object(A)))), base)
        eqs = eqs .+ (a .* coeffs) 
    end
    

    M = zero(matrix_space(F,length(eqs),n))

    for (i,e) ∈ zip(1:length(eqs),eqs)
        M[i,:] = [coeff(e, a) for a ∈ poly_basis]
    end

    N = nullspace(M)[2]

    _,cols = size(N)

    basis_coeffs = [N[:,i] for i ∈ 1:cols]

    T = ModuleMorphism{RightModuleObject}
    module_hom_basis = T[morphism(X,Y,sum(b .* B)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,module_hom_basis)
end

function Hom(X::BiModuleObject, Y::BiModuleObject)
    @assert parent(X) == parent(Y)
    H = Hom(right_module(X), right_module(Y))
    B_H = morphism.(basis(H))
    F = base_ring(X)
    n = length(basis(H))

    A = left_algebra(parent(X))
    B = right_algebra(parent(X))

    if n == 0 
        return HomSpace(X,Y, ModuleMorphism[])
    end 

    Fx,poly_basis = polynomial_ring(F,n)

    base_r = basis(Hom(object(X)⊗object(B), object(Y)))
    base_l = basis(Hom(object(A)⊗object(X), object(Y)))

    eqs_r = [zero(Fx) for _ ∈ 1:length(base_r)]
    eqs_l = [zero(Fx) for _ ∈ 1:length(base_l)]

    r_X = right_action(X)
    r_Y = right_action(Y)
    l_X = left_action(X)
    l_Y = left_action(Y)

    for (f,a) ∈ zip(B_H,poly_basis)
        # Condition for right module morphism
        coeffs = express_in_basis((f ∘ r_X) - (r_Y ∘ (f⊗id(object(B)))), base_r)
        eqs_r = eqs_r .+ (a .* coeffs) 

        # Condition for left module morphism
        coeffs = express_in_basis((f ∘ l_X) - (l_Y ∘ (id(object(A))⊗f)), base_l)
        eqs_l = eqs_l .+ (a .* coeffs) 
    end
    
    eqs = [eqs_l; eqs_r]

    M = zero(matrix_space(F,length(eqs),n))

    for (i,e) ∈ zip(1:length(eqs),eqs)
        M[i,:] = [coeff(e, a) for a ∈ poly_basis]
    end

    N = nullspace(M)[2]

    _,cols = size(N)

    basis_coeffs = [N[:,i] for i ∈ 1:cols]

    T = ModuleMorphism{BiModuleObject}
    module_hom_basis = T[morphism(X,Y,sum(b .* B_H)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,module_hom_basis)
end


function free_adjunction(X::Object, M::RightModuleObject, hom::Vector{<:Morphism} = basis(Hom(X, object(M))))

    A = algebra(parent(M))
    base = [right_action(M) ∘ (f ⊗ id(object(A))) for f ∈ hom]
    IX = free_right_module(X,A, parent(M))
    base = [morphism(IX,M, f) for f ∈ base]
    HomSpace(IX , M, base)
end

function free_adjunction(X::RightModuleObject, M::BiModuleObject, hom::Vector{<:Morphism} = basis(Hom(X, right_module(M))))

    A = left_algebra(parent(M))
    base = [left_action(M) ∘ (id(object(A)) ⊗ morphism(f)) for f ∈ hom]
    IX = free_bimodule(X,A, parent(M))
    base = [morphism(IX,M, f) for f ∈ base]
    HomSpace(IX , M, base)
end


function end_of_free_bimodule(X::Object, A::AlgebraObject, B::AlgebraObject)
    M = free_bimodule(X, A, B)
    H = Hom(X,object(M))

    after = compose(
        left_action(M) ⊗ id(object(B)),
        right_action(M)
    )
    base = [after ∘ ((id(object(A)) ⊗ f) ⊗ id(object(B))) for f ∈ H]
    base = [morphism(M,M, f) for f ∈ base]

    HomSpace(M,M, base)
end

function end_of_free_bimodule(X::Object, A::AlgebraObject)
    end_of_free_bimodule(X,A,A)
end


function (R::Ring)(f::ModuleMorphism)
    R(morphism(f))
end
#=----------------------------------------------------------
    isomorphic 
----------------------------------------------------------=#

function is_isomorphic(M::ModuleObject, N::ModuleObject)
    if !is_isomorphic(object(M), object(N))[1]
        return false, nothing
    end
    H = Hom(M,N)

    for f ∈ H
        if is_invertible(f)
            return true, f
        end
    end
    if is_simple(M) && is_simple(N)
        return int_dim(H) ≥ 1 ? (true, basis(H)[1]) : (false,nothing)
    elseif int_dim(End(M)) != int_dim(End(N))
        return false, nothing
    else
        error("not implemented yet")
    end
end


#=----------------------------------------------------------
    Compute Module Categories 
----------------------------------------------------------=#

is_semisimple(M::LeftModuleCategory) = is_separable(left_algebra(M))
is_semisimple(M::RightModuleCategory) = is_separable(right_algebra(M))
is_semisimple(M::BiModuleCategory) = is_separable(right_algebra(M)) && is_separable(left_algebra(M))

@doc raw""" 

    category_of_right_modules(A::AlgebraObject)

Return the category of right ``A`` modules in parent(A).
"""
function category_of_right_modules(A::AlgebraObject)
    M = RightModuleCategory(parent(object(A)), A)
end

@doc raw""" 

    category_of_left_modules(A::AlgebraObject)

Return the category of left ``A`` modules in parent(A)
"""
function category_of_left_modules(A::AlgebraObject)
    M = LeftModuleCategory(parent(object(A)), A)
end

@doc raw""" 

    category_of_bimodules(A::AlgebraObject, B::AlgebraObject)

Return the category of ``A-B`` bimodules in parent(A)
"""
function category_of_bimodules(A::AlgebraObject, B::AlgebraObject)
    M = BiModuleCategory(parent(object(A)), A,B)
end

@doc raw""" 

    category_of_bimodules(A::AlgebraObject)

Return the category of ``A-A`` bimodules in parent(A)
"""
function category_of_bimodules(A::AlgebraObject)
    category_of_bimodules(A,A)
end

function simples(M::ModuleCategory)
    if isdefined(M, :simples)
        return M.simples
    end

    C = category(M)

    if typeof(M) == BiModuleCategory && is_semisimple(M)
        simpls = bimodule_simples(M)
        M.simples = simpls
        return simpls
    end

    free_objects = [free_module(s,M) for s ∈ simples(C)]

    simpls = unique_simples(vcat([minimal_subquotients(x) for x ∈ free_objects]...))

    M.simples = simpls
end

function bimodule_simples(M::BiModuleCategory)
    A = left_algebra(M)
    B = right_algebra(M)

    right_modules = category_of_right_modules(B)
    simple_right_modules = simples(right_modules)

    free_bimodules = [free_bimodule(x, A, M) for x ∈ simple_right_modules]
    homs = [free_adjunction(x, Ix) for (x,Ix) ∈ zip(simple_right_modules, free_bimodules)]

    simpls = vcat([simple_subobjects(m, H) for (m,H) in zip(free_bimodules, homs)]...)

    if A == B
        return unique_simples([one(M); simpls])
    else
        return unique_simples(simpls)
    end
end
#=----------------------------------------------------------
    action matrices 
----------------------------------------------------------=#


function action_matrix(X::Object, M::ModuleCategory)
    @assert is_semisimple(M)
    S = simples(M)
    end_dims = [int_dim(End(m)) for m ∈ S]
    
    acts = [module_action(X,m) for m ∈ S]

    [div(int_dim(Hom(Y,m)), d) for Y ∈ acts, (m,d) ∈ zip(S,end_dims)]
end



#=----------------------------------------------------------
    Checks 
----------------------------------------------------------=#

function is_right_module(X::ModuleObject)
    A = right_algebra(parent(X))
    m = multiplication(A)
    a = associator(object(X), object(A), object(A))
    r_X = right_action(X)

    r_X ∘ (r_X ⊗ id(object(A))) == r_X ∘ (id(object(X)) ⊗ m) ∘ a
end

function is_left_module(X::ModuleObject)
    A = left_algebra(parent(X))
    m = multiplication(A)
    a = associator(object(A), object(A), object(X))
    l_X = left_action(X)

    l_X ∘ (m ⊗ id(object(X))) == l_X ∘ (id(object(A)) ⊗ l_X ) ∘ a
end

function is_bimodule(X::BiModuleObject)
    !is_left_module(X) && return false
    !is_right_module(X) && return false

    A = left_algebra(parent(X))
    B = right_algebra(parent(X))

    l_X = left_action(X)
    r_X = right_action(X)

    first = compose(
        associator(object(A), object(X), object(B)),
        id(object(A)) ⊗ r_X,
        l_X
    )

    second = compose(
        l_X ⊗ id(object(B)),
        r_X
    )

    first == second
end
#=----------------------------------------------------------
    Pretty printing
----------------------------------------------------------=#

function show(io::IO, C::ModuleCategory)
    typeof(C) == LeftModuleCategory && 
        print(io, """Category of left modules over $(algebra(C))""")
    typeof(C) == RightModuleCategory && 
        print(io, """Category of right modules over $(algebra(C))""")
    typeof(C) == BiModuleCategory && 
        print(io, """Category of bimodules over $(left_algebra(C))-$(right_algebra(C))""")
end

function show(io::IO, X::ModuleObject)
    typeof(X) == LeftModuleObject && 
        print(io, """Left module: $(object(X))""")
    typeof(X) == RightModuleObject && 
        print(io, """Right module: $(object(X))""")
    typeof(X) == BiModuleObject && 
        print(io, """Bimodule: $(object(X))""")
end


