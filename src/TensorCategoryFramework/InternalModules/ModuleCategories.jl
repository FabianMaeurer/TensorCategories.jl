#=----------------------------------------------------------
    Structures for left, right and bimodules over 
    algebra objects in multitensor categories
----------------------------------------------------------=#

abstract type ModuleCategory <: Category end
abstract type ModuleObject <: Object end

struct LeftModuleObject <: ModuleObject
    parent::ModuleCategory
    object::Object
    left_action::Morphism
end

struct RightModuleObject <: ModuleObject
    parent::ModuleCategory
    object::Object
    right_action::Morphism 
end

struct BiModuleObject <: ModuleObject
    parent::ModuleCategory
    object::Object
    left_action::Morphism
    right_action::Morphism
end

mutable struct LeftModuleCategory <: ModuleCategory
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

mutable struct RightModuleCategory <: ModuleCategory
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

mutable struct BiModuleCategory <: ModuleCategory
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

function Morphism(X::T, Y::T, m::Morphism) where T <: ModuleObject
    ModuleMorphism(X,Y,m)
end


#=----------------------------------------------------------
    Getter/Setter 
----------------------------------------------------------=#

function left_action(M::ModuleObject) 
    !isdefined(M, :left_action) && error("Not a left module")

    M.left_action
end


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

is_multitensor(C::BiModuleCategory) = true

is_tensor(C::BiModuleCategory) = int_dim(End(one(C))) == 1

is_weak_multifusion(C::BiModuleCategory) = left_algebra(C) == right_algebra(C) && is_seperable(left_algebra(C))

is_weak_fusion(C::BiModuleCategory) = is_weak_multifusion(C) && int_dim(End(one(C))) == 1
#=----------------------------------------------------------
    Constructors 
----------------------------------------------------------=#

function left_module(A::AlgebraObject)
    C = LeftModuleCategory(parent(A), A)
    LeftModuleObject(C, object(A), multiplication(A))
end

function right_module(A::AlgebraObject)
    C = RightModuleCategory(parent(A), A)
    RightModuleObject(C, object(A), multiplication(A))
end

function bimodule(A::AlgebraObject)
    C = BiModuleCategory(parent(A), A,A)
    BiModuleObject(C, object(A), multiplication(A), multiplication(A))
end

function right_module(M::BiModuleObject)
    C = RightModuleCategory(category(parent(M)), right_algebra(parent(M)))
    RightModuleObject(C, object(M), right_action(M))
end

function left_module(M::BiModuleObject)
    C = LeftModuleCategory(category(parent(M)), left_algebra(parent(M)))
    LeftModuleObject(C, object(M), left_action(M))
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
    X = domain(f)
    A = object(algebra(parent(X)))
    l = left_action(X)
    C,c = cokernel(morphism(f))
    inv_c = right_inverse(c)
    coker_action = c ∘ l ∘ (id(A) ⊗ inv_c)
    coker = LeftModuleObject(parent(f), C, coker_action)
    return coker, ModuleMorphism(X, coker, c)
end


function right_module_cokernel(f::ModuleMorphism)
    X = domain(f)
    A = object(algebra(parent(X)))
    r = right_action(X)
    C,c = cokernel(morphism(f))
    inv_c = right_inverse(c)
    coker_action = c ∘ r ∘ (inv_c ⊗ id(A))
    coker = RightModuleObject(parent(f), C, coker_action)
    return coker, ModuleMorphism(X, coker, c)
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
#=----------------------------------------------------------
    Free modules   
----------------------------------------------------------=#

function free_left_module(X::Object, A::AlgebraObject, parent = LeftModuleCategory(parent(X), A))
    AA = object(A)
    M = AA ⊗ X
    l = compose(
        inv_associator(AA, AA, X),
        multiplication(A) ⊗ id(X)
    )
    LeftModuleObject(parent, M, l)
end

function free_right_module(X::Object, A::AlgebraObject, parent = RightModuleCategory(parent(X), A))
    AA = object(A)
    M = X ⊗ AA
    l = compose(
        associator(X, AA, AA),
        id(X) ⊗ multiplication(A)
    )
    RightModuleObject(parent, M, l)
end

function free_bimodule(X::Object, A::AlgebraObject, parent = BiModuleCategory(parent(X), A, A))
    free_bimodule(X,A,A,parent)
end

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

function free_module(M::ModuleCategory, X::Object)
    if typeof(M) == LeftModuleCategory
        return free_left_module(X,left_algebra(M), M)
    elseif typeof(M) == RightModuleCategory
        return free_right_module(X, right_algebra(M), M)
    else
        return free_bimodule(X, left_algebra(M), right_algebra(M), M)
    end
end


@doc raw""" 

    left_module(M::RightModuleObject)

Construct the left module ``(∗M, q)`` from ``(M,p)``. 
"""
function left_module(N::RightModuleObject)
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

    C = RightModuleCategory(
        parent(M),
        algebra(parent(N))
    )

    LeftModuleObject(C, dM, q)
end

function right_module(N::LeftModuleObject)
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

    C = LeftModuleCategory(
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

function tensor_product(f::ModuleMorphism{BiModuleObject}, g::ModuleMorphism{BiModuleObject})
    dom, p_dom = bimodule_tensor_product(domain(f), domain(g))
    cod, p_cod = bimodule_tensor_product(codomain(f), codomain(g))

    inv_p_dom = right_inverse(p_dom)

    ModuleMorphism(dom, cod, p_cod ∘ (morphism(f) ⊗ morphism(g)) ∘ inv_p_dom)
end

function associator(X::BiModuleObject, Y::BiModuleObject, Z::BiModuleObject)
    XY, p_XY = bimodule_tensor_product(X,Y)
    XY_Z, p_XY_Z = bimodule_tensor_product(XY, Z)

    YZ, p_YZ = bimodule_tensor_product(Y, Z)
    X_YZ, p_X_YZ = bimodule_tensor_product(X, YZ)

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
    Morphism(XY_Z, X_YZ, a)
end

function one(C::BiModuleCategory)
    @assert left_algebra(C) == right_algebra(C)
    A = left_algebra(C)
    bimodule(A)
end

function dual(M::BiModuleObject)
    @assert left_algebra(parent(M)) == right_algebra(parent(M))

    lM = left_module(right_module(M))
    rM = right_module(left_module(M))
    
    BiModuleObject(
        parent(M),
        object(lM),
        left_action(lM),
        right_action(rM)
    )
end

function ev(M::BiModuleObject)
    dM = dual(M)
    e = ev(right_module(M), right_module(one(parent(M))))
    dMM,p = bimodule_tensor_product(dM,M)
    inv_p = right_inverse(p)
    Morphism(dMM, one(parent(M)), morphism(e) ∘ inv_p)
end

function coev(M::BiModuleObject)
    dM = dual(M)
    MdM, p = bimodule_tensor_product(M,dM)

    base = basis(Hom(one(parent(M)), MdM))
    M_basis = basis(End(M))

    m = (id(M)⊗ev(M))∘associator(M,dM,M)

    M = hcat([express_in_basis(m∘(f ⊗ id(M)), M_basis) for f ∈ base]...)
    one_coeffs = express_in_basis(id(one(parent(dM))), M_basis)

    s = solve(transpose(matrix(base_ring(dM), size(M,1), size(M,2), M)), one_coeffs)

    sum(s .* base)
end

function spherical(M::BiModuleObject)
    ddM = dual(dual(M))
    Morphism(M, ddM, spherical(object(M)))
end

#=----------------------------------------------------------
    Internal Hom 
----------------------------------------------------------=#

function internal_hom(M::RightModuleObject, N::RightModuleObject)
    dual(M ⊗ left_module(N))
end

function internal_hom_adjunction(X::Object, M::RightModuleObject, N::RightModuleObject, f::Morphism)

    A = algebra(parent(M))
    p = right_action(M)
    dN = left_module(N)
    q = left_action(dN)
    
    # Build the internal hom with the projection 
    # M ⊗ ∗N → M ⊗ₐ ∗N 
    H,h = _tensor_product(M, left_module(N))

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
    dN = left_module(N)
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

⊗(X::Object, M::RightModuleObject) = left_module_action(X, M)
⊗(M::LeftModuleObject, X::Object) = right_module_action(X, M)

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

    module_hom_basis = [Morphism(X,Y,sum(b .* B)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,module_hom_basis)
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

    module_hom_basis = [Morphism(X,Y,sum(b .* B)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,module_hom_basis)
end

function Hom(X::BiModuleObject, Y::BiModuleObject)
    H = Hom(object(X), object(Y))
    B_H = basis(H)
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

    module_hom_basis = ModuleMorphism[Morphism(X,Y,sum(b .* B_H)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,module_hom_basis)
end

#=----------------------------------------------------------
    isomorphic 
----------------------------------------------------------=#

function is_isomorphic(M::ModuleObject, N::ModuleObject)
    if is_simple(M) && is_simple(N)
        return int_dim(Hom(M,N)) > 0
    else
        error("not implemented yet")
    end
end


#=----------------------------------------------------------
    Compute Module Categories 
----------------------------------------------------------=#

is_semisimple(M::LeftModuleCategory) = is_separable(left_algebra(M))
is_semisimple(M::RightModuleCategory) = is_separable(right_algebra(M))

function category_of_right_modules(A::AlgebraObject)
    M = RightModuleCategory(parent(object(A)), A)
end

function category_of_left_modules(A::AlgebraObject)
    M = LeftModuleCategory(parent(object(A)), A)
end

function category_of_bimodules(A::AlgebraObject, B::AlgebraObject)
    M = BiModuleCategory(parent(object(A)), A,B)
end

function category_of_bimodules(A::AlgebraObject)
    category_of_bimodules(A,A)
end

function simples(M::ModuleCategory)
    if isdefined(M, :simples)
        return M.simples
    end

    C = category(M)

    free_objects = [free_module(M,s) for s ∈ simples(C)]

    simpls = unique_simples(vcat([simple_subobjects(x) for x ∈ free_objects]...))

    M.simples = simpls
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