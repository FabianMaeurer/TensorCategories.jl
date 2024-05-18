#=----------------------------------------------------------
    Structures for left, right and bimodules over 
    algebra objects in multitensor categories
----------------------------------------------------------=#

abstract type ModuleCategory <: Category end
abstract type ModuleObject <: Object end

struct LeftModuleCategory <: ModuleCategory
    category::Category
    left_algebra::AlgebraObject
end

struct RightModuleCategory <: ModuleCategory
    category::Category
    right_algebra::AlgebraObject
end

struct BiModuleCategory <: ModuleCategory
    categegory::Category
    left_algebra::AlgebraObject
    right_algebra::AlgebraObject
end

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

#=----------------------------------------------------------
    Functionality 
----------------------------------------------------------=#

id(X::ModuleObject) = ModuleMorphism(X,X, id(object(X)))

function kernel(f::ModuleMorphism)
    X = domain(f)
    typeof(X) == LeftModuleObject  && return left_module_kernel(f)
    typeof(X) == RightModuleObject && return right_module_kernel(f)
    typeof(X) == BiModuleObject    && return bi_module_kernel(f)
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

function cokernel(f::ModuleMorphism)
    X = domain(f)
    typeof(X) == LeftModuleObject  && return left_module_cokernel(f)
    typeof(X) == RightModuleObject && return right_module_cokernel(f)
    typeof(X) == BiModuleObject    && return bi_module_cokernel(f)
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

#=----------------------------------------------------------
    Internal Hom 
----------------------------------------------------------=#

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


function tensor_product(M::RightModuleObject, N::LeftModuleObject)
    @assert algebra(parent(M)) == algebra(parent(N))

    A = algebra(parent(M))
    p = right_action(M)
    q = left_action(N)

    C,c = coequilizer(
        p ⊗ id(object(N)),
        (id(object(M)) ⊗ q) ∘ associator(object(M), object(A), object(N))
    )

    return C
end


function internal_hom(M::RightModuleObject, N::RightModuleObject)
    dual(M ⊗ left_module(N))
end

function internal_hom_adjunction(X::Object, M::RightModuleObject, N::RightModuleObject, f::Morphism)

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

⊗(X::Object, M::RightModuleObject) = left_module_action(X, M)

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

    base = basis(Hom(object(A)⊗object(X), object(Y)))

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

#=----------------------------------------------------------
    Pretty printing
----------------------------------------------------------=#

function show(io::IO, C::ModuleCategory)
    typeof(C) == LeftModuleCategory && 
        print(io, """Category of left modules over $(algebra(C))""")
    typeof(C) == RightModuleCategory && 
        print(io, """Category of right modules over $(algebra(C))""")
    typeof(C) == BiModuleCategory && 
        print(io, """Category of bimodules over $(algebra(C))""")
end

function show(io::IO, X::ModuleObject)
    typeof(X) == LeftModuleObject && 
        print(io, """Left module: $(object(X))""")
    typeof(X) == RightModuleObject && 
        print(io, """Right module: $(object(X))""")
    typeof(X) == BiModuleObject && 
        print(io, """Bimodule: $(object(X))""")
end