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
    categegory::Category
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


#=----------------------------------------------------------
    Constructors 
----------------------------------------------------------=#

function LeftModule(A::AlgebraObject)
    C = LeftModuleCategory(parent(A), A)
    LeftModuleObject(C, object(A), multiplication(A))
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