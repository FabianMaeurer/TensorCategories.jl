#=----------------------------------------------------------
    Construct algebra objects in multitensor categories
    
    Reference: https://math.mit.edu/~etingof/egnobookfinal.pdf
----------------------------------------------------------=#

@attributes mutable struct AlgebraObject <: Object
    parent::Category
    object::Object
    multiplication::Morphism
    unit::Morphism

    function AlgebraObject(C::Category, X::Object, m::Morphism, u::Morphism)
        new(C,X,m,u)
    end
end

struct AlgebraMorphism <: Morphism 
    domain::AlgebraObject
    codomain::AlgebraObject
    map::Morphism
end

function morphism(X::AlgebraObject, Y::AlgebraObject, f::Morphism)
    AlgebraMorphism(X,Y,f)
end
#=----------------------------------------------------------
    Getter 
----------------------------------------------------------=#

multiplication(A::AlgebraObject) = A.multiplication
unit(A::AlgebraObject) = A.unit
object(A::AlgebraObject) = A.object

#=----------------------------------------------------------
    Check 
----------------------------------------------------------=#

function is_algebra(X::Object, m::Morphism, u::Morphism)
    m ∘ (id(X) ⊗ u) == m ∘ (u ⊗ id(X)) == id(X) &&
    m ∘ (id(X) ⊗ m) ∘ associator(X,X,X) == m ∘ (m ⊗ id(X))
end

is_algebra(X::AlgebraObject) = is_algebra(object(X), multiplication(X), unit(X))

function is_separable(A::AlgebraObject)
    get_attribute!(A, :is_separable) do
        C = parent(A)
        m = multiplication(A)

        # Define the multiplication as a bimodule morphism
        AA = free_bimodule(one(C), A)
        m = morphism(AA, bimodule(A), multiplication(A))

        # A is seperable if m has a right inverse as bimodule morphism
        has_right_inverse(m)
    end
end

function is_commutative(A::AlgebraObject)
    get_attribute!(A, :is_commutative) do 
        m = multiplication(A)
        m == m ∘ braiding(object(A),object(A))
    end
end
#=----------------------------------------------------------
    Group Algebras
----------------------------------------------------------=#

function group_algebra(C::Category, G::Group)
    @assert is_multitensor(C)
    KG = QQ[G]
    S,i,p = direct_sum([one(C) for _ ∈ 1:dim(KG)])

    mult = multiplication_table(KG)

    m = sum([i[mult[m,n]] ∘ (p[m]⊗p[n]) for m ∈ 1:dim(KG), n ∈ 1:dim(KG)])
    u = i[1]

    AlgebraObject(C, S, m, u)
end

#=----------------------------------------------------------
    Morita equivalence 
----------------------------------------------------------=#

function is_morita_equivalent(A::AlgebraObject, B::AlgebraObject)
    M = category_of_right_modules(A)
    N = category_of_right_modules(B)
    is_equivalent(M,N)
end

#=----------------------------------------------------------
    Generic Algebras 
----------------------------------------------------------=#

function generic_algebra(X::Object)
    dX = dual(X)
    A = X ⊗ dX
    m = compose(
        associator(X, dX, A),
        id(X) ⊗ inv_associator(dX,X,dX),
        id(X) ⊗ (ev(X) ⊗ id(dX))
    )
    return AlgebraObject(parent(X), A, m, coev(X))
end

#=----------------------------------------------------------
    extension of scalars  
----------------------------------------------------------=#

function extension_of_scalars(A::AlgebraObject, K::Ring)
    C = extension_of_scalars(parent(A), K)
    X = extension_of_scalars(object(A), K)
    m = extension_of_scalars(multiplication(A), K)
    u = extension_of_scalars(unit(A), K)
    AlgebraObject(C,X,m,u)
end
#=----------------------------------------------------------
    Pretty printing
----------------------------------------------------------=#

function show(io::IO, A::AlgebraObject)
    print(io, """Algebra object $(object(A))""")
end