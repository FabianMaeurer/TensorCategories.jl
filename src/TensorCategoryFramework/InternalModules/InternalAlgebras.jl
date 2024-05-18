#=----------------------------------------------------------
    Construct algebra objects in multitensor categories
    
    Reference: https://math.mit.edu/~etingof/egnobookfinal.pdf
----------------------------------------------------------=#

struct AlgebraObject <: Object
    parent::Category
    object::Object
    multiplication::Morphism
    unit::Morphism
end

struct AlgebraMorphism <: Morphism 
    domain::AlgebraObject
    codomain::AlgebraObject
    map::Morphism
end

function Morphism(X::AlgebraObject, Y::AlgebraObject, f::Morphism)
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
    m ∘ (id(X) ⊗ m) == m ∘ (m ⊗ id(X)) &&
    m ∘ (id(X) ⊗ u) == m ∘ (u ⊗ id(X)) == id(X) &&
    m ∘ (id(X) ⊗ m) ∘ associator(X,X,X) == m ∘ (m ⊗ id(X))
end

is_algebra(X::AlgebraObject) = is_algebra(object(X), multiplication(X), unit(X))

#=----------------------------------------------------------
    Group Algebras
----------------------------------------------------------=#

function group_algebra(C::Category, G::GAPGroup)
    @assert is_multitensor(C)
    KG = QQ[G]
    S,i,p = direct_sum([one(C) for _ ∈ 1:dim(KG)])

    mult = multiplication_table(KG)

    m = sum([i[mult[m,n]] ∘ (p[m]⊗p[n]) for m ∈ 1:dim(KG), n ∈ 1:dim(KG)])
    u = i[1]

    AlgebraObject(C, S, m, u)
end

#=----------------------------------------------------------
    Generic Algebras 
----------------------------------------------------------=#

function generic_algebra(X::Object)
    dX = dual(X)
    A = X ⊗ dX
    m = compose(
        associator(X ⊗ dX, X, dX),
        inv_associator(X,dX,X) ⊗ id(dX),
        (id(X) ⊗ ev(X)) ⊗ id(dX)
    )
    return AlgebraObject(parent(X), A, m, coev(X))
end
#=----------------------------------------------------------
    Pretty printing
----------------------------------------------------------=#

function show(io::IO, A::AlgebraObject)
    print(io, """Algebra object $(object(A))""")
end