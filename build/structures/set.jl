
struct Sets <: Category end

struct SetCategoryObject <: CategoryObject
    set::T where T <: AbstractSet
end

struct SetCategoryMorphism <: CategoryMorphism
    m::Dict
    domain::SetCategoryObject
    codomain::SetCategoryObject
end

#--------------------------------------------------
#   Constructors
#--------------------------------------------------

SetCategoryObject(S::Array) = SetCategoryObject(Set(S))

function SetCategoryMorphism(D::SetCategoryObject, C::SetCategoryObject, m::Dict)
    if keys(m) ⊆ D && values(m) ⊆ C
        return SetCategoryMorphism(m,D,C)
    else
        throw(ErrorException("Mismatching (co)domain"))
    end
end

SetCategoryMorphism(D::SetCategoryObject, C::SetCategoryObject, m::Function) = SetCategoryMorphism(D,C, Dict(x => m(x) for x ∈ D))

#--------------------------------------------------
#   Functionality
#--------------------------------------------------

in(item,S::SetCategoryObject) = in(item,S.set)

issubset(item, S::SetCategoryObject) = issubset(item, S.set)

iterate(X::SetCategoryObject) = iterate(X.set)
iterate(X::SetCategoryObject,state) = iterate(X.set,state)

length(X::SetCategoryObject) = length(X.set)

(f::SetCategoryMorphism)(item) = f.m[item]

==(X::SetCategoryObject,Y::SetCategoryObject) = X.set == X.set

parent(X::SetCategoryObject) = Sets()

#--------------------------------------------------
#   Functionality: CategoryMorphisms
#--------------------------------------------------

function compose(f::SetCategoryMorphism...)
    if length(f) == 1 return f[1] end
    if [domain(f[i]) == codomain(f[i-1]) for i ∈ 2:length(f)] != trues(length(f)-1)
        throw(ErrorException("Morphisms not compatible"))
    end
    m = f[1]
    for g in f[2:end]
        m = Dict(x => g(m(x)) for x ∈ keys(m.m))
    end
    return SetCategoryMorphism(domain(f[1]), codomain(f[end]),m)
end


function inv(f::SetCategoryMorphism)
    if length(values(f.m)) == length(keys(f.m))
        SetCategoryMorphism( codomain(f),domain(f), Dict(v => k for (k,v) ∈ f.m))
    else
        throw(ErrorException("Not invertible"))
    end
end

id(X::SetCategoryObject) = SetCategoryMorphism(X,X, x->x)
==(f::SetCategoryMorphism, g::SetCategoryMorphism) = f.m == g.m

#--------------------------------------------------
#   Product
#--------------------------------------------------

function product(X::SetCategoryObject, Y::SetCategoryObject, projections = false)
    Z = SetCategoryObject(Set([(x,y) for x ∈ X, y ∈ Y]))
    pX = SetCategoryMorphism(Z,X, x -> x[1])
    pY = SetCategoryMorphism(Z,Y, x -> x[2])
    return projections ? (Z,[pX,pY]) : Z
end

function coproduct(X::SetCategoryObject, Y::SetCategoryObject, injections = false)
    if length(X.set ∩ Y.set) != 0
        Z = SetCategoryObject(union([(x,0) for x ∈ X],[(y,1) for y ∈ Y]))
        ix = SetCategoryMorphism(X,Z, x -> (x,0))
        iy = SetCategoryMorphism(Y,Z, y -> (y,1))
    else
        Z = SetCategoryObject(union(X.set,Y.set))
        ix = SetCategoryMorphism(X,Z, x -> x)
        iy = SetCategoryMorphism(Y,Z, y -> y)
    end
    return injections ? (Z, [ix,iy]) : Z
end

#--------------------------------------------------
#   CategoryHomSets
#--------------------------------------------------

struct SetCategoryHomSet <: CategoryHomSet
    X::SetCategoryObject
    Y::SetCategoryObject
end

Hom(X::SetCategoryObject, Y::SetCategoryObject) = SetCategoryHomSet(X,Y)

#--------------------------------------------------
#   Pretty printing
#--------------------------------------------------

function show(io::IO, X::Sets)
    print(io,"Category of finte sets")
end

function show(io::IO, X::SetCategoryObject)
    print(X.set)
end
