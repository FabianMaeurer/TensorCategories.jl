
struct Sets <: Category end

struct SetObject <: Object
    set::T where T <: AbstractSet
end

struct SetMorphism <: Morphism
    m::Dict
    domain::SetObject
    codomain::SetObject
end

#--------------------------------------------------
#   Constructors
#--------------------------------------------------

SetObject(S::Array) = SetObject(Set(S))

function SetMorphism(D::SetObject, C::SetObject, m::Dict)
    if keys(m) ⊆ D && values(m) ⊆ C
        return SetMorphism(m,D,C)
    else
        throw(ErrorException("Mismatching (co)domain"))
    end
end

SetMorphism(D::SetObject, C::SetObject, m::Function) = SetMorphism(D,C, Dict(x => m(x) for x ∈ D))

#--------------------------------------------------
#   Functionality
#--------------------------------------------------

in(item,S::SetObject) = in(item,S.set)

issubset(item, S::SetObject) = issubset(item, S.set)

iterate(X::SetObject) = iterate(X.set)
iterate(X::SetObject,state) = iterate(X.set,state)

length(X::SetObject) = length(X.set)

(f::SetMorphism)(item) = f.m[item]

==(X::SetObject,Y::SetObject) = X.set == X.set

parent(X::SetObject) = Sets()

#--------------------------------------------------
#   Functionality: Morphisms
#--------------------------------------------------

function compose(f::SetMorphism...)
    if length(f) == 1 return f[1] end
    if [domain(f[i]) == codomain(f[i-1]) for i ∈ 2:length(f)] != trues(length(f)-1)
        throw(ErrorException("Morphisms not compatible"))
    end
    m = f[1]
    for g in f[2:end]
        m = Dict(x => g(m(x)) for x ∈ keys(m.m))
    end
    return SetMorphism(domain(f[1]), codomain(f[end]),m)
end


function inv(f::SetMorphism)
    if length(values(f.m)) == length(keys(f.m))
        SetMorphism( codomain(f),domain(f), Dict(v => k for (k,v) ∈ f.m))
    else
        throw(ErrorException("Not invertible"))
    end
end

id(X::SetObject) = SetMorphism(X,X, x->x)
==(f::SetMorphism, g::SetMorphism) = f.m == g.m

#--------------------------------------------------
#   Product
#--------------------------------------------------

function product(X::SetObject, Y::SetObject)
    Z = SetObject(Set([(x,y) for x ∈ X, y ∈ Y]))
    pX = SetMorphism(Z,X, x -> x[1])
    pY = SetMorphism(Z,Y, x -> x[2])
    return Z,[pX,pY]
end

function coproduct(X::SetObject, Y::SetObject)
    if length(X.set ∩ Y.set) != 0
        Z = SetObject(union([(x,0) for x ∈ X],[(y,1) for y ∈ Y]))
        ix = SetMorphism(X,Z, x -> (x,0))
        iy = SetMorphism(Y,Z, y -> (y,1))
    else
        Z = SetObject(union(X.set,Y.set))
        ix = SetMorphism(X,Z, x -> x)
        iy = SetMorphism(Y,Z, y -> y)
    end
    return Z, [ix,iy]
end

#--------------------------------------------------
#   HomSets
#--------------------------------------------------

struct SetHomSet <: HomSet
    X::SetObject
    Y::SetObject
end

Hom(X::SetObject, Y::SetObject) = SetHomSet(X,Y)
