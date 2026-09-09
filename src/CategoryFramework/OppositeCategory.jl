# EGNO, Tensor Categories (2015), §1.1: Hom_{C^op}(X,Y)=Hom_C(Y,X).
# compose(f,g) means g ∘ f in this package, so op reverses the argument list.
struct OppositeCategory <: Category
    C::Category
end

struct OppositeObject <: Object
    parent::OppositeCategory
    X::Object
end

struct OppositeMorphism <: Morphism
    domain::OppositeObject
    codomain::OppositeObject
    m::Morphism
end

op(C::Category) = OppositeCategory(C)
op(C::OppositeCategory) = C.C
op(X::Object) = op(parent(X))(X)
op(X::OppositeObject) = object(X)
op(f::Morphism) = op(parent(f))(f)
op(f::OppositeMorphism) = morphism(f)

"Construct the opposite category, reversing the direction of every arrow."
opposite_category(C::Category) = op(C)
"Regard X as an object in the opposite category."
opposite_object(X::Object) = op(X)
"Regard f : X → Y as an arrow Y → X in the opposite category."
opposite_morphism(f::Morphism) = op(f)

function (C::OppositeCategory)(X::Object)
    parent(X) == C.C || throw(ArgumentError("object belongs to a different category"))
    OppositeObject(C,X)
end
(C::OppositeCategory)(f::Morphism) = morphism(C(codomain(f)), C(domain(f)), f)

function morphism(X::OppositeObject, Y::OppositeObject, f::Morphism)
    parent(X) == parent(Y) || throw(ArgumentError("mismatching opposite categories"))
    domain(f) == object(Y) && codomain(f) == object(X) ||
        throw(ArgumentError("underlying morphism must reverse the specified endpoints"))
    OppositeMorphism(X,Y,f)
end

base_ring(C::OppositeCategory) = base_ring(C.C)
parent(X::OppositeObject) = X.parent
object(X::OppositeObject) = X.X
morphism(f::OppositeMorphism) = f.m
object_type(::OppositeCategory) = OppositeObject
morphism_type(::OppositeCategory) = OppositeMorphism
==(C::OppositeCategory,D::OppositeCategory) = C.C == D.C
==(X::OppositeObject,Y::OppositeObject) = parent(X) == parent(Y) && object(X) == object(Y)
==(f::OppositeMorphism,g::OppositeMorphism) = domain(f) == domain(g) && codomain(f) == codomain(g) && morphism(f) == morphism(g)

function compose(fs::OppositeMorphism...)
    isempty(fs) && throw(ArgumentError("at least one morphism is required"))
    all(codomain(fs[i]) == domain(fs[i+1]) for i in 1:length(fs)-1) ||
        throw(ArgumentError("morphisms are not composable"))
    morphism(domain(first(fs)), codomain(last(fs)), compose(reverse(morphism.(fs))...))
end
id(X::OppositeObject) = morphism(X,X,id(object(X)))
inv(f::OppositeMorphism) = morphism(codomain(f),domain(f),inv(morphism(f)))
zero_morphism(X::OppositeObject,Y::OppositeObject) = morphism(X,Y,zero_morphism(object(Y),object(X)))
is_zero(f::OppositeMorphism) = is_zero(morphism(f))
is_zero(X::OppositeObject) = is_zero(object(X))
*(a,f::OppositeMorphism) = morphism(domain(f),codomain(f),a*morphism(f))
function +(f::OppositeMorphism,g::OppositeMorphism)
    domain(f) == domain(g) && codomain(f) == codomain(g) || throw(ArgumentError("mismatching endpoints"))
    morphism(domain(f),codomain(f),morphism(f)+morphism(g))
end
matrix(f::OppositeMorphism) = transpose(matrix(morphism(f)))
express_in_basis(f::OppositeMorphism,B::Vector{<:OppositeMorphism}) = express_in_basis(morphism(f),morphism.(B))

function product(X::OppositeObject,Y::OppositeObject)
    Z,incs = coproduct(object(X),object(Y))
    parent(X)(Z), parent(X).(incs)
end
function coproduct(X::OppositeObject,Y::OppositeObject)
    Z,projs = product(object(X),object(Y))
    parent(X)(Z), parent(X).(projs)
end
function direct_sum(X::OppositeObject,Y::OppositeObject)
    Z,incs,projs = direct_sum(object(X),object(Y))
    parent(X)(Z), parent(X).(projs), parent(X).(incs)
end
function kernel(f::OppositeMorphism)
    X,p = cokernel(morphism(f))
    parent(f)(X),parent(f)(p)
end
function cokernel(f::OppositeMorphism)
    X,i = kernel(morphism(f))
    parent(f)(X),parent(f)(i)
end

function tensor_product(X::OppositeObject,Y::OppositeObject)
    parent(X) == parent(Y) || throw(ArgumentError("mismatching opposite categories"))
    parent(X)(object(X) ⊗ object(Y))
end
tensor_product(f::OppositeMorphism,g::OppositeMorphism) = morphism(domain(f)⊗domain(g),codomain(f)⊗codomain(g),morphism(f)⊗morphism(g))
# The associator must still point (X⊗Y)⊗Z → X⊗(Y⊗Z) in C^op.
associator(X::OppositeObject,Y::OppositeObject,Z::OppositeObject) = parent(X)(inv(associator(object(X),object(Y),object(Z))))
one(C::OppositeCategory) = C(one(C.C))
zero(C::OppositeCategory) = C(zero(C.C))
simples(C::OppositeCategory) = C.(simples(C.C))
indecomposables(C::OppositeCategory) = C.(indecomposables(C.C))
decompose(X::OppositeObject) = [(parent(X)(Y),m) for (Y,m) in decompose(object(X))]
function is_isomorphic(X::OppositeObject,Y::OppositeObject)
    parent(X) == parent(Y) || return false,nothing
    ok,f = is_isomorphic(object(Y),object(X))
    ok ? (true,morphism(X,Y,f)) : (false,nothing)
end

is_semisimple(C::OppositeCategory) = is_semisimple(C.C)
is_abelian(C::OppositeCategory) = is_abelian(C.C)
is_additive(C::OppositeCategory) = is_additive(C.C)
is_linear(C::OppositeCategory) = is_linear(C.C)
is_monoidal(C::OppositeCategory) = is_monoidal(C.C)
is_finite(C::OppositeCategory) = is_finite(C.C)
is_locally_finite(C::OppositeCategory) = is_locally_finite(C.C)
is_tensor(C::OppositeCategory) = is_tensor(C.C)
is_multitensor(C::OppositeCategory) = is_multitensor(C.C)
is_fusion(C::OppositeCategory) = is_fusion(C.C)
is_multifusion(C::OppositeCategory) = is_multifusion(C.C)
is_krull_schmidt(C::OppositeCategory) = is_krull_schmidt(C.C)

function Hom(X::OppositeObject,Y::OppositeObject)
    parent(X) == parent(Y) || throw(ArgumentError("mismatching opposite categories"))
    HomSpace(X,Y,OppositeMorphism[parent(X)(f) for f in basis(Hom(object(Y),object(X)))])
end
