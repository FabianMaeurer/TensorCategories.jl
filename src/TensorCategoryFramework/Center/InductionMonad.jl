#=----------------------------------------------------------
    Tools to compute the center as algebras over the 
    induction monad. 
----------------------------------------------------------=#

mutable struct InductionMonad <: Monad 
    codomain::Category 
    multiplication::Dict{<:Object, <:Morphism}

    function InductionMonad(C::Category)
        T = New()
        T.category = C
        return T
    end 
end
 

function (T::InductionMonad)(X::Object)
    @assert parent(X) == codomain(T)
    induction_restriction(X)
end

function (T::InductionMonad)(f::Morphism)
    @assert parent(f) == codomain(T)
end

function domain(T::InductionMonad) 
    C = codomain(T)
    OppositeCategory(C) × C
end


function multiplication(T::InductionMonad, X::Object)
    if X ∈ keys(T.multiplication)
        return T.multiplication[X]
    end

    TX = induction(X)

    horizontal_direct_sum([dim(x) * compose(
        (inv(half_braiding(S[1],x)) ⊗ id(dual(x))),
        associator(object(S[1]),x,dual(x)),
        (id(object(S[1])) ⊗ ev(dual(x)))
    ) for x ∈ simples(parent(X))])
end

#=----------------------------------------------------------
    Maths 
----------------------------------------------------------=#    

function induction_monad(V::Object, TV = induction(V))
    C = category(parent(TV))

    μ = horizontal_direct_sum([
        compose(
            inv(half_braiding(TV,Y)) ⊗ id(dual(Y)),
            associator(object(TV), Y, dual(Y)),
            id(object(TV)) ⊗ (ev(dual(Y)) ∘ (spherical(Y) ⊗ id(dual(Y))))
        ) for Y ∈ simples(C)
    ])
    object(TV), μ
end

mutable struct Wedge 
    object::Object 
    simples_maps::Vector{Morphism}
    tensor_product_maps::Dict{Tuple{Object,Object}, Morphism}
end

