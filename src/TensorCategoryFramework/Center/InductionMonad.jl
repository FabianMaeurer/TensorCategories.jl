#=----------------------------------------------------------
    Tools to compute the center as algebras over the 
    induction monad. 
----------------------------------------------------------=#

mutable struct InductionMonad <: Monad 
    codomain::Category 
    multiplication::Dict{<:Object, <:Morphism}

    function InductionMonad(C::Category)
        T = new()
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
        (inv(half_braiding(TX,x)) ⊗ id(dual(x))),
        associator(object(TX),x,dual(x)),
        (id(object(TX)) ⊗ ev(dual(x)))
    ) for x ∈ simples(parent(X))])
end

#=----------------------------------------------------------
    Maths 
----------------------------------------------------------=#    

@doc raw""" 

    induction_monad(V::Object)

Compute the induction monad together with the structure morphism
`μ : TT(V) → T(V)``.

Formula according to  
https://doi.org/10.48550/arXiv.1203.4180
"""
function induction_monad(V::Object)
    C = parent(V)

    S = simples(C)
    DS = dual.(S)

    SdS = zip(S,DS)

    dom = direct_sum([(J ⊗ ((I ⊗ V) ⊗ I̅)) ⊗ J̅ for (I,I̅) ∈ SdS, (J,J̅) ∈ SdS][:])[1]

    cod = direct_sum([(K ⊗ V) ⊗ dK for (K,dK) ∈ SdS])[1]
    
    component_maps = Morphism[]

    for (J,J̅) ∈ SdS, (I,I̅) ∈ SdS
        dom_JI = (J ⊗ ((I ⊗ V) ⊗ I̅)) ⊗ J̅
        JI_component = zero_morphism(dom_JI, zero(C))

        _,_,i,p = direct_sum_decomposition(J⊗I)

        a_IJ = compose(
            inv_associator(J, I⊗V, I̅) ⊗ id(J̅),
            associator(J⊗(I⊗V), I̅, J̅),
            inv_associator(J,I,V) ⊗ id(I̅⊗J̅)
        )

        for (K,dK) ∈ SdS
            base = [f for f ∈ p if codomain(f) == K]
            dual_base = dual.([f for f ∈ i if domain(f) == K])

            if length(base) == 0 
                JIK = zero_morphism(dom_JI, K ⊗ V ⊗ dK)
            else
                JIK = sum((f ⊗ id(V)) ⊗ g for (f,g) ∈ zip(base, dual_base))
            end

            JI_component = vertical_direct_sum(JI_component, JIK)
        end

        push!(component_maps, JI_component ∘ a_IJ)
    end

    μ = horizontal_direct_sum(component_maps)

    IVI̅ = [I ⊗ V ⊗ I̅ for (I,I̅) ∈ SdS]
    _JIVI̅ = [[J ⊗ (I ⊗ V ⊗ I̅) for (I,I̅) ∈ SdS] for J ∈ S]


    distr_before = compose(
        direct_sum([distribute_right(J, IVI̅) ⊗ id(J̅) for (J,J̅) ∈ SdS]),
        direct_sum([distribute_left(JIVI̅, J̅) for (JIVI̅,J̅) ∈ zip(_JIVI̅,DS)])
    )
    
    return μ ∘ distr_before
    # μ = horizontal_direct_sum([
    #     compose(
    #         inv(half_braiding(TV,Y)) ⊗ id(dual(Y)),
    #         associator(object(TV), Y, dual(Y)),
    #         id(object(TV)) ⊗ (ev(dual(Y)) ∘ (spherical(Y) ⊗ id(dual(Y))))
    #     ) for Y ∈ simples(C)
    # ])
    # object(TV), μ
end

mutable struct Wedge 
    object::Object 
    simples_maps::Vector{Morphism}
    tensor_product_maps::Dict{Tuple{Object,Object}, Morphism}
end


#=----------------------------------------------------------
    Induction Coend 
----------------------------------------------------------=#

function induction_wedge(V::Object)
    indecs = indecomposables(parent(V))

end    