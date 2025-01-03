#=----------------------------------------------------------
    Define G-actions on Tensor Categories
----------------------------------------------------------=#


struct GTensorAction 
    C::Category
    G::GAPGroup
    images::Vector{<:AbstractFunctor}
end

group(T::GTensorAction) = T.G

function gtensor_action(C::Category, G::GAPGroup, images::Vector{<:Functor})
    GTensorAction(C,G,images)
end

function gtensor_action(C::Category, G::GAPGroup, images::Vector{<:Object})
    GTensorAction(C,G,[inner_autoequivalence(C,X) for X ∈ images])
end

function (T::GTensorAction)(g::GroupElem)
    i = findfirst(==(g), elements(group(T)))
    return T.images[i]
end


#=----------------------------------------------------------
    Cannonical G-action
----------------------------------------------------------=#

function gtensor_action(C::Category, G::GAPGroup)

    C₁ = invertibles(C)
    n = length(C₁)

    mult = multiplication_table(C₁)
    M = [findfirst(!iszero, mult[i,j,:]) for i ∈ 1:n, j ∈ 1:n]

    H = MultTableGroup(M)

    m,r = divrem(order(G), n)

    if r != 0 
        return gcrossed_product(C, trivial_gtensor_action(C,G))
    end

    subs = representative.(subgroup_classes(G, order = m))

    i = findfirst(N -> is_isomorphic(quo(G,N)[1], H), subs)

    Q,p = quo(G,subs[i])
     
    proj = compose(p, is_isomorphic_with_map(Q, H)[2])

    els = elements(H)

    images = [findfirst(==(proj(g)), els) for g ∈ elements(G)]
    images = [C₁[i] for i ∈ images]

    action = gtensor_action(C, G, images)
end

#=----------------------------------------------------------
    Trivial G-action 
----------------------------------------------------------=#

