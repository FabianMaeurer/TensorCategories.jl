#=----------------------------------------------------------
    Define G-actions on Tensor Categories
----------------------------------------------------------=#


struct GTensorAction <: AbstractMonoidalFunctor
    category::Category
    G::Group
    elements::Vector{<:GroupElem}
    images::Vector{<:AbstractFunctor}
    monoidal_structure::Dict
end

group(T::GTensorAction) = T.G
domain(T::GTensorAction) = group(T) 


@doc raw""" 

    gtensor_action(C::Category, elems::Vector{<:GroupElem}, images::Vector{<:AbstractFunctor}, monoidal_structure::Dict)

Define an action of ```G``` on ```C``` by providing an image functor for every element in ```G```.
"""
function gtensor_action(C::Category, elems::Vector{<:GroupElem}, images::Vector{<:AbstractFunctor}, monoidal_structure::Dict)
    GTensorAction(C, parent(elems[1]), elems, images, monoidal_structure)
end

# function gtensor_action(C::Category, G::Group, images::Vector{<:Object})
#     GTensorAction(C,G,[inner_autoequivalence(C,X) for X ∈ images])
# end

function (T::GTensorAction)(g::GroupElem)
    i = findfirst(==(g), T.elements)
    return T.images[i]
end

images(T::GTensorAction) = T.images
category(T::GTensorAction) = T.category

function monoidal_structure(T::GTensorAction, g::GroupElem, h::GroupElem) 
    i = findfirst(==(g), T.elements)
    j = findfirst(==(h), T.elements)
    
    T.monoidal_structure[(i,j)]
end
#=----------------------------------------------------------
    Cannonical G-action
----------------------------------------------------------=#

function action_by_inner_autoequivalences(C::Category)

    # Build Group of inner autoequivalences
    inner_autos = inner_autoequivalences(C)

    nat_trafos = [[monoidal_natural_transformations(F∘G,e) for e ∈ inner_autos] for F ∈ inner_autos, G ∈ inner_autos]

    mult = [findfirst(e -> length(e) > 0, n) for n ∈ nat_trafos]

    M = MultTableGroup(mult)

    n = length(inner_autos)
    monoidal_structure = Dict()

    indecs = indecomposables(inner_autos[1])

    for i ∈ 1:n, j ∈ 1:n 
        F,G = inner_autos[[i,j]]
        Y = G.object 
        X = F.object 
        dY = G.dual_object 
        dX = F.dual_object 

        dual_iso = inv(dual_monoidal_structure(X,Y))
        m = [compose( 
            inv_associator(Y,X⊗V,dX) ⊗ id(dY),
            (inv_associator(Y,X,V) ⊗ id(dX)) ⊗ id(dY),
            associator((Y⊗X)⊗V,dX,dY),
            id((Y⊗X)⊗V) ⊗ dual_iso
        ) for V ∈ indecs]

        monoidal_structure[(i,j)] = AdditiveNaturalTransformation(
            F∘G,
            inner_autos[mult[i,j]],
            indecs,
            m
        ) 
    end

    action = gtensor_action(C, elements(M), inner_autos, monoidal_structure)
end

@doc raw""" 

    inner_automorphisms(C::Category)

Return a vector with all non-equivalent inner automorphisms
"""
function inner_autoequivalences(C::Category) 
    invertibls = invertibles(C)

    inner_autos = inner_autoequivalence.(invertibls)

    unique_autos = [inner_autos[1]]

    for F ∈ inner_autos[2:end]
        if sum(length.([monoidal_natural_transformations(F,e) for e ∈ unique_autos])) == 0
            push!(unique_autos, F)
        end
    end
    unique_autos
end

#=----------------------------------------------------------
    Trivial G-action 
----------------------------------------------------------=#


#=----------------------------------------------------------
    Is Action? 
----------------------------------------------------------=#

function is_tensor_action(F::GTensorAction)
    S = F.elements
    
    for X ∈ S, Y ∈ S, Z ∈ S 
        left = compose(
            monoidal_structure(F,X,Y) ⊗ id(F(Z)),
            monoidal_structure(F, X * Y, Z)
        )

        right = compose(
            id(F(X)) ⊗ monoidal_structure(F,Y,Z),
            monoidal_structure(F, X, Y * Z)
        )

        if right != left 
            return false 
        end
    end

    true
end
#=----------------------------------------------------------
    Printing 
----------------------------------------------------------=#

function show(io::IO, T::GTensorAction)
    print(io, "Action of $(group(T)) on $(category(T))")
end