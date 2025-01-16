#=----------------------------------------------------------
    Compute Natrural Transformations in finitary categories. 
----------------------------------------------------------=#
abstract type NaturalTransformation <: Morphism end

mutable struct AdditiveNaturalTransformation <: NaturalTransformation
    domain::AbstractFunctor 
    codomain::AbstractFunctor 
    indecomposables::Vector{<:Object}
    maps::Vector{<:Morphism} 
end

function AdditiveNaturalTransformation(F::AbstractFunctor, G::AbstractFunctor, indecs::Vector{<:Object}, maps::Vector{<:Pair})
    @show [x ∈ keys(maps) for x in indecs]
    maps = [x ∈ keys(maps) ? maps[x] : zero_morphism(F(x),G(x)) for x ∈ indecs]
    AdditiveNaturalTransformation(F,G,indecs,maps)
end

struct NaturalTransformations <: AbstractHomSpace
    X::AbstractFunctor
    Y::AbstractFunctor
    basis::Vector{<:NaturalTransformation}
    parent::VectorSpaces
end

#=----------------------------------------------------------
    Functionality   
----------------------------------------------------------=#

function compose(η::AdditiveNaturalTransformation...)
    length(η) == 1 && return η[1]
    AdditiveNaturalTransformation(
        domain(η[1]),
        codomain(η[end]),
        indecomposables(η[1]),
        [compose([e.maps[i] for e ∈ η]...) for i ∈ 1:length(indecomposables(η[1]))]
    )
end

function inv(η::AdditiveNaturalTransformation)
    try 
        return AdditiveNaturalTransformation(
            domain(η),
            codomain(η),
            indecomposables(η),
            inv.(η.maps)
        )
    catch
        error("Not invertible")
    end
end

# function getindex(η::AdditiveNaturalTransformation, X::Object)
#     if X ∈ keys(η.maps)
#         return η.maps[X]
#     end

#     return zero_morphism(domain(η)(X), codomain(η(X)))
# end

function indecomposables(η::AdditiveNaturalTransformation)
    η.indecomposables
end

function id(F::AbstractFunctor)
    indecs = indecomposables(domain(F))
    AdditiveNaturalTransformation(
        F,
        F,
        indecs,
        [id(F(x)) for x ∈ indecs]
    )
end

function (η::AdditiveNaturalTransformation)(X::Object)
    indecs = indecomposables(η)
 
    i = findfirst(x -> x == X, indecs)

    i !== nothing && return η.maps[i]

    _,_,i,p = direct_sum_decomposition(X, indecs)
    F,G = domain(η), codomain(η)
    
    sum([G(iᵢ) ∘ η(domain(iᵢ)) ∘ F(pᵢ) for (iᵢ,pᵢ) ∈ zip(i,p)])
end

function *(x, η::AdditiveNaturalTransformation)
    AdditiveNaturalTransformation(
        domain(η),
        codomain(η),
        indecomposables(η),
        x .* η.maps
    )
end

function +(η::AdditiveNaturalTransformation, ν::AdditiveNaturalTransformation)
    S = indecomposables(η)
    T = indecomposables(ν)

    if S == T
        return AdditiveNaturalTransformation(
            domain(η),
            codomain(η),
            S,
            η.maps .+ ν.maps
        )
    end

    error("Not implemented")
end

#=----------------------------------------------------------
    Compute natural transformations 
----------------------------------------------------------=#

function Nat(F::AbstractFunctor, G::AbstractFunctor)
    @assert domain(F) == domain(G)
    @assert codomain(F) == codomain(G)

    if is_additive(F) && is_additive(G)
        nats = additive_natural_transformations(F,G)
        return NaturalTransformations(F,G,nats, VectorSpaces(base_ring(domain(F))))
    end

    @error("Cannot compute natural transformations")
end

function additive_natural_transformations(F::AbstractFunctor, G::AbstractFunctor, indecs = nothing)
    
    C = domain(F)
    K = base_ring(C)
    if indecs === nothing 
        indecs = indecomposables(C)
    end

    nats = NaturalTransformation[]

    for g ∈ group_indecomposables(indecs)
        # Bases for f: X → Y
        bases = [basis(Hom(s,t)) for s ∈ g, t ∈ g]

        # Bases for natural transform ηₓ: F(X) → G(X)
        nat_bases = [basis(Hom(F(s),G(s))) for s ∈ g]

        # Bases for coefficient comparison
        comp_bases = [basis(Hom(F(s), G(t))) for s ∈ g, t ∈ g]

        n = length(g)
        ns = length.(nat_bases)
        
        sum(ns) == 0 && continue

        Kx,x = polynomial_ring(K,sum(ns))

        x_blocks = [[popfirst!(x) for _ ∈ b] for b ∈ nat_bases]

        eqs = eltype(Kx)[]

        for (i,j) ∈ zip(1:n,1:n)
            B = bases[i,j]
            comp_B = comp_bases[i,j]
    
            xa,xb = x_blocks[[i,j]]
            nata,natb = nat_bases[[i,j]]

            e = [zero(Kx) for _ ∈ comp_B]
            
            for f ∈ B
                Gf,Ff = G(f),F(f)
                for (a,η) ∈ zip(xa, nata)
                    coeffs_a = express_in_basis(Gf ∘ η, comp_B)
                    e = e .+ (a .* coeffs_a) 
                end
                for (b,σ) ∈ zip(xb, natb)
                    coeffs_b = express_in_basis(σ ∘ Ff, comp_B)
                    e = e .- (b .* coeffs_b)
                end
            end

            eqs = eltype(Kx)[eqs; e]
        end
        
        x = union(x_blocks...)
        m = [coeff(e,a) for a ∈ x, e ∈ eqs]
        M = matrix(K, size(m,1), size(m,2), m)

        d,N = nullspace(M)

        solution_cols = [collect(c) for c ∈ eachcol(collect(N))]
        solutions = []

        for c ∈ solution_cols
            c_blocks = [[popfirst!(c) for _ ∈ b] for b ∈ nat_bases]
            maps = [sum(a .* B) for (a,B) ∈ zip(c_blocks, nat_bases)]
            push!(solutions, maps)
        end
        
        new_nats = [AdditiveNaturalTransformation(
            F,
            G,
            indecs,
            [x ∈ g ? popfirst!(T) : zero_morphism(F(x),G(x)) for x ∈ indecs]
        ) for T ∈ solutions]
        nats = [nats; new_nats]
    end

    return nats
end

function group_indecomposables(indecs::Vector{T}) where T <: Object
    
    G = graph_from_adjacency_matrix(Directed, [int_dim(Hom(x,y)) > 0 for x ∈ indecs, y ∈ indecs])

    groups = weakly_connected_components(G)

    return [indecs[g] for g ∈ groups]
end


#=----------------------------------------------------------
    Monoidal structures     
----------------------------------------------------------=#

function lax_monoidal_structures(F::AbstractFunctor)
    @assert is_additive(F)

    C = domain(F)
    D = codomain(F)

    S = indecomposables(C)

    F1 = Functor(
        C×C,
        D,
        X -> F(X[1]) ⊗ F(X[2]),
        f -> F(f[1]) ⊗ F(f[2])
    )

    F2 = Functor(
        C × C,
        C,
        X -> F(X[1] ⊗ X[2]),
        f -> F(f[1] ⊗ f[2])
    )

    indecs = [ProductObject(s,t) for s ∈ S, t ∈ S][:]
    nats = additive_natural_transformations(F1,F2,indecs)
        
end
#=----------------------------------------------------------
    Pretty Printing 
----------------------------------------------------------=#

function show(io::IO, η::NaturalTransformation)
    print(io, "Natural transformation")
end





