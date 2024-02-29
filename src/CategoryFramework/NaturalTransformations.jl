#=----------------------------------------------------------
    Compute Natrural Transformations in finitary categories. 
----------------------------------------------------------=#
abstract type NaturalTransformation <: Morphism end

struct AdditiveNaturalTransformation <: NaturalTransformation
    domain::AbstractFunctor 
    codomain::AbstractFunctor 
    maps::Vector{<:Morphism} 
end

function AdditiveNaturalTransformation(F::AbstractFunctor, G::AbstractFunctor, indecs::Vector{<:Object}, maps::Vector{<:Pair})

    maps = [x ∈ keys(maps) ? maps[x] : zero_morphism(F(x),G(x)) for x ∈ indecs]
    AdditiveNaturalTransformation(F,G,maps)
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

function compose(η::NaturalTransformation, σ::NaturalTransformation)

end

#=----------------------------------------------------------
    Compute natural transformations 
----------------------------------------------------------=#

function Nat(F::AbstractFunctor, G::AbstractFunctor)
    @assert domain(F) == domain(G)
    @assert codomain(F) == codomain(G)

    if is_additive(F) && is_additive(G)
        return additive_natural_transformations(F,G)
    end

    @error("Cannot compute natural transformations")
end

function additive_natural_transformations(F::AbstractFunctor, G::AbstractFunctor)
    
    C = domain(F)
    K = base_ring(C)
    indecs = indecomposables(C)

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

        Kx,x = PolynomialRing(K,sum(ns))

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

        nats = [nats; [AdditiveNaturalTransformation(F,G,indecs, [g[i] => T[i] for i ∈ 1:n]) for T ∈ solutions]]
    end

    return NaturalTransformations(F,G,nats, VectorSpaces(K))
end

function group_indecomposables(indecs::Vector{T}) where T <: Object
    
    graph = SimpleDiGraph([int_dim(Hom(x,y)) > 0 for x ∈ indecs, y ∈ indecs])

    groups = weakly_connected_components(graph)

    return [indecs[g] for g ∈ groups]
end


#=----------------------------------------------------------
    Monoidal structures     
----------------------------------------------------------=#

function lax_monoidal_structures(F::AbstractFunctor)
    @assert is_additive(F)

    C = domain(F)

    T = Functor(
        C×C,

    )
        
end
#=----------------------------------------------------------
    Pretty Printing 
----------------------------------------------------------=#

function show(io::IO, η::NaturalTransformation)
    print(io, "Natural transformation")
end





