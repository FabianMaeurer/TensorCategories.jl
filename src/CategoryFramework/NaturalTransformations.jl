#=----------------------------------------------------------
    Compute Natrural Transformations in finitary categories. 
----------------------------------------------------------=#


struct NaturalTransformation <: Morphism
    domain::AbstractFunctor 
    codomain::AbstractFunctor 
    maps::Dict{<:Object, <:Morphism} 
end

function NaturalTransformation(F::AbstractFunctor, G::AbstractFunctor, maps::Vector{<:Pair})
    NaturalTransformation(F,G,Dict(maps))
end

struct NaturalTransformations <: AbstractHomSpace
    X::AbstractFunctor
    Y::AbstractFunctor
    basis::Vector{NaturalTransformation}
    parent::VectorSpaces
end

function Nat(F::AbstractFunctor, G::AbstractFunctor)
    @assert domain(F) == domain(G)
    @assert codomain(F) == codomain(G)
    
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

        nats = [nats; [NaturalTransformation(F,G,[g[i] => T[i] for i ∈ 1:n]) for T ∈ solutions]]
    end

    return NaturalTransformations(F,G,nats, VectorSpaces(K))
end

function group_indecomposables(indecs::Vector{T}) where T <: Object
    
    graph = SimpleDiGraph([int_dim(Hom(x,y)) > 0 for x ∈ indecs, y ∈ indecs])

    groups = weakly_connected_components(graph)

    return [indecs[g] for g ∈ groups]
end



        
#=----------------------------------------------------------
    Pretty Printing 
----------------------------------------------------------=#

function show(io::IO, η::NaturalTransformation)
    show(io, """Natural transformation""")
end





