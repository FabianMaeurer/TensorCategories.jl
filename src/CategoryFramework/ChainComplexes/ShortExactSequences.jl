#=----------------------------------------------------------
    Structures for the category of short exact sequences
    in an abelian category. 
----------------------------------------------------------=#

@attributes mutable struct ShortExactSequences <: AbstractChainComplexCategory 
    category::Category

    ShortExactSequences(C::Category) = new(C)
end
   
struct ShortExactSequence <: AbstractChainComplex
    parent::ShortExactSequences
    X::Object
    Y::Object
    Z::Object
    mono::Morphism
    epi::Morphism
end

struct ShortExactSequenceMorphism <: AbstractChainComplexMorphism
    domain::ShortExactSequence
    codomain::ShortExactSequence
    left::Morphism
    middle::Morphism
    right::Morphism
end

objects(X::ShortExactSequence) = [X.X, X.Y, X.Z]
morphisms(X::ShortExactSequence) = [X.mono, X.epi]
morphisms(f::ShortExactSequenceMorphism) = [f.left, f.middle, f.right]

function ShortExactSequence(f::Morphism, g::Morphism)
    C = parent(f)
    @assert kernel(f)[1] == zero(C) 
    @assert cokernel(g)[1] == zero(C)
    @assert image(f) == kernel(g)

    par = ShortExactSequences(C)
    X = domain(f)
    Y = codomain(f)
    Z = codomain(g)
    ShortExactSequence(par, X,Y,Z, f,g)
end

morphism_type(C::ShortExactSequences) = ShortExactSequenceMorphism
#=----------------------------------------------------------
    Additive Structures 
----------------------------------------------------------=#

function direct_sum(S::ShortExactSequence...)
    C = parent(S[1])

    objs = objects.(S)
    mors = morphisms.(S)

    sum_X, i_X, p_X = direct_sum([o[1] for o ∈ objs]...)
    sum_Y, i_Y, p_Y = direct_sum([o[2] for o ∈ objs]...)
    sum_Z, i_Z, p_Z = direct_sum([o[3] for o ∈ objs]...)

    sum_mono = direct_sum([m[1] for m ∈ mors]...)
    sum_epi  = direct_sum([m[2] for m ∈ mors]...)

    sum_sequence = ShortExactSequence(C, sum_X, sum_Y, sum_Z, sum_mono, sum_epi)

    sum_incl = [Morphism(s, sum_sequence, i, j, k) for (s,i,j,k) ∈ zip(S,i_X,i_Y,i_Z)]

    sum_proj = [Morphism(sum_sequence, s, i, j, k) for (s,i,j,k) ∈ zip(S,p_X,p_Y,p_Z)]

    return sum_sequence, sum_incl, sum_proj
end

function direct_sum(f::ShortExactSequenceMorphism...)
    domains = domain.(f)
    codomain = codomain.(f)

    sum_domain = direct_sum(domains...)[1]
    sum_codomain = direct_sum(codomains...)[1]

    morphism_sum = [direct_sum([morphisms(fi)[j] for fi ∈ f]...) for j ∈ 1:3]

    return Morphism(sum_domain, sum_codomain, morphism_sum...)
end

function +(f::ShortExactSequenceMorphism, g::ShortExactSequenceMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g)
    Morphism(domain(f), codomain(f), (morphism(f) .+ morphisms(g))...)
end

function *(λ, f::ShortExactSequenceMorphism) 
    Morphism(domain(f), codomain(f), (λ .* morphisms(f))...)
end

#=----------------------------------------------------------
    Monoidal Structures 
----------------------------------------------------------=#

function tensor_product(S::ShortExactSequence, T::ShortExactSequence)
    X,Y,Z = [x⊗y for (x,y) ∈ zip(objects(S), objects(Y))]
    mono,epi = [f⊗g for (f,g) ∈ zip(morphism(f), morphisms(g))]
    ShortExactSequence(parent(X), X,Y,Z, mono, epi)
end

function tensor_product(f::ShortExactSequenceMorphism, g::ShortExactSequenceMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)
    mors = [h⊗l for (h,l) ∈ zip(morphisms(f), morphisms(g))]
    Morphism(dom,cod, mors...)
end


#=----------------------------------------------------------
    Morphisms
----------------------------------------------------------=#

function compose(f::ShortExactSequenceMorphism, g::ShortExactSequenceMorphism)
    Morphism(domain(f), codomain(g), [compose(h,l) for (h,l) ∈ zip(morphisms(f), morphisms(g))]...)
end

function Hom(S::ShortExactSequence, T::ShortExactSequence)
    base = ShortExactSequenceMorphism[]

    S_X, S_Y, S_Z = objects(S)
    T_X, T_Y, T_Z = objects(T)

    s1,s2 = morphisms(S)
    t1,t2 = morphisms(T)

    H_X = Hom(S_X, T_X)
    H_Y = Hom(S_Y, T_Y)
    H_Z = Hom(S_Z, T_Z)

    if int_dim.([H_X,H_Y,H_Z]) == 0 
        return HomSpace(S,T,base)
    end

    base_X = [zero_morphism(S_X,T_X); basis(H_X)]
    base_Y = basis(H_Y)
    base_Z = [zero_morphism(S_Z,T_Z); basis(H_Z)]

    base_1 = basis(Hom(S_X, T_Y))
    base_2 = basis(Hom(S_Y, T_Z))

    F = base_ring(S)

    mat_1 = hcat([express_in_basis(h∘s1, base_1) for h ∈ base_Y]...)
    mat_2 = hcat([express_in_basis(s2∘h, base_2) for h ∈ base_Y]...)

    M = matrix(F, length(base_Y), length(base_1) + length(base_2), [mat_1; mat_2])

    for f ∈ H_X, g ∈ H_Z
        N_1 = matrix(F, length(base_1), 1, express_in_basis(t1∘f, base_1))
        N_2 = matrix(F, length(base_2), 1, express_in_basis(g∘s2, base_2))
        N = [N_1; N_2]

        sols = solve_left(transpose(M),transpose(N))

        n,m = size(sols)

        if n*m == 0 
            if f == 0 && g == 0 
                continue
            else
                push!(base, ShortExactSequenceMorphism(S,T, f, zero_morphism(S_Y,T_Y), g))
            end
        end

        B = collect(eachrow(collect(sols)))

        base = [base; [ShortExactSequenceMorphism(S, T, f, sum(b .* base_Y), g) for b ∈ B]]
    end

    return HomSpace(S, T, base)
end



        