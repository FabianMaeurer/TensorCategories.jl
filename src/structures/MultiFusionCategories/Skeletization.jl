function SixJCategory(C::Category, names::Vector{String} = ["X$i" for i ∈ 1:length(simples(C))])
    @assert is_multifusion(C)

    S = simples(C)
    n = length(S)
    F = base_ring(C)

    mult = Array{Int,3}(undef,n,n,n)

    for i ∈ 1:n, j ∈ 1:n
        mult[i,j,:] = [int_dim(Hom(S[k],S[i]⊗S[j])) for k ∈ 1:n]
    end

    # Define SixJCategory
    skel_C = SixJCategory(F,names)

    set_tensor_product!(skel_C,mult)

    set_associator!(skel_C, six_j_symbols(C, S))

    # Try to set spherical
    try 
        set_spherical!(skel_C,[F(spherical(s)) for s ∈ S])
    catch
    end

    set_one!(skel_C, [int_dim(Hom(one(C), s)) for s ∈ S])

    set_name!(skel_C, "Skeletization of $C")

    return skel_C, SkeletizationFunctor(C, skel_C)
end

function six_j_symbols(C::Category, S = simples(C))
    @assert is_multifusion(C)

    n = length(S)
    ass = Array{MatElem,4}(undef,n,n,n,n)

    objects = []
    isos = []

    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        X = (S[i] ⊗ S[j]) ⊗ S[k]
        Y = S[i] ⊗ (S[j] ⊗ S[k])

        _, ij_decomposition_morphism, ij_incl, ij_proj = direct_sum_decomposition(S[i] ⊗ S[j], S)
        _, jk_decomposition_morphism, jk_incl, jk_proj =  direct_sum_decomposition(S[j] ⊗ S[k], S)

        distribute_before = (inv(ij_decomposition_morphism) ⊗ id(S[k])) ∘ inv(distribute_left([domain(i) for i ∈ ij_incl], S[k]))

        distribute_after = distribute_right(S[i], [domain(i) for i ∈ jk_incl]) ∘ (id(S[i]) ⊗ jk_decomposition_morphism)

        # @show codomain(distribute_after) == domain(distribute_before)

        summands = vcat([[x for _ ∈ 1:k] for (x,k) ∈ decompose(X, S)]...)

        Z, before, incl, _ = direct_sum_decomposition(X,S)
        Z2, after, _, proj = direct_sum_decomposition(Y,S)
        
        before = inv(before) ∘ distribute_before ∘ before
        after = after ∘ distribute_after ∘ inv(after)
        # if X ∈ objects
        #     before = isos[findfirst(e -> e == X, objects)]
        # else
        #     _,before = is_isomorphic(Z,X)
        #     isos = [isos; before]
        #     objects = [objects; X]
        # end

        # if Y ∈ objects
        #     after = isos[findfirst(e -> e == Y, objects)]
        # else
        #     _,after = is_isomorphic(Z2,Y)
        #     isos = [isos; after]
        #     objects = [objects; Y]
        # end

    
        ass_mor =  after ∘ associator(S[i],S[j],S[k]) ∘ before 

        for l ∈ 1:n
            before_incl_l = filter(f -> domain(f) == S[l], incl)
            after_proj_l = filter(f -> codomain(f) == S[l], proj)
            
            if length(before_incl_l) == 0
                ass[i,j,k,l] = zero_matrix(base_ring(C),0,0)
                continue
            end

            m = length(before_incl_l)
            F = base_ring(C)
            ass[i,j,k,l] = zero_matrix(F,m,m)
            for q ∈ 1:m, p ∈ 1:m
                ass[i,j,k,l][p,q] = F(after_proj_l[q] ∘ ass_mor ∘ before_incl_l[p])
            end
    
        end
    end
    return ass
end

struct SkeletizationFunctor <: AbstractFunctor
    domain::Category
    codomain::Category
end

function SkeletizationFunctor(C::Category)
    SkeletizationFunctor(C,SixJCategory(C))
end

function (F::SkeletizationFunctor)(X::Object)
    return SixJObject(codomain(F), [int_dim(Hom(X,s)) for s ∈ simples(C)])
end

function (F::SkeletizationFunctor)(f::Morphism)
    X = domain(F)
    Y = codomain(F)

    S = simples(domain(f))
    n = length(S)

    _, before, before_incl, before_proj = direct_sum_decomposition(X, S)
    _, after,  after_incl,  after_proj  = direct_sum_decomposition(Y, S)

    m = Vector{MatElem}(undef,n)

    for l ∈ 1:n
        before_incl_l = filter(f -> domain(f) == S[l], before_incl)
        before_proj_l = filter(f -> codomain(f) == S[l], before_proj)
        after_incl_l = filter(f -> domain(f) == S[l], after_incl)
        after_proj_l = filter(f -> codomain(f) == S[l], after_proj)
        
        if length(before_incl_l) == 0
            m[l] = zero_matrix(base_ring(C),0,0)
            continue
        end
        m[l] = matrix(direct_sum(after_proj_l) ∘ ass_mor ∘ direct_sum(before_incl_l))
    end

    return Morphism(F(X),F(Y),m)
end

function show(io::IO, F::SkeletizationFunctor)
    print(io,"Skeletal projection of $(domain(F))")
end


