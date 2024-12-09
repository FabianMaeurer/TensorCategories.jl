#=----------------------------------------------------------
    Methods for decomposition of objects in abelian 
    categories. 
    
    In the cases cases where it is possible we provide 
    methods to compute simple subobjects. 
----------------------------------------------------------=#

function simple_subobjects(X::Object, E = End(X), is_simple = false)
    #=  Compute all simple subobjects in a tensor category

        The approach is the MeatAxe algorithm. 
    =#


    # Over QQBar it's easier
    K = base_ring(X)
    if K == QQBar || typeof(K) == CalciumField || typeof(K) <: Union{PadicField,QadicField}
        return simple_subobjects_over_qqbar(X,E)
    end

    # Just a single Endomorphism -> Done
    if length(basis(E)) == 1
        return [X]
    end

    #A simple algebra of squarefree dimension is a division algebra
    if is_simple && is_squarefree(dim(E)) 
        return [X]
    end

    R = endomorphism_ring(X,E)

    i = findfirst(f -> !is_irreducible(f), minpoly.(basis(E)))

    if !is_simple && is_semisimple(endomorphism_ring(X,E))
        img = [image(i)[1] for i ∈ central_primitive_idempotents(E)]
        if length(img) == int_dim(E)
            return img
        end
        return vcat([simple_subobjects(i, End(i), true) for i in img]...)
    end

    if i !== nothing && is_semisimple(R)
        f = E[i]
        l = roots(minpoly(f))[1]
        K = kernel(f - l*id(X))[1]
        I = image(f - l*id(X))[1]
        K = simple_subobjects(K, End(K), true)
        I = simple_subobjects(I, End(I), true)

        return unique_simples([K;I])
    end


    M = regular_module(R)
    M.action_of_gens = [representation_matrix(g) for g ∈ gens(R)]
  
    submods = minimal_submodules(M)
    fi = [sum(collect(s)[1,:] .* basis(E)) for s in submods]

    unique_simples([image(f)[1] for f ∈ fi])
    # indecomposables = indecomposable_subobjects(X, E)
    # if is_semisimple(parent(X))
    #     return indecomposables
    # else
    #     return indecomposables[is_simple.(indecomposables)]
    # end
end
