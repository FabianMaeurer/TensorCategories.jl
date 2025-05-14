#=----------------------------------------------------------
    Methods for decomposition of objects in abelian 
    categories. 
    
    In the cases cases where it is possible we provide 
    methods to compute simple subobjects. 
----------------------------------------------------------=#

function simple_subobjects(X::Object, E = End(X), is_simple = false, is_indecomposable = false)
    #=  Compute all simple subobjects in a tensor category

        The approach is the MeatAxe algorithm. 
    =#
    # Over QQBar it's easier
    K = base_ring(X)
    if K == QQBarField() || typeof(K) == CalciumField || typeof(K) <: Union{PadicField,QadicField}
        return simple_subobjects_over_qqbar(X,E)
    end

    # Just a single Endomorphism -> Done
    if length(basis(E)) == 1
        return [X]
    end

    R = endomorphism_ring(X,E)

    #A simple algebra of squarefree dimension is a division algebra

    if is_simple && is_squarefree(div(dim(E),dim(center(R)[1]))) 
        return [X]
    end

    i = findfirst(f -> degree(f) > 1 && length(roots(f)) > 0, minpoly.(basis(E)))

    if i !== nothing #&& is_semisimple(R)
        f = E[i]
        min_f = minpoly(f)
        rs = roots(min_f)
        is_irreducible(min_f)
        if length(rs) == int_dim(E) && is_semisimple(parent(X))
            return unique_simples([kernel(f - r*id(X))[1] for r ∈ rs])
        end

        l = rs[1]
        K = kernel(f - l*id(X))[1]
        I = image(f - l*id(X))[1]
        K = simple_subobjects(K, End(K), false)
        I = simple_subobjects(I, End(I), false)

        return unique_simples([K;I])
    end

   
    if !is_indecomposable #&& is_semisimple(endomorphism_ring(X,E))
        img = [image(i)[1] for i ∈ central_primitive_idempotents(E)]
        if length(img) == int_dim(E)
            return img
        end
        return vcat([simple_subobjects(i, End(i), is_simple, true) for i in img]...)
    end

 

    if is_squarefree(div(dim(E),dim(center(R)[1]))) && issimple(R)
        return [X]
    end


    if schur_index_over_center(R)^2 == dim(E)//dim(center(R)[1])
        return [X]
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
