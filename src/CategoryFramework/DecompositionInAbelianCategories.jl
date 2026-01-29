#=----------------------------------------------------------
    Methods for decomposition of objects in abelian 
    categories. 
    
    In the cases cases where it is possible we provide 
    methods to compute simple subobjects. 
----------------------------------------------------------=#

function simple_subobjects(X::Object, E = End(X), _is_simple = false, is_indecomposable = false)
    #=  Compute all simple subobjects in a tensor category

        The approach is the MeatAxe algorithm. 
    =#
    # Over QQBar it's easier
    K = base_ring(X)
    if K == QQBarField() || typeof(K) <: Union{CalciumField, PadicField,QadicField, AcbField}
        return simple_subobjects_over_qqbar(X,E)
    end

    # Just a single Endomorphism -> Done
    if length(basis(E)) == 1
        return [X]
    end

    R = endomorphism_ring(X,E)

    #A simple algebra of squarefree dimension is a division algebra

    if _is_simple && is_squarefree(div(dim(E),dim(center(R)[1]))) 
        return [X]
    end

    if !is_indecomposable #&& is_semisimple(endomorphism_ring(X,E))
        img = [image(i)[1] for i ∈ central_primitive_idempotents(E)]
        if length(img) == int_dim(E)
            return img
        end
        return vcat([simple_subobjects(i, End(i), _is_simple, true) for i in img]...)
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



    if is_squarefree(div(dim(E),dim(center(R)[1]))) && is_simple(R)
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


function simple_subobjects_of_free_object(X::Object, _free_construction::Function, _free_adjunction::Function)

    free_X = _free_construction(X)

    H = _free_adjunction(X, free_X)

    idems = central_primitive_idempotents(H)

    S = typeof(free_X)[]

    # compute the direct sum components and their projections
    components = [image(e) for e in idems]

    # For each component compute the endomorphism space 
    # which is a simple matrix algebra. Find the corresponding
    # simple.
    for (Y, incl) ∈ components
        E = end_of_free_subobject(X, incl, _free_adjunction, H)
        
        append!(S, simple_subobjects(Y, E, false, true))
    end
    S
end

@doc raw""" 

    end_of_free_subobject(X::Object, incl::Morphism, free_adjunction::function, E::HomSpace)

Compute the endomorphism space ``End(domain(incl))`` via ``E = End(free(X))``. 
"""
function end_of_free_subobject(X::Object, incl::Morphism, _free_adjunction::Function, E::HomSpace = _free_adjunction(X,codomain(incl)))

    free_X = codomain(incl)
    M = domain(incl)

    # Compute a projection onto the subobject
    H = _free_adjunction(X, M)
    i = findfirst(f -> cokernel(f)[1] == 0, basis(H))
    
    proj = if i === nothing 
        while i === nothing
            c = rand(-1:1, length(basis(H)))
            global f = sum(c[j]*basis(H)[j] for j ∈ 1:length(basis(H)))
            i = is_zero(cokernel(f)[1]) ? 1 : nothing
        end
        f
    else
        basis(H)[i]
    end

    # collect morphisms in End(M)
    mors = [proj ∘ f ∘ incl for f ∈ basis(E)]

    return HomSpace(M,M, basis(mors))
end