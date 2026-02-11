#=----------------------------------------------------------
    MeatAxe analog for module objects in tensor categories 
----------------------------------------------------------=#

function spin_submodule(M::RightModuleObject, u::Morphism)
    # Take a monomorphism u: U ↪ M in the parent category 
    # and compute the image of U⊗A → M. Repeat the process 
    # until the image is isomorphic to U
    # (Analog to the spinning algorithm for modules over an
    # algebra)

    A = object(algebra(parent(M)))
    
    action_on_U = compose(
        u ⊗ id(A),
        right_action(M)
    )

    image_of_U, inclusion = image(action_on_U)

    if is_isomorphic(image_of_U, domain(u))[1]
        return RightModuleObject(
            parent(M),
            image_of_U, 
            compose(
                inclusion ⊗ id(A),
                right_action(M),
                left_inverse(inclusion)
            )
        ), inclusion
    else
        return spin_submodule(M, inclusion)
    end
end

function spin_submodule(M::LeftModuleObject, u::Morphism)
    # Take a monomorphism u: U ↪ M in the parent category 
    # and compute the image of A ⊗ U → M. Repeat the process 
    # until the image is isomorphic to U
    # (Analog to the spinning algorithm for modules over an
    # algebra)
    A = object(algebra(parent(M)))
    
    action_on_U = compose(
        id(A) ⊗ u,
        left_action(M)
    )

    image_of_U, inclusion = image(action_on_U)

    if is_isomorphic(image_of_U, domain(u))[1]
        return LeftModuleObject(
            parent(M),
            image_of_U, 
            compose(
                id(A) ⊗ inclusion,
                left_action(M),
                left_inverse(inclusion)
            )
        ), inclusion 
    else
        return spin_submodule(M, inclusion)
    end
end

function meataxe(M::RightModuleObject) 
    _,_,incl,proj = direct_sum_decomposition(object(algebra(parent(M))))
    meataxe(M, incl, proj)
end

function meataxe(M::RightModuleObject, incl::Vector{T}, proj::Vector{T}) where T <: Morphism
    # Apply the MeatAxe approach to a module object
    # over an algebra object in a fusion category. 
    # Let μ: M⊗A → M be the right action on M. 
    
    # Find a random "element" a: A'↪ A such that 
    # ker(id(M)⊗a ∘ μ) ≠ 0.
    C = category(parent(M))
    A = algebra(parent(M))

    if is_simple(object(M))
        return true, M, id(M)
    end
    
    kernel_of_a = zero(C)
    kernel_inclusion = zero_morphism(zero(C),zero(C))
    # Dummy for now 
    
    n = length(incl)
    K = base_ring(C)
    θ = sum(K.(rand(Int,n)) .* [i ∘ p for (i,p) ∈ zip(incl, proj)])
    kernel_non_zero = false

    while !kernel_non_zero
        polys = [k for (k,v) ∈ collect(factor(minpoly(θ)))]
        polys = sort(polys, by = degree)
        for p ∈ polys 
            p(θ) == 0 && continue 
            ξ = p(θ)  
            @show rank(matrix(ξ))
            action_of_a = compose(
                id(object(M)) ⊗ ξ,
                right_action(M)
            )
            action_of_a = compose(
                id(object(M)) ⊗ image(ξ)[2],
                right_action(M)
            )

            kernel_of_a, kernel_inclusion = kernel(action_of_a)
            @show int_dim(Hom(kernel_of_a, object(M)))
            if !is_zero(int_dim(Hom(kernel_of_a, object(M))))
                kernel_non_zero = true
                break
            end
        end
        θ = sum(rand(Bool,n) .* [i ∘ p for (i,p) ∈ zip(incl, proj)])
    end

    # for (i,p) ∈ zip(incl,proj)
    #     action_of_a = compose(
    #         id(object(M)) ⊗ i,
    #         right_action(M)
    #     )
    #     @show kernel_of_a, kernel_inclusion = kernel(action_of_a)
    #     if !is_zero(kernel_of_a)
    #         break
    #     end
    # end

    # if is_zero(kernel_of_a)
    #     return true, M, id(M)
    # end

    for f ∈ basis(Hom(kernel_of_a, object(M)))
        N, inclusion = spin_submodule(M,image(f)[2])
        if !is_invertible(inclusion) && int_dim(End(N)) != 0
            return (false, N, morphism(N,M,inclusion))
        end 
    end

    f = Hom(dual(kernel_of_a), dual(object(M)))[1]
    N, inclusion = spin_submodule(transposed_module(M), image(f)[2]) 
    
    if !is_invertible(inclusion) && int_dim(End(N)) != 0
        N = transposed_module(N)
        return false, N, morphism(N,M,right_inverse(dual(inclusion)))
    end

    return true, M, id(M)
end




function minimal_subquotients_with_multiplicity(M::ModuleObject)

    if is_semisimple(parent(M))
        return decompose(M)
    end

    M == zero(parent(M)) && return Tuple{typeof(M), Int}[]
    if has_attribute(M, :is_simple)
        get_attribute(M, :is_simple) && return [(M,1)]
    end
    
    # Check M for simplicity
    irred, N, incl = meataxe(M)

    if irred
        set_attribute!(M, :is_simple, true)
        return [(M,1)]
    end

    # Get the quotient object Q = M/N
    Q,_ = cokernel(incl)

    # Compute subquotients of N and Q
    sub_list = minimal_subquotients_with_multiplicity(N)
    quo_list = minimal_subquotients_with_multiplicity(Q)

    # Collect non-isomorphic subquotients
    
    collected = falses(length(quo_list))
    for i ∈ 1:length(sub_list)
        for j ∈ 1:length(quo_list)
            if !collected[j] && is_isomorphic(sub_list[i][1], quo_list[j][1])[1]
                sub_list[i] = (sub_list[i][1], sub_list[i][2] + quo_list[j][2])
                collected[j] = true
                break
            end
        end
    end

    for i ∈ 1:length(quo_list)
        if !collected[i]
            push!(sub_list, quo_list[i])
        end
    end

    return sub_list

end

function minimal_subquotients(M::ModuleObject)
    [x for (x,_) ∈ minimal_subquotients_with_multiplicity(M)]
end


