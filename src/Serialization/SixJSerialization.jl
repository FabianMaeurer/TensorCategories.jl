#=----------------------------------------------------------
    Serialize Categories of Type SixJCategory 
----------------------------------------------------------=#

const database_path = joinpath(@__DIR__,"src/SixJCategoryDatabase/")

function save_object(s::SerializerState, C::SixJCategory)
    save_data_dict(s) do 

        save_typed_object(s, base_ring(C), :base_ring)

        save_object(s, C.simples, :simples)

        save_object(s, Tuple(simples_names(C)), :simples_names)
        # save_data_array(s, :simples_names) do   
        #     for n ∈ simples_names(C) 
        #         save_object(s,n)
        #     end
        # end

        save_object(s, Tuple(C.tensor_product), :tensor_product)
        # save_data_array(s, :tensor_product) do 
        #     for x ∈ C.tensor_product 
        #         save_object(s,x)
        #     end
        # end

        #save_typed_object(s, Dict(Tuple(k) => v for (k,v) ∈ F_symbols(C)), :F_symbols)
        save_data_array(s, :ass) do 
            for m ∈ C.ass[:]
                save_object(s, m)
            end
        end

        
        #save_typed_object(s,C.spherical, :spherical)

        if isdefined(C, :pivotal)
            save_typed_object(s, Tuple(C.pivotal), :pivotal)
        end

        if isdefined(C, :embedding)
            save_typed_object(s, getfield(C, :embedding), :embedding)
        end

        if isdefined(C, :one)
            save_object(s, Tuple(C.one), :one)
        end

        if isdefined(C, :name)
            save_object(s, C.name, :name)
        end

        if is_braided(C)
            save_data_array(s, :braiding) do 
                for m ∈ C.braiding[:]
                    save_object(s, m)
                end
            end
        end
    end
end

function load_object(s::DeserializerState, ::Type{SixJCategory})
    C = SixJCategory()
    
    C.base_ring = load_typed_object(s, :base_ring)
    
    C.simples = load_object(s, Int64, :simples)
    
    n = C.simples

    C.simples_names = collect(load_object(s, NTuple{n, String}, :simples_names))

    # C.simples_names = load_array_node(s, :simples_names) do (i,n)
    #     load_object(s, String)
    # end



    C.tensor_product = reshape(
        collect(load_object(s, NTuple{n^3 ,Int}, :tensor_product)),
        n,n,n
    )

    m = maximum(C.tensor_product)
    _n = m == 1 ? 6 : 10

    # F_symb = Dict(collect(k) => v for (k,v) ∈ load_typed_object(s, :F_symbols))

    # set_attribute!(C, :F_symbols, F_symb)
    # ass = dict_to_associator(n, base_ring(C), F_symb)

    # C.ass = ass
    C.ass = reshape(
        load_array_node(s, :ass) do (i,m)
            m = load_object(s, Matrix{elem_type(base_ring(C))}, base_ring(C))
            a = size(m,1)
            b = length(size(m)) == 2 ? size(m,2) : 0
            matrix(base_ring(C), a, b, m)
        end,
        n,n,n,n
    )


    if haskey(s, :pivotal)
       C.pivotal = collect(load_typed_object(s, :pivotal))
    end

    if haskey(s, :one)
        C.one = collect(load_object(s, NTuple{n,Int64}, :one))
    end

    if haskey(s, :embedding)
        C.embedding = load_typed_object(s, :embedding)
    end

    if haskey(s, :name)
        C.name = load_object(s, String, :name)
    end

    if haskey(s, :braiding)
        C.braiding = reshape(
            load_array_node(s, :braiding) do (i,m)
                m = load_object(s, Matrix{elem_type(base_ring(C))}, base_ring(C))
                a = size(m,1)
                b = length(size(m)) == 2 ? size(m,2) : 0
                matrix(base_ring(C), a, b, m)
            end,
        n,n,n
    )
    end

    
    return C
end

#=----------------------------------------------------------
    Serialize SixJMorphism 
----------------------------------------------------------=#

type_params(X::SixJObject) = TypeParams(X.parent, parent(X))

function save_object(s::SerializerState, X::SixJObject)
    save_data_dict(s) do 
        save_object(s, X.components, :components)
    end
end

function load_object(s::DeserializerState, ::Type{SixJObject},params::SixJCategory)
    components = load_object(s, Vector{Int64}, :components)
    return SixJObject(parent, components)
end

function save_object(s::SerializerState, m::SixJMorphism)
    save_data_dict(s) do 
        save_object(s, m.domain, :domain)
        save_object(s, m.codomain, :codomain)
        save_object(s, m.m, :mats)
    end
end

function load_object(s::DeserializerState, ::Type{SixJMorphism}, parent::SixJCategory)
    domain = load_object(s, SixJobject, :domain, parent)
    codomain = load_object(s, SixJObject, :codomain, parent)
    m = load_object(s, Array{MatElem}, :mats)
    return SixJMorphism(domain, codomain, m)
end