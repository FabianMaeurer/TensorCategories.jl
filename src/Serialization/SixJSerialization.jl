#=----------------------------------------------------------
    Serialize Categories of Type SixJCategory 
----------------------------------------------------------=#

const database_path = joinpath(@__DIR__,"src/SixJCategoryDatabase/")

@register_serialization_type SixJCategory

function save_object(s::SerializerState, C::SixJCategory)
    save_data_dict(s) do 

        save_typed_object(s, base_ring(C), :base_ring)

        save_object(s, C.simples, :simples)

        save_data_array(s, :simples_names) do   
            for n ∈ simples_names(C) 
                save_object(s,n)
            end
        end

        save_data_array(s, :tensor_product) do 
            for x ∈ C.tensor_product 
                save_object(s,x)
            end
        end

        save_data_array(s, :ass) do 
            for m ∈ C.ass[:]
                save_object(s, m)
            end
        end

        
        #save_typed_object(s,C.spherical, :spherical)

        try 
            save_object(s, C.one, :one)
        catch
        end

        try 
            save_object(s, C.name, :name)
        catch 
        end

        if is_braided(C)
            save_data_array(s, :braiding) do 
                for m ∈ C.braiding 
                    save_object(s,m)
                end
            end
        end
    end
end

function load_object(s::DeserializerState, ::Type{SixJCategory})
    C = SixJCategory()
    
    C.base_ring = load_typed_object(s, :base_ring)
    
    C.simples = load_object(s, Int64, :simples)
    
    C.simples_names = load_array_node(s, :simples_names) do (i,n)
        load_object(s, String)
    end

    n = C.simples

    C.tensor_product = reshape(
        load_array_node(s, :tensor_product) do (i,a) 
            load_object(s,Int64)
        end,
        n,n,n
    )

    C.ass = reshape(
        load_array_node(s, :ass) do (i,m)
            m = load_object(s, Matrix, (elem_type(base_ring(C)), base_ring(C)))
            matrix(base_ring(C), size(m,1),size(m,2), m)
        end,
        n,n,n,n
    )


    try
        C.one = load_object(s, Vector{Int64}, :one)
    catch
    end

    try
        C.name = load_object(s, String, :name)
    catch 
    end

    try
        C.braiding = reshape(
            load_array_node(s, :braiding) do (i,m)
                m = load_object(s, Matrix, (elem_type(base_ring(C)), base_ring(C)))
                matrix(base_ring(C), size(m,1),size(m,2), m)
            end,
            n,n,n
        )
    catch 
    end

    try 
        set_canonical_spherical!(C)
    catch 
    end
    
    return C
end