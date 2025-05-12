#=----------------------------------------------------------
    Serialize The Datatype CenterCategory. 
----------------------------------------------------------=#

@register_serialization_type CenterCategory

function save_object(s::SerializerState, C::CenterCategory)

    save_data_dict(s) do 

        save_typed_object(s, base_ring(C), :base_ring)

        save_object(s, C.category, :category)

        save_data_array(s, :simples) do 
            for x ∈ simples(C)
                save_object(s, (object(x), half_braiding(x)))
            end        
        end
    end

end
