#=----------------------------------------------------------
    Set Up Serialization 
----------------------------------------------------------=#

using Serialization
import Oscar: @register_serialization_type,
                save_data_dict,
                save_data_array,
                save_object,
                save_type_params,
                save_typed_object,
                load_object,
                load_type_params,
                load_typed_object,
                SerializerState,
                DeserializerState,
                encode_type,
                register_serialization_type,
                load_array_node,
                serialize_with_params