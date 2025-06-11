#=----------------------------------------------------------
    Check the properties of AnyonWiki categories
----------------------------------------------------------=#




function test_cyclotomic_anyonwiki(properties = nothing)
    # For every cyclotomic category check if the properties
    # align with the provided ones

    # Return a list with misaligned properties of the form 
    # n => [property => [expected, actual]]
    errors = Dict{Int, Dict{Symbol, Vector{Bool}}}()

    cyclic_n = parse.(Int, readlines(joinpath(@__DIR__, "AnyonWikiData/cyclopos.txt")))

    for n ∈ cyclic_n

        # Check properties 
        all_passed, n_errors = test_anyonwiki_cyclotomic(n, properties)

        if !all_passed 
            push!(errors, n => n_errors)
        end
    end

    return sort(errors)
end


function test_cyclotomic_anyonwiki(n::Int, properties = nothing)
    # The attributes as in the file as a Dict
    # (Braided, Unitary, Spherical, Ribbon, Modular) 
    expected_properties = TensorCategories.load_anyonwiki_attributes(n)
    push!(expected_properties, :monoidal => true)
    push!(expected_properties, :pivotal => true)


    C = TensorCategories.anyonwiki(n)

    # Get actual properties
    actual_properties = anyonwiki_checks(n)

    errors = Dict{Symbol, Vector{Bool}}()

    properties = properties === nothing ? [:monoidal, :pivotal, :spherical, :braided, :modular] : properties

    # Compare and collect errors 
    for prop ∈ properties

        if (ex = expected_properties[prop]) != (ac = actual_properties[prop])
            # Push error data
            push!(errors, prop => [ex, ac])
        end
    end
    return isempty(errors), errors
end


function anyonwiki_checks(n::Int)
    # For now we only check pentagon, hexagon,
    # pivotal and spherical 

    # Construct the nth list from Wiki
    C = TensorCategories.anyonwiki(n)

    # the result will be a Dict property/structure => bool
    result = Dict(
        :monoidal   => false,
        :braided    => false,
        :pivotal    => false,
        :spherical  => false,
        :modular    => false
    )

    # Check the pentagon axiom
    result[:monoidal] = pentagon_axiom(C)

    # If pentagon failed quit 
    !result[:monoidal] && return result

    # Check the pivotal structure
    result[:pivotal] = is_pivotal(C) 
    
    # If the pivotal structure is correct check if it is spherical
    if result[:pivotal] 
        result[:spherical] = is_spherical(C)
    end

    # if there is a braiding provided we check that
    if isdefined(C, :braiding)
        
        # Check hexagon axiom
        result[:braided] = hexagon_axiom(C) 
            
        # If the hexagon axiom fails quit
        ! result[:braided] && return result

        # If the braiding is correct and C is spherical check for madularity
        if result[:spherical] 
            result[:modular] = is_modular(C)
        end
    end

    return result
end