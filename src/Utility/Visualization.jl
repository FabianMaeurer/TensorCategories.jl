#=----------------------------------------------------------
    Some functionality to visualize F-symbols
----------------------------------------------------------=#

function complex_matrix_to_hsv(M::MatElem)
    complex_matrix_to_hsv(collect(M))
end

function complex_matrix_to_hsv(M::Matrix{T}) where T
    complex_matrix_to_hsv(ComplexF64.(M))
end

function complex_matrix_to_hsv(M::Matrix{ComplexF64})
    [complex_float_to_hsv(z) for z ∈ M]
end

function complex_float_to_hsv(z::ComplexF64)
    r = abs(z)
    θ = angle(z)*180/pi
    HSV(θ, sqrt(r), 1)
end

function morphism_to_image(f::Morphism)
    complex_matrix_to_hsv(matrix(f))
end

function finite_field_matrix_to_hsv(M::MatElem)
    n = Int(order(base_ring(M)))
    elems = collect(base_ring(M))
    color_map = Dict([x => HSV(360*i/n,1,1) for (i,x) ∈ enumerate(elems)])
    [color_map[x] for x ∈ M]
end