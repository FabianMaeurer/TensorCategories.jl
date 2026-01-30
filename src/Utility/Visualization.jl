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