function kernel(Mat::AcbMatrix; side = :left)
    if side == :right 
        return transpose(kernel(transpose(Mat)))
    else 
        M = transpose(Mat) 
    end
    F = base_ring(M)

    if overlaps(Mat, zero(parent(Mat)))
        return identity_matrix(F, size(Mat,1))
    end

    if size(Mat,2) == 0 
        return identity_matrix(F,size(M,2))
    elseif size(Mat,1) == 0
        return zero_matrix(F,0,0)
    end

    acc = minimum([Int(floor(minimum([a for a in Oscar.accuracy_bits.(M) if a > 0], init = precision(F)))), Int(floor(precision(F)))])

    ker = Vector{AcbFieldElem}[] 

    done = false 
    _M = [overlaps(m, F(0)) ? F(0) : m for m in collect(M)]
    taken = []
    while !done 
        coeffs = complex_lindep(_M, acc)
        n = size(M,1)
        i = findfirst(!=(0), coeffs)

        if i === nothing 
            return zero_matrix(F,0,size(M,2))
        end

        # coeffs[i]
        # coeffs = 1//coeffs[i] .* coeffs 
        
        if sum(coeffs .!= 0, init = 0) <= 1 && !all(overlaps(k,F(0)) for k ∈ _M[:,i])
            done = true 
        else 
            taken = [taken; findfirst(!iszero, coeffs)]
            for j in reverse(taken[1:end-1])
                insert!(coeffs, j, F(0))
            end

            ker = [ker; [coeffs]]
            _M = collect(_M[:,setdiff(1:size(_M,2), taken[end])])
            isempty(_M) && break
        end
    end
    
    ker = length(ker) == 0 ? Matrix{AcbFieldElem}(undef,size(M,2),0) : hcat(ker...)

    ker = [overlaps(k, F(0)) ? F(0) : k for k in ker]

    ker = transpose(matrix(F,size(M,2), size(ker,2),ker))

    if size(ker,1) == 0 
        return ker
    end

    if overlaps(ker * Mat, zero_matrix(F, size(ker,1),size(Mat,2)))
        return ker 
    end


    # Add the error to make the kernel satisfy ker*M = 0
    err_mat = ker * Mat 

    real_err = maximum(abs.((err_mat)))
    imag_err = maximum(abs.((err_mat)))

    # add the errors 
    real_ker = real.(ker)
    imag_ker = imag.(ker)
    for i ∈ eachindex(real_ker) 
        r = real_ker[i]
        is_exact(r) && continue
        add_error!(r, real_err)
        real_ker[i] = r
    end

    for i ∈ eachindex(imag_ker)
        c = imag_ker[i] 
        is_exact(c) && continue
        add_error!(c, imag_err)
        imag_ker[i] = c
    end
    k = matrix(F, size(ker,1), size(ker,2), [F(r,c) for (r,c) ∈ zip(real_ker, imag_ker)])
    overlaps(k * Mat, zero_matrix(F, size(ker,1),size(Mat,2)))
    side == :right && return transpose(k)
    k
end

cokernel(M::AcbMatrix, side = :left) = transpose(kernel(transpose(M), side = side))

function complex_lindep(M::Matrix{AcbFieldElem}, m::Int)
    @assert m > 0 

    if isempty(M) 
        return AcbFieldElem[]
    end

    F = parent(M[1,1])
    R = parent(real(F(1)))

    if all(overlaps(x,F(0)) for x in M)
        return AcbFieldElem[F(0) for _ ∈ size(M,2)]
    end

    try 
        coeffs = R.(Oscar.lindep([M (F(im) .*M)], m))
    catch
        return complex_lindep(M, m - 2)
    end
    n = div(length(coeffs), 2)

    #normalise 
    j = findfirst(!=(0), coeffs)
    if j !== nothing 
        coeffs = 1//coeffs[j] .* coeffs
    end

    if sum(!=(0).(coeffs[1:n] .+ F(im) .* coeffs[n+1:end])) < 2
        return coeffs[1:n] .+ (F(im) .* coeffs[n+1:end])
    end

    if all(overlaps(x,F(0)) for x ∈ [sum(c * m for (c,m) in zip(coeffs[1:n] .+ F(im) .* coeffs[n+1:end],mm)) for mm in eachrow(M)])
        return coeffs[1:n] .+ (F(im) .* coeffs[n+1:end])
    end
   
    err = maximum(abs.([sum(c * m for (c,m) in zip(coeffs[1:n] .+ F(im) .* coeffs[n+1:end],mm)) for mm in eachrow(M)]))

    if abs(err) > 10^(-8)
       coeffs = coeffs[1:n] .+ (F(im) .* coeffs[n+1:end])
       i = findfirst(!=(0), coeffs)
       return insert!(complex_lindep(M[:, setdiff(1:n, i)], m), i, F(0))
    end

    for i ∈ eachindex(coeffs)
        r = R(coeffs[i])
        r == 0 && continue
        add_error!(r, err)#R("+/- 1e-$err"))
        coeffs[i] = r 
    end

    return coeffs[1:n] .+ (F(im) .* coeffs[n+1:end])
end

function nullspace(M::AcbMatrix)
    k = kernel(M)
    (size(k,1), k)
end

function add_error(x::AcbFieldElem, err::ArbFieldElem)
    r = real(x)
    if !is_exact(r) 
        add_error!(r, err)
    end

    c = imag(x) 
    if !is_exact(c)
        add_error!(c, err)
    end

    return parent(x)(r,c)
end
