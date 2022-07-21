function duplicate_row(M::MatElem, i::Int, k::Int)
    duplicated = vcat([M[i,:] for _ ∈ 1:k]...)
    return [M[1:i-1, :]; duplicated; M[i+1:end, :]]
end

function duplicate_column(M::MatElem, i::Int, k::Int)
    duplicated = hcat(M[:,i] for _ ∈ 1:k)
    return [M[:, 1:i-1] duplicated M[:,i+1:end]]
end

