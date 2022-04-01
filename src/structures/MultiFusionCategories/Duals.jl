function ev_coev(X::Object)
    Y = dual(X)
    ev_dom = X⊗Y
    coev_cod = Y⊗X
    C = parent(X)

    F = base_ring(X)

    H_ev = Hom(ev_dom, one(C))
    H_coev = Hom(one(C), coev_cod)

    r,s = dim(H_coev), dim(H_ev)

    R,x = PolynomialRing(F,r+s)

    f, g = basis(H_coev), basis(H_ev)
    EndX = End(X)
    base = basis(EndX)
    eqs = [zero(R) for b ∈ base]

    for i ∈ 1:r, j ∈ 1:s
        eqs = eqs .+ (x[i]*x[r+j]) .*express_in_basis((g[i]⊗id(X))∘inv(associator(X,dual(X),X))∘(id(X)⊗f[j]),base)
    end

    I = ideal(eqs .- express_in_basis(id(X),base))

    @show dim(I)
    coeffs = recover_solutions(msolve(I),F)[1]
    return sum(coeffs[r+1:end] .* g), sum(coeffs[1:r] .* f)
end

ev(X::Object) = ev_coev(X)[1]
coev(X::Object) = ev_coev(X)[2]
