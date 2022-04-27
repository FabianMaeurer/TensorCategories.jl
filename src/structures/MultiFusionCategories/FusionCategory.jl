
mutable struct RingCategory <: Category
    base_ring::Field
    simples::Int64
    simples_names::Vector{String}
    ass::Array{<:MatElem,4}
    braiding::Function
    tensor_product::Array{Int,3}
    spherical::Vector
    twist::Vector

    function RingCategory(F::Field, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1])])
        C = New(F, length(mult[1]), names)
        C.tensor_product = mult
        #C.ass = [id(⊗(X,Y,Z)) for X ∈ simples(C), Y ∈ simples(C), Z ∈ simples(C)]
        #C.dims = [1 for i ∈ 1:length(names)]
        return C
    end

    function RingCategory(F::Field, names::Vector{String})
        C = new(F,length(names), names)
        #C.dims = [1 for i ∈ 1:length(names)]
        return C
    end

end


struct RingObject <: Object
    parent::RingCategory
    components::Vector{Int}
end

struct RingMorphism <: Morphism
    domain::RingObject
    codomain::RingObject
    m::Vector{<:MatElem}
end


#-------------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------------

RingCategory(x...) = RingCategory(x...)

Morphism(X::RingObject, Y::RingObject, m::Vector) = RingMorphism(X,Y,m)

#-------------------------------------------------------------------------------
#   Setters/Getters
#-------------------------------------------------------------------------------

function set_tensor_product!(F::RingCategory, tensor::Array{Int,3})
    F.tensor_product = tensor
    n = size(tensor,1)
    F.ass = Array{MatElem,4}(undef,n,n,n,n)
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        F.ass[i,j,k,:] = matrices(id(F[i]⊗F[j]⊗F[k]))
    end
end

function set_braiding!(F::RingCategory, braiding::Function)
    F.braiding = braiding
end

function set_associator!(F::RingCategory, i::Int, j::Int, k::Int, ass::Vector{<:MatElem})
    F.ass[i,j,k,:] = ass
end

function set_ev!(F::RingCategory, ev::Vector)
    F.evals = ev
end

function set_coev!(F::RingCategory, coev::Vector)
    F.coevals = coev
end

function set_spherical!(F::RingCategory, sp::Vector)
    F.spherical = sp
end

function set_duals!(F::RingCategory, d::Vector)
    F.duals = d
end

function set_ribbon!(F::RingCategory, r::Vector)
    F.ribbon = r
end

function set_dims!(F::RingCategory, d::Vector)
    F.dims = d
end

function set_twist!(F::RingCategory, t::Vector)
    F.twist = t
end

# function set_ev!(F::RingCategory, ev::Vector)
#     F.ev = ev
# end
#
# function set_coev!(F::RingCategory, coev::Vector)
#     F.coev = coev
# end

dim(X::RingObject) = base_ring(X)(tr(id(X)))

(::Type{Int})(x::fmpq) = Int(numerator(x))


braiding(X::RingObject, Y::RingObject) = parent(X).braiding(X,Y)

function associator(X::RingObject, Y::RingObject, Z::RingObject)
    @assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"
    C = parent(X)
    F = base_ring(C)
    n = C.simples
    dom = X⊗Y⊗Z
    m = zero_morphism(zero(C),zero(C))

    table = C.tensor_product
    associator = C.ass

    # Order of summands in domain
    dom_order_temp = [(k, m1*m2*table[i,j,k],[i,j]) for k ∈ 1:n, (i,m1) ∈ zip(1:n,X.components), (j,m2) ∈ zip(1:n,Y.components)][:]
    filter!(e -> e[2] != 0, dom_order_temp)
    sort!(dom_order_temp, by = e -> e[1])
    dom_order = [(k,m1*m2*table[i,j,k], [id; j]) for k ∈ 1:n, (i,m1,id) ∈ dom_order_temp, (j,m2) ∈ zip(1:n,Z.components)][:]
    filter!(e -> e[2] != 0, dom_order)
    sort!(dom_order, by = e -> e[1])

    # Order of summands in codomain
    cod_order_temp = [(k, m1*m2*table[i,j,k], [i,j]) for k ∈ 1:n, (i,m1) ∈ zip(1:n,Y.components), (j,m2) ∈ zip(1:n,Z.components)][:]
    filter!(e -> e[2] != 0, cod_order_temp)
    sort!(cod_order_temp, by = e -> e[1])
    cod_order = [(k,m1*m2*table[i,j,k], [i; id]) for k ∈ 1:n, (i,m2) ∈ zip(1:n,X.components), (j,m1, id) ∈ cod_order_temp][:]
    filter!(e -> e[2] != 0, cod_order)
    sort!(cod_order)

    # Associator
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        for i2 ∈ 1:X[i], j2 ∈ 1:Y[j], k2 ∈ 1:Z[k]
            T = C[i]⊗C[j]⊗C[k]
            m = m ⊕ Morphism(T,T,associator[i,j,k,:])
        end
    end

    # Order of summands in associator
    ass_order_temp = [(k, m1*m2*table[i,j,k],[i,j]) for k ∈ 1:n, (i,m1) ∈ zip(1:n,X.components), (j,m2) ∈ zip(1:n,Y.components)][:]
    filter!(e -> e[2] != 0, ass_order_temp)
    ass_order = [(k,m1*m2*table[i,j,k], [id; j]) for k ∈ 1:n, (i,m1,id) ∈ ass_order_temp, (j,m2) ∈ zip(1:n,Z.components)][:]
    filter!(e -> e[2] != 0, ass_order)

    comp_maps = matrices(m)

    # Permutation matrices
    for i ∈ 1:n
        dom_i = filter(e -> e[1] == i, dom_order)
        cod_i = filter(e -> e[1] == i, cod_order)
        ass_i = filter(e -> e[1] == i, ass_order)

        c_ass = vector_permutation(dom_i,ass_i)

        dom_dims = [k for (_,k,_) ∈ dom_i]
        ass_dims = [k for (_,k,_) ∈ ass_i]
        cod_dims = [k for (_,k,_) ∈ cod_i]

        # Permutation dom -> associator
        ass_perm = zero(MatrixSpace(F,sum(dom_dims),sum(dom_dims)))
        j = 0

        for (k,d) ∈ zip(c_ass,dom_dims)

            nk = sum(ass_dims[1:k-1])

            for i ∈ 1:d
                ass_perm[j+i,nk+i] = F(1)
            end
            j = j+d
        end

        # Permutation associator -> cod
        cod_perm = zero(MatrixSpace(F,sum(dom_dims),sum(dom_dims)))

        c_cod = vector_permutation(ass_i,cod_i)
        j = 0
        for (k,d) ∈ zip(c_cod,ass_dims)
            nk = sum(cod_dims[1:k-1])

            for i ∈ 1:d
                cod_perm[j+i,nk+i] = F(1)
            end
            j = j+d
        end
        comp_maps[i] = ass_perm*comp_maps[i]*cod_perm

    end

    return Morphism(dom,dom, comp_maps)
end

function vector_permutation(A::Vector,B::Vector)
    temp = deepcopy(B)
    perm = Int[]
    for a ∈ A
        i = findall(e -> e == a, temp)
        j = filter(e -> !(e ∈ perm), i)[1]
        perm = [perm; j]
    end
    return perm
end



#-------------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------------
issemisimple(::RingCategory) = true

issimple(X::RingObject) = sum(X.components) == 1

==(X::RingObject, Y::RingObject) = parent(X) == parent(Y) && X.components == Y.components
==(f::RingMorphism, g::RingMorphism) = domain(f) == domain(g) && codomain(f) == codomain(g) && f.m == g.m

decompose(X::RingObject) = [(x,k) for (x,k) ∈ zip(simples(parent(X)), X.components) if k != 0]

inv(f::RingMorphism) = RingMorphism(codomain(f),domain(f), inv.(f.m))

id(X::RingObject) = RingMorphism(X,X, [one(MatrixSpace(base_ring(X),d,d)) for d ∈ X.components])

function compose(f::RingMorphism, g::RingMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"
    return RingMorphism(domain(f), codomain(g), [m*n for (m,n) ∈ zip(f.m,g.m)])
end

function +(f::RingMorphism, g::RingMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    RingMorphism(domain(f), codomain(f), [m + n for (m,n) ∈ zip(f.m,g.m)])
end

"""
    dual(X::RingObject)

Return the dual object of ``X``. An error is thrown if ``X`` is not rigid.
"""
function dual(X::RingObject)
    C = parent(X)

    # Dual of simple Object
    if issimple(X)
        # Check for rigidity
        i = findfirst(e -> e == 1, X.components)
        j = findall(e -> C.tensor_product[i,e,1] >= 1, 1:C.simples)
        if length(j) != 1
            throw(ErrorException("Object not rigid."))
        end
        return RingObject(C,[i == j[1] ? 1 : 0 for i ∈ 1:C.simples])
    end

    # Build dual from simple objects
    return dsum([dual(Y)^(X.components[i]) for (Y,i) ∈ zip(simples(C), 1:C.simples)])
end

function coev(X::RingObject) where T
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    if sum(X.components) == 0 return zero_morphism(one(C), X) end

    m = []

    for (x,k) ∈ zip(simples(C),X.components), y ∈ simples(C)

        if x == dual(y)
            c = [F(a==b) for a ∈ 1:k, b ∈ 1:k][:]
            m = [m; c]
        else
            c = [0 for _ ∈ 1:(x⊗y).components[1]]
            m = [m; c]
        end
    end

    mats = matrices(zero_morphism(one(C), X⊗DX))
    M = parent(mats[1])
    mats[1] = M(F.(m))
    return Morphism(one(C), X⊗DX, mats)
end

function ev(X::RingObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    # Simple Objects
    if issimple(X)
        # If X is simple
        e = basis(Hom(DX⊗X, one(C)))[1]
        # Scale ev
        f = (id(X)⊗e)∘associator(X,DX,X)∘(coev(X)⊗id(X))
        return inv(F(f))*e
    end

    m = elem_type(F)[]
    #Arbitrary Objects
    for (x,k) ∈ zip(simples(C),DX.components), y ∈ simples(C)
        if x == dual(y)
            c = F(ev(y)[1]).*([F(a==b) for a ∈ 1:k, b ∈ 1:k][:])
            m = [m; c]
        else
            c = [0 for _ ∈ 1:(x⊗y).components[1]]
            m = [m; c]
        end
    end

    mats = matrices(zero_morphism(X⊗DX, one(C)))
    M = parent(mats[1])
    mats[1] = M(F.(m))
    return Morphism(X⊗DX,one(C),mats)
end

function spherical(X::RingObject)
    C = parent(X)
    sp = C.spherical
    return dsum([x^k for (x,k) ∈ zip(sp, X.components)])
end


*(λ,f::RingMorphism) = RingMorphism(domain(f), codomain(f), λ .*f.m)

# function tr(f::RingMorphism)
#     sum(tr.(f.m))
# end

# function smatrix(C::RingCategory)
#     θ = C.twist
#     #[inv(θ(i))*inv(θ(j))*sum() i ∈ simples(C), j ∈ simples(C)]
# end

function getindex(f::RingMorphism, i)
    m = zero_morphism(domain(f),codomain(f)).m
    m[i] = f.m[i]
    simple = simples(parent(domain(f)))
    dom = simple[i]^domain(f).components[i]
    cod = simple[i]^codomain(f).components[i]
    return RingMorphism(dom,cod,m)
end

getindex(X::RingObject, i) = X.components[i]

function matrices(f::RingMorphism)
    f.m
end


function (F::Field)(f::RingMorphism)
    if !(domain(f) == codomain(f) && issimple(domain(f)))
        throw(ErrorException("Cannot convert Morphism to $F"))
    end
    i = findfirst(e -> e == 1, domain(f).components)
    return F(f.m[i][1,1])
end
#-------------------------------------------------------------------------------
#   Tensor Product
#-------------------------------------------------------------------------------

function tensor_product(X::RingObject, Y::RingObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    C = parent(X)
    n = C.simples
    T = [0 for i ∈ 1:n]

    Xc = X.components
    Yc = Y.components

    for (i,j) ∈ Base.product(1:n, 1:n)
        if (c = Xc[i]) != 0 && (d = Yc[j]) != 0
            coeffs = C.tensor_product[i,j,:]
            T = T .+ ((c*d) .* coeffs)
        end
    end

    return RingObject(C,T)
end

function tensor_product(f::RingMorphism, g::RingMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)
    C = parent(dom)

    h = zero_morphism(zero(C), zero(C))

    table = C.tensor_product
    simpl = simples(C)

    for i ∈ 1:C.simples, j ∈ 1:C.simples
        A = kronecker_product(f.m[i],g.m[j])
        d1,d2 = size(A)
        #if d1*d2 == 0 continue end
        for k ∈ 1:C.simples
            if table[i,j,k] > 0
                m = zero_morphism(simpl[k]^d1,simpl[k]^d2).m
                m[k] = A

                for _ ∈ 1:table[i,j,k]
                    h = h ⊕ RingMorphism(simpl[k]^d1,simpl[k]^d2, m)
                end

            end
        end
    end
    #dom_left = dom.components - domain(h).components
    #cod_left = cod.components - codomain(h).components
    return h #⊕ zero_morphism(RingObject(C,dom_left), RingObject(C,cod_left))
end


one(C::RingCategory) = simples(C)[1]

#-------------------------------------------------------------------------------
#   Direct sum
#-------------------------------------------------------------------------------

function dsum(X::RingObject, Y::RingObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    return RingObject(parent(X), X.components .+ Y.components)
end

function dsum(f::RingMorphism, g::RingMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    F = base_ring(dom)
    m = zero_morphism(dom,cod).m
    for i ∈ 1:parent(dom).simples
        mf,nf = size(f.m[i])
        mg,ng = size(g.m[i])
        z1 = zero(MatrixSpace(F,mf,ng))
        z2 = zero(MatrixSpace(F,mg,nf))
        m[i] = [f.m[i] z1; z2 g.m[i]]
    end
    return RingMorphism(dom,cod, m)
end


zero(C::RingCategory) = RingObject(C,[0 for i ∈ 1:C.simples])

function zero_morphism(X::RingObject, Y::RingObject)
    return RingMorphism(X,Y,[zero(MatrixSpace(base_ring(X), cX, cY)) for (cX,cY) ∈ zip(X.components, Y.components)])
end
#-------------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------------

function simples(C::RingCategory)
    n = C.simples
    [RingObject(C, [i == j ? 1 : 0 for j ∈ 1:n]) for i ∈ 1:n]
end

function getindex(C::RingCategory, i)
    RingObject(C,[i == j ? 1 : 0 for j ∈ 1:C.simples])
end

#-------------------------------------------------------------------------------
#   Examples
#-------------------------------------------------------------------------------

function Ising()
    Qx,x = QQ["x"]
    F,a = NumberField(x^2-2, "√2")
    C = RingCategory(F,["1", "χ", "X"])
    M = zeros(Int,3,3,3)

    M[1,1,:] = [1,0,0]
    M[1,2,:] = [0,1,0]
    M[1,3,:] = [0,0,1]
    M[2,1,:] = [0,1,0]
    M[2,2,:] = [1,0,0]
    M[2,3,:] = [0,0,1]
    M[3,1,:] = [0,0,1]
    M[3,2,:] = [0,0,1]
    M[3,3,:] = [1,1,0]

    set_tensor_product!(C,M)

    set_associator!(C,2,3,2, matrices(-id(C[3])))
    set_associator!(C,3,1,3, matrices(id(C[1])⊕(-id(C[2]))))
    set_associator!(C,3,2,3, matrices((-id(C[1]))⊕id(C[2])))
    z = zero(MatrixSpace(F,0,0))
    set_associator!(C,3,3,3, [z, z, inv(a)*matrix(F,[1 1; 1 -1])])

    set_spherical!(C, [id(s) for s ∈ simples(C)])

    a,b,c = simples(C)

    return C
end

#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct RingCatHomSpace<: HomSpace
    X::RingObject
    Y::RingObject
    basis::Vector{RingMorphism}
    parent::VectorSpaces
end

function Hom(X::RingObject, Y::RingObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Xi, Yi = X.components, Y.components
    F = base_ring(X)

    d = sum([x*y for (x,y) ∈ zip(Xi,Yi)])

    if d == 0 return RingCatHomSpace(X,Y,RingMorphism[], VectorSpaces(F)) end

    basis = [zero_morphism(X,Y).m for i ∈ 1:d]
    next = 1
    for k ∈ 1:parent(X).simples

        for i ∈ 1:Xi[k], j ∈ 1:Yi[k]
            basis[next][k][i,j] = 1
            next = next + 1
        end
    end
    basis_mors = [RingMorphism(X,Y,m) for m ∈ basis]
    return RingCatHomSpace(X,Y,basis_mors, VectorSpaces(F))
end

function express_in_basis(f::RingMorphism, base::Vector)
    F = base_ring(domain(f))
    A = Array{elem_type(F),2}(undef,length(base),0)
    b = []
    for g ∈ base
        y = []
        for m ∈ g.m
            y = [y; [x for x ∈ m][:]]
        end
        A = [A y]
    end
    for m ∈ f.m
        b = [b; [x for x ∈ m][:]]
    end

    return [i for  i ∈ solve_left(transpose(matrix(F,A)), MatrixSpace(F,1,length(b))(F.(b)))][:]
end


#-------------------------------------------------------------------------------
#   Pretty Printing
#-------------------------------------------------------------------------------

function show(io::IO, C::RingCategory)
    print(io, "Fusion Category with $(C.simples) simple objects")
end

function show(io::IO, X::RingObject)
    coeffs = X.components

    if sum(coeffs) == 0
        print(io,"0")
        return
    end

    strings = parent(X).simples_names
    non_zero_coeffs = coeffs[coeffs .> 0]
    non_zero_strings = strings[coeffs .> 0]

    disp = non_zero_coeffs[1] == 1 ? "$(non_zero_strings[1])" : "$(non_zero_coeffs[1])⋅$(non_zero_strings[1])"

    for (Y,d) ∈ zip(non_zero_strings[2:end], non_zero_coeffs[2:end])
        disp = d == 1 ? disp*" ⊕ $Y" : disp*" ⊕ $(d)⋅$Y"
    end
    print(io,disp)
end

function show(io::IO, f::RingMorphism)
    print(io, """Morphism with
Domain: $(domain(f))
Codomain: $(codomain(f))
Matrices: """)
print(io, join(["$(m)" for m ∈ f.m], ", "))
end

#-------------------------------------------------------------------------------
#   Utility
#-------------------------------------------------------------------------------
