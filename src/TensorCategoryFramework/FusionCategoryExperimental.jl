
mutable struct NonStrictSixJCategory <: Category
    base_ring::Field
    simples::Int64
    simples_names::Vector{String}
    ass::Array{<:MatElem,4}
    braiding::Array{MatElem,3}
    tensor_product::Array{Int,3}
    spherical::Vector
    twist::Vector

    function NonStrictSixJCategory(F::Field, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1])])
        C = New(F, length(mult[1]), names)
        C.tensor_product = mult
        #C.ass = [id(⊗(X,Y,Z)) for X ∈ simples(C), Y ∈ simples(C), Z ∈ simples(C)]
        #C.dims = [1 for i ∈ 1:length(names)]
        return C
    end

    function NonStrictSixJCategory(F::Field, names::Vector{String})
        C = new(F,length(names), names)
        #C.dims = [1 for i ∈ 1:length(names)]
        return C
    end

end


struct NonStrictSixJObject <: Object
    parent::NonStrictSixJCategory
    components::Vector{Int}
end

struct NonStrictSixJMorphism <: Morphism
    domain::NonStrictSixJObject
    codomain::NonStrictSixJObject
    m::MatElem
end


#-------------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------------

NonStrictSixJCategory(x...) = NonStrictSixJCategory(x...)

Morphism(X::NonStrictSixJObject, Y::NonStrictSixJObject, m::MatElem) = NonStrictSixJMorphism(X,Y,m)


#-------------------------------------------------------------------------------
#   Setters/Getters
#-------------------------------------------------------------------------------

function set_tensor_product!(F::NonStrictSixJCategory, tensor::Array{Int,3})
    F.tensor_product = tensor
    n = size(tensor,1)
    F.ass = Array{MatElem,4}(undef,n,n,n,n)
    for i ∈ 1:n, j ∈ 1:n, k ∈ 1:n
        F.ass[i,j,k,:] = matrices(id(F[i]⊗F[j]⊗F[k]))
    end
end

function set_braiding!(F::NonStrictSixJCategory, braiding::Function)
    F.braiding = braiding
end

function set_associator!(F::NonStrictSixJCategory, i::Int, j::Int, k::Int, ass::Vector{<:MatElem})
    F.ass[i,j,k,:] = ass
end


function set_spherical!(F::NonStrictSixJCategory, sp::Vector)
    F.spherical = sp
end


function set_ribbon!(F::NonStrictSixJCategory, r::Vector)
    F.ribbon = r
end

function set_twist!(F::NonStrictSixJCategory, t::Vector)
    F.twist = t
end

dim(X::NonStrictSixJObject) = base_ring(X)(tr(id(X)))

(::Type{Int})(x::fmpq) = Int(numerator(x))


braiding(X::NonStrictSixJObject, Y::NonStrictSixJObject) = parent(X).braiding(X,Y)



function associator(X::NonStrictSixJObject, Y::NonStrictSixJObject, Z::NonStrictSixJObject)
    C = parent(X)
    if is_simple(X) && is_simple(Y) && is_simple(Z)
        dom = (X⊗Y)⊗Z
        cod = X⊗(Y⊗Z)
        mat = matrix(zero_morphism(dom,cod))
        for i ∈ eachindex(simples(C))
            dom_index = findall(e -> e == i, dom.components)
            cod_index = findall(e -> e == i, cod.components)
            mat[cod_index, dom_index] = C.ass[X.components[1], Y.components[1], Z.components[1], i]
        end
        ass = Morphism(NonStrictSixJObject(C, dom.components), NonStrictSixJObject(C, cod.components ), mat)
        return ass
    end

    return  Morphism((X ⊗ Y) ⊗ Z, X ⊗ (Y ⊗ Z), diagonal_matrix([matrix(associator(C[i],C[j],C[k])) for k ∈ Z.components, j ∈ Y.components, i ∈ X.components][:]))
end




function vector_permutation(A::Vector,B::Vector)
    perm = Int[]
    for a ∈ A
        i = findall(e -> e == a, B)
        j = filter(e -> !(e ∈ perm), i)[1]
        perm = [perm; j]
    end
    return perm
end



#-------------------------------------------------------------------------------
#   Functionality
#-------------------------------------------------------------------------------
is_semisimple(::NonStrictSixJCategory) = true

is_simple(X::NonStrictSixJObject) = length(X.components) == 1

==(X::NonStrictSixJObject, Y::NonStrictSixJObject) = parent(X) == parent(Y) && X.components == Y.components
==(f::NonStrictSixJMorphism, g::NonStrictSixJMorphism) = domain(f) == domain(g) && codomain(f) == codomain(g) && f.m == g.m

#decompose(X::NonStrictSixJObject) = [(x,k) for (x,k) ∈ zip(simples(parent(X)), X.components) if k != 0]

inv(f::NonStrictSixJMorphism) = NonStrictSixJMorphism(codomain(f),domain(f), inv(f.m))

id(X::NonStrictSixJObject) = NonStrictSixJMorphism(X,X, one(MatrixSpace(base_ring(X),length(X.components),length(X.components))))

function compose(f::NonStrictSixJMorphism, g::NonStrictSixJMorphism)
    @assert codomain(f) == domain(g) 
    
    NonStrictSixJMorphism(domain(f), codomain(g), matrix(f)*matrix(g))
end

function +(f::NonStrictSixJMorphism, g::NonStrictSixJMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    NonStrictSixJMorphism(domain(f), codomain(f), matrix(f) + matrix(g))
end

"""
    dual(X::NonStrictSixJObject)

Return the dual object of ``X``. An error is thrown if ``X`` is not rigid.
"""
function dual(X::NonStrictSixJObject)
    C = parent(X)

    # Dual of simple Object
    if issimple(X)
        # Check for rigidity
        i = findfirst(e -> e == 1, X.components)
        j = findall(e -> C.tensor_product[i,e,1] >= 1, 1:C.simples)
        if length(j) != 1
            throw(ErrorException("Object not rigid."))
        end
        return NonStrictSixJObject(C,[i == j[1] ? 1 : 0 for i ∈ 1:C.simples])
    end

    # Build dual from simple objects
    return direct_sum([dual(Y)^(X.components[i]) for (Y,i) ∈ zip(simples(C), 1:C.simples)])
end

function coev(X::NonStrictSixJObject) where T
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

function ev(X::NonStrictSixJObject)
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




*(λ,f::NonStrictSixJMorphism) = NonStrictSixJMorphism(domain(f), codomain(f), λ *f.m)


# function getindex(f::NonStrictSixJMorphism, i)
#     m = zero_morphism(domain(f),codomain(f)).m
#     m[i] = f.m[i]
#     simple = simples(parent(domain(f)))
#     dom = simple[i]^domain(f).components[i]
#     cod = simple[i]^codomain(f).components[i]
#     return NonStrictSixJMorphism(dom,cod,m)
# end

getindex(X::NonStrictSixJObject, i) = X.components[i]

function matrix(f::NonStrictSixJMorphism)
    return f.m
end

function matrices(f::NonStrictSixJMorphism)
    n = parent(domain(f)).simples
    F = base_ring(f)

    mats = MatElem[]
    for i ∈ 1:n
        dom_index = findall(e -> e == i, domain(f).components)
        cod_index = findall(e -> e == i, codomain(f).components)
        mats = [mats; matrix(f)[cod_index, dom_index]]
    end

    return mats
end


function (F::Field)(f::NonStrictSixJMorphism)
    if !(domain(f) == codomain(f) && issimple(domain(f)))
        throw(ErrorException("Cannot convert Morphism to $F"))
    end
    return F(f.m[1,1])
end

spherical(X::NonStrictSixJObject) = id(X)
#-------------------------------------------------------------------------------
#   Tensor Product
#-------------------------------------------------------------------------------

function tensor_product(X::NonStrictSixJObject, Y::NonStrictSixJObject)
    #@assert parent(X) == parent(Y) "Mismatching parents"
    C = parent(X)
    n = C.simples
    T = Int[]

    table = C.tensor_product
    for i ∈ X.components, j ∈ Y.components
        for k ∈ 1:n
            T = [T; [k for _ ∈ 1:table[i,j,k]]]
        end
    end

    return NonStrictSixJObject(C,T)
end
function tensor_product(f::NonStrictSixJMorphism, g::NonStrictSixJMorphism)
    dom = domain(f) ⊗ domain(g)
    cod = codomain(f) ⊗ codomain(g)
    C = parent(dom)

    table = C.tensor_product

    mat = zero(MatrixSpace(base_ring(C), 0, length(cod.components)))

    mat_f, mat_g = matrix(f), matrix(g)
    nf,mf = size(mat_f)
    ng,mg = size(mat_g)

    for row_f ∈ 1:nf, row_g ∈ 1:ng
        dom_k = sum(table[domain(f).components[row_f], domain(g).components[row_g], :])
        temp_mat = zero(MatrixSpace(base_ring(C), dom_k, 0))

        for col_f ∈ 1:mf, col_g ∈ 1:mg
            cod_k = sum(table[codomain(f).components[col_f], codomain(g).components[col_g], :])
            temp_mat = [temp_mat diagonal_matrix(mat_f[row_f,col_f]*mat_g[row_g,col_g], dom_k, cod_k)]
        end
        size(mat)
        size(temp_mat)
        mat = [mat; temp_mat]
    end

    return Morphism(dom, cod, mat)
 
end


one(C::NonStrictSixJCategory) = simples(C)[1]

#-------------------------------------------------------------------------------
#   Direct sum
#-------------------------------------------------------------------------------

function direct_sum(X::NonStrictSixJObject, Y::NonStrictSixJObject)
    S = NonStrictSixJObject(parent(X), [X.components; Y.components])
    ix_mat = matrix(zero_morphism(X,S))
    iy_mat = matrix(zero_morphism(Y,S))
    px_mat = matrix(zero_morphism(S,X))
    py_mat = matrix(zero_morphism(S,Y))

    for i ∈ 1:length(X.components)
        ix_mat[i,i] = 1
        px_mat[i,i] = 1
    end

    for i ∈ 1:length(Y.components)
        iy_mat[i, length(X.components) + i] = 1
        py_mat[length(X.components) + i,i] = 1
    end

    ix = Morphism(X,S, ix_mat)
    px = Morphism(S,X, px_mat)
    iy = Morphism(Y,S, iy_mat)
    py = Morphism(S,Y, py_mat)

    return S,[ix,iy],[px,py]
end

function direct_sum(f::NonStrictSixJMorphism, g::NonStrictSixJMorphism)
    dom = domain(f) ⊕ domain(g)
    cod = codomain(f) ⊕ codomain(g)
    m = diagonal_matrix(f.m, g.m)
    return Morphism(dom,cod,m)
end


zero(C::NonStrictSixJCategory) = NonStrictSixJObject(C,[])

function zero_morphism(X::NonStrictSixJObject, Y::NonStrictSixJObject)
    return NonStrictSixJMorphism(X,Y,zero(MatrixSpace(base_ring(X), length(X.components), length(Y.components))))
end

function is_isomorphic(X::NonStrictSixJObject, Y::NonStrictSixJObject)
    if sort(X.components) != sort(Y.components)
        return false, nothing
    else
    F = base_ring(X)
    σ₁ = sortperm(X.components)
    σ₂ = sortperm(Y.components)
    permutation = permutation_matrix(F,σ₁)*inv(permutation_matrix(F,σ₂))
        return true, Morphism(X,Y,permutation)
    end
end
#-------------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------------

function simples(C::NonStrictSixJCategory)
    n = C.simples
    [NonStrictSixJObject(C, [i]) for i ∈ 1:n]
end

function getindex(C::NonStrictSixJCategory, i::Int)
    NonStrictSixJObject(C,[i])
end

#-------------------------------------------------------------------------------
#   Kernel and Cokernel
#-------------------------------------------------------------------------------

function kernel(f::NonStrictSixJMorphism)
    C = parent(domain(f))
    kernels = [kernel(Morphism(m)) for m ∈ f.m]
    mats = [matrix(m) for (k,m) ∈ kernels]
    ker = NonStrictSixJObject(C,[dim(k) for (k,m) ∈ kernels])

    return ker, Morphism(ker, domain(f), mats)
end


function left_inverse(f::NonStrictSixJMorphism)
    inverses = [left_inverse(Morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return Morphism(codomain(f), domain(f), mats)
end

#-------------------------------------------------------------------------------
#   Examples
#-------------------------------------------------------------------------------

function NonStrictIsing()
    Qx,x = QQ["x"]
    F,a = NumberField(x^2-2, "√2")
    C = NonStrictSixJCategory(F,["𝟙", "χ", "X"])
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
    set_associator!(C,3,1,3, matrices(id(C[1])⊕(id(C[2]))))
    set_associator!(C,3,2,3, matrices((id(C[1]))⊕(-id(C[2]))))
    z = zero(MatrixSpace(F,0,0))
    set_associator!(C,3,3,3, [z, z, inv(a)*matrix(F,[1 1; 1 -1])])

    set_spherical!(C, [F(1) for s ∈ simples(C)])


    return C
end

#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct NonStrictSixJHomSpace<: AbstractHomSpace
    X::NonStrictSixJObject
    Y::NonStrictSixJObject
    basis::Vector{NonStrictSixJMorphism}
    parent::VectorSpaces
end

function Hom(X::NonStrictSixJObject, Y::NonStrictSixJObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Xi, Yi = X.components, Y.components
    F = base_ring(X)

    d = sum([x*y for (x,y) ∈ zip(Xi,Yi)])

    if d == 0 return SixJHomSpace(X,Y,NonStrictSixJMorphism[], VectorSpaces(F)) end

    basis = [zero_morphism(X,Y).m for i ∈ 1:d]
    next = 1
    for k ∈ 1:parent(X).simples

        for i ∈ 1:Xi[k], j ∈ 1:Yi[k]
            basis[next][k][i,j] = 1
            next = next + 1
        end
    end
    basis_mors = [NonStrictSixJMorphism(X,Y,m) for m ∈ basis]
    return SixJHomSpace(X,Y,basis_mors, VectorSpaces(F))
end

function express_in_basis(f::NonStrictSixJMorphism, base::Vector)
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

function show(io::IO, C::NonStrictSixJCategory)
    print(io, "Fusion Category with $(C.simples) simple objects")
end

function show(io::IO, X::NonStrictSixJObject)
    x_comps = X.components
    coeffs = [length(x_comps[x_comps .== k]) for k ∈ 1:parent(X).simples]

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

function show(io::IO, f::NonStrictSixJMorphism)
    print(io, """Morphism with
Domain: $(domain(f))
Codomain: $(codomain(f))
Matrices: """)
print(io, join(["$(m)" for m ∈ matrices(f)], ", "))
end

#-------------------------------------------------------------------------------
#   Utility
#-------------------------------------------------------------------------------
