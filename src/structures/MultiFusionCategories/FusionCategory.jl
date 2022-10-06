
mutable struct RingCategory <: Category
    base_ring::Field
    simples::Int64
    simples_names::Vector{String}
    ass::Array{<:MatElem,4}
    braiding::Function
    tensor_product::Array{Int,3}
    spherical::Vector
    twist::Vector
    one::Vector{Int}

    function RingCategory(F::Field, mult::Array{Int,3}, names::Vector{String} = ["X$i" for i ∈ 1:length(mult[1,1,:])])
        C = new(F, length(mult[1,1,:]), names)
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


struct RingCatObject <: Object
    parent::RingCategory
    components::Vector{Int}
end

struct RingCatMorphism <: Morphism
    domain::RingCatObject
    codomain::RingCatObject
    m::Vector{<:MatElem}
end

#-------------------------------------------------------------------------------
#   Constructors
#-------------------------------------------------------------------------------

#RingCategory(x...) = RingCategory(x...)

Morphism(X::RingCatObject, Y::RingCatObject, m::Vector) = RingCatMorphism(X,Y,m)

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

set_associator!(F::RingCategory, ass::Array{M,4}) where M <: MatElem = F.ass = ass
function set_associator!(F::RingCategory, i::Int, j::Int, k::Int, ass::Vector{<:MatElem})
    F.ass[i,j,k,:] = ass
end

function set_spherical!(F::RingCategory, sp::Vector)
    F.spherical = sp
end

function set_one!(F::RingCategory, v::Vector{Int}) 
    F.one = v
end 

function set_ribbon!(F::RingCategory, r::Vector)
    F.ribbon = r
end

function set_twist!(F::RingCategory, t::Vector)
    F.twist = t
end


dim(X::RingCatObject) = base_ring(X)(tr(id(X)))

(::Type{Int})(x::fmpq) = Int(numerator(x))


braiding(X::RingCatObject, Y::RingCatObject) = parent(X).braiding(X,Y)


# function associator(X::RingCatObject, Y::RingCatObject, Z::RingCatObject)
#     @assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"
#
#     C = parent(X)
#     F = base_ring(C)
#     n = C.simples
#     dom = X⊗Y⊗Z
#
#     table = C.tensor_product
#     C_associator = C.ass
#
#     #---------------------------------
#     # associators on simple objects
#     #---------------------------------
#     if issimple(X) && issimple(Y) && issimple(Z)
#         i = findfirst(e -> e ≠ 0, X.components)
#         j = findfirst(e -> e ≠ 0, Y.components)
#         k = findfirst(e -> e ≠ 0, Z.components)
#         return Morphism(X⊗Y⊗Z, X⊗Y⊗Z, C_associator[i,j,k,:])
#     end
#
#     #---------------------------------
#     # associators for arbitrary objects
#     #---------------------------------
#     simple_objects = simples(parent(X))
#
#     #-------------------------------------
#     # Order of summands in domain
#     #-------------------------------------
#     ids = ones(Int,n)
#     domain_order = [[] for _ ∈ 1:n]
#
#     for i ∈ (X⊗Y).components, j ∈ Z.components
#
#
# end

# function associator(X::RingCatObject, Y::RingCatObject, Z::RingCatObject)
#     C = parent(X)
#     n = C.simples
#     simple_objects = simples(C)
#     #---------------------------------
#     # associators on simple objects
#     #---------------------------------
#     if issimple(X) && issimple(Y) && issimple(Z)
#         i = findfirst(e -> e ≠ 0, X.components)
#         j = findfirst(e -> e ≠ 0, Y.components)
#         k = findfirst(e -> e ≠ 0, Z.components)
#         return Morphism(X⊗Y⊗Z, X⊗Y⊗Z, C.ass[i,j,k,:])
#     end

#     X_summands = vcat([[(s,[k,l]) for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
#     Y_summands = vcat([[(s,[k,l]) for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
#     Z_summands = vcat([[(s,[k,l]) for l ∈ 1:Z.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)

#     m = zero_morphism(zero(C),zero(C))
#     for (x,_) ∈ X_summands, (y,_) ∈ Y_summands, (z,_) ∈ Z_summands
#         m = m ⊕ associator(x,y,z)
#     end

#     #=------------------------------------------------
#         Correct permutations
#     ------------------------------------------------=#
#     mats = matrices(id(X⊗Y⊗Z))

#     left_mors = [(f⊗g)⊗h  for f ∈ basis(End(X)), g ∈ basis(End(Y)), h ∈ basis(End(Z))][:]
#     right_mors = [f⊗(g⊗h)  for f ∈ basis(End(X)), g ∈ basis(End(Y)), h ∈ basis(End(Z))][:]
#     k = 1
#     for q in left_mors
#         right_mats = [matrices(w)[k] for w in right_mors]
#         j = findall(e -> e == matrices(q)[k], right_mats)
#         @show j
#         @show matrices(q)[k]
#     end

#     for f ∈ basis(End(X)), g ∈ basis(End(Y)), h ∈ basis(End(Z))
#         left_morphism  = m ∘ ((f ⊗ g) ⊗ h) ∘ inv(m)
#         right_morphism = (f ⊗ (g ⊗ h))
#         left_mats  = matrices(left_morphism)
#         right_mats = matrices(right_morphism)
#         for i ∈ 1:n
#             left_m  = left_mats[i]
#             right_m = right_mats[i]
#             if prod(size(left_m)) == 0 
#                 continue 
#             end

#             @show corr_mat = similarity_matrix(right_m, left_m)
            
#             mats[i] =  mats[i] * corr_mat
#             left_mats[i] = left_m * corr_mat
#             right_mats[i] = corr_mat * right_m
#         end
#     end
#     correction = Morphism(domain(m), codomain(m), mats)
#     for f ∈ basis(End(X)), g ∈ basis(End(Y)), h ∈ basis(End(Z))
#         @show correction ∘ ((f⊗g)⊗h) == (f⊗(g⊗h)) ∘ correction
#     end

#     return  correction ∘ m
# end


"""
    similarity_matrix(m::MatElem, n::MatElem)

Return matrix ```P``` such that ```PmP^{-1} = n```.
"""
function similarity_matrix(m::MatElem, n::MatElem)
    J₁, S = jordan_normal_form(m)
    J₂, T = jordan_normal_form(n)
    @assert J₁ == J₂ "Not same Jordan form"
    return inv(T)*S
end


#=-------------------------------------------------
    best associator so far 
-------------------------------------------------=#

function associator(X::RingCatObject, Y::RingCatObject, Z::RingCatObject)
    @assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"

    C = parent(X)
    F = base_ring(C)
    n = C.simples
    dom = X⊗Y⊗Z

    C_associator = C.ass

    #---------------------------------
    # associators on simple objects
    #---------------------------------
    if issimple(X) && issimple(Y) && issimple(Z)
        i = findfirst(e -> e ≠ 0, X.components)
        j = findfirst(e -> e ≠ 0, Y.components)
        k = findfirst(e -> e ≠ 0, Z.components)
        return Morphism(X⊗Y⊗Z, X⊗Y⊗Z, C_associator[i,j,k,:])
    end

    #---------------------------------
    # associators for arbitrary objects
    #---------------------------------
    simple_objects = simples(parent(X))

    X_summands = vcat([[s for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Y_summands = vcat([[s for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    Z_summands = vcat([[s for l ∈ 1:Z.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)

    #=-------------------------------------------------
        Distribution 
    -------------------------------------------------=#

    # Before
    distr_before = distribute_left(X_summands, Y) ⊗ id(Z)
    distr_before = (dsum([distribute_right(Xᵢ,Y_summands) for Xᵢ ∈ X_summands]...)⊗id(Z)) ∘ distr_before
    distr_before = distribute_left([Xᵢ⊗Yⱼ for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:], Z) ∘ distr_before
    distr_before = dsum([distribute_right(Xᵢ⊗Yⱼ,Z_summands) for Yⱼ ∈ Y_summands, Xᵢ ∈ X_summands][:]...) ∘ distr_before

    # After
    distr_after = id(X)⊗distribute_right(Y, Z_summands)
    distr_after = (id(X)⊗dsum([distribute_left(Y_summands,Zₖ) for Zₖ ∈ Z_summands]...)) ∘ distr_after
    distr_after = distribute_right(X, [Yⱼ⊗Zₖ for  Zₖ ∈ Z_summands, Yⱼ ∈ Y_summands][:]) ∘ distr_after
    distr_after = dsum([distribute_left(X_summands, Yⱼ⊗Zₖ) for  Zₖ ∈ Z_summands, Yⱼ ∈ Y_summands][:]...) ∘ distr_after


    #-----------------------------------
    # Associator morphism
    #-----------------------------------
    m = zero_morphism(zero(C),zero(C))
    for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands
        m = m ⊕ associator(x,y,z)
    end
    @show m
    return inv(distr_after) ∘ m ∘ distr_before


    # X_summands = vcat([[(s,[k,l]) for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    # Y_summands = vcat([[(s,[k,l]) for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    # Z_summands = vcat([[(s,[k,l]) for l ∈ 1:Z.components[k]] for (k,s) ∈ zip(1:n, simple_objects)]...)
    # #-------------------------------------
    # # Order of summands in domain
    # #-------------------------------------
    # domain_order_temp = []
    # for (x, x_id) ∈ X_summands, (y, y_id) ∈ Y_summands
    #      for (s,k) ∈ zip(simple_objects, (x⊗y).components)
    #          append!(domain_order_temp, [(s, [x_id; y_id]) for l ∈ 1:k])
    #      end
    # end
    # sort!(domain_order_temp, by = e -> findfirst(k -> k != 0, e[1].components))
    # domain_order = []
    # for (x, x_id) ∈ domain_order_temp, (z, z_id) ∈ Z_summands
    #     for (s,k) ∈ zip(simple_objects, (x⊗z).components)
    #         append!(domain_order, [(s, [x_id; z_id]) for l ∈ 1:k])
    #     end
    # end


    # #-----------------------------------
    # # Order of summands in codomain
    # #-----------------------------------
    # codomain_order_temp = []
    # for (y, y_id) ∈ Y_summands, (z, z_id) ∈ Z_summands
    #     for (s,k) ∈ zip(simple_objects, (y⊗z).components)
    #         append!(codomain_order_temp, [(s, [y_id; z_id]) for l ∈ 1:k])
    #     end
    # end
    # sort!(codomain_order_temp, by = e -> findfirst(k -> k != 0, e[1].components))
    # codomain_order = []
    # for (x, x_id) ∈ X_summands, (z, z_id) ∈ codomain_order_temp
    #     for (s,k) ∈ zip(simple_objects, (x⊗z).components)
    #         append!(codomain_order, [(s, [x_id; z_id]) for l ∈ 1:k])
    #     end
    # end

    # #-----------------------------------
    # # Order of summands in associator
    # #-----------------------------------
    # associator_order = []
    # for (x, x_id) ∈ X_summands, (y, y_id) ∈ Y_summands, (z, z_id) ∈ Z_summands
    #     for (s,k) ∈ zip(simple_objects, ((x⊗y)⊗z).components)
    #         append!(associator_order, [(s, [x_id; y_id; z_id]) for i ∈ 1:k])
    #     end
    # end

    # #-----------------------------------
    # # Associator morphism
    # #-----------------------------------
    # m = zero_morphism(zero(C),zero(C))
    # for (x,_) ∈ X_summands, (y,_) ∈ Y_summands, (z,_) ∈ Z_summands
    #     m = m ⊕ associator(x,y,z)
    # end


    # #-----------------------------------
    # # permutations
    # #-----------------------------------
    # comp_maps = matrices(m)

    # for i ∈ 1:n
    #     dom_i = filter(e -> e[1] == C[i], domain_order)
    #     cod_i = filter(e -> e[1] == C[i], codomain_order)
    #     ass_i = filter(e -> e[1] == C[i], associator_order)

    #     if length(dom_i) == 0 continue end
        
    #     c_ass = vector_permutation(dom_i, ass_i)

    #     # Permutation dom -> associator
    #     ass_perm = zero(MatrixSpace(F,length(dom_i),length(dom_i)))

    #     for (i,k) ∈ zip(1:length(c_ass), c_ass)
    #         ass_perm[i,k] = F(1)
    #     end
        
    #     # Permutation associator -> cod
    #     cod_perm = zero(MatrixSpace(F,length(cod_i),length(cod_i)))

    #     c_cod = vector_permutation(ass_i, cod_i)

    #     for (i,k) ∈ zip(1:length(c_cod), c_cod)
    #         cod_perm[i,k] = F(1)
    #     end
        
    #     comp_maps[i] = ass_perm*comp_maps[i]*cod_perm

    # end
    # return Morphism(dom,dom, comp_maps)

end

#=-------------------------------------------------
    Experimental 
-------------------------------------------------=#

# function associator(X::RingCatObject, Y::RingCatObject, Z::RingCatObject)
#     @assert parent(X) == parent(Y) == parent(Z) "Mismatching parents"

#     C = parent(X)
#     F = base_ring(C)
#     n = C.simples
#     dom = X⊗Y⊗Z

#     C_associator = C.ass

#     #---------------------------------
#     # associators on simple objects
#     #---------------------------------
#     if issimple(X) && issimple(Y) && issimple(Z)
#         i = findfirst(e -> e ≠ 0, X.components)
#         j = findfirst(e -> e ≠ 0, Y.components)
#         k = findfirst(e -> e ≠ 0, Z.components)
#         return Morphism(X⊗Y⊗Z, X⊗Y⊗Z, C_associator[i,j,k,:])
#     end

#     #=-------------------------------------------------
#         Associators on non-simple objects
#     -------------------------------------------------=#
#     R,x = PolynomialRing(F, dim(End(X⊗Y⊗Z)))
#     D = RingCategory(R, C.tensor_product)
#     set_associator!(D, [matrix(R,collect(m)) for m ∈ C_associator])

#     poly_simples = simples(D)
#     X_summands = vcat([[(s) for l ∈ 1:X.components[k]] for (k,s) ∈ zip(1:n, poly_simples)]...)
#     Y_summands = vcat([[(s) for l ∈ 1:Y.components[k]] for (k,s) ∈ zip(1:n, poly_simples)]...)
#     Z_summands = vcat([[(s) for l ∈ 1:Z.components[k]] for (k,s) ∈ zip(1:n, poly_simples)]...)

#     @show direct_associator = dsum([associator(x,y,z) for x ∈ X_summands, y ∈ Y_summands, z ∈ Z_summands][:])

#     PX = dsum(X_summands)
#     PY = dsum(Y_summands)
#     PZ = dsum(Z_summands)

#     P_XYZ = PX⊗PY⊗PZ
#     P_XYZ_dims = P_XYZ.components
#     Q_mats = MatElem[]
#     i = 1
#     for k ∈ P_XYZ_dims
#         m = matrix(R,k,k, reshape(x[i:i+k^2-1],k,k))
#         Q_mats = [Q_mats; m]
#         i = i + k^2
#     end
#     Q = Morphism(P_XYZ, P_XYZ, Q_mats)

#     # Set up basis for End(X), End(Y), End(Z) in d
#     PX_basis = [Morphism(PX,PX, [matrix(R, collect(m)) for m ∈ matrices(f)]) for f ∈ basis(End(X))]
#     PY_basis = [Morphism(PY,PY, [matrix(R, collect(m)) for m ∈ matrices(f)]) for f ∈ basis(End(Y))]
#     PZ_basis = [Morphism(PZ,PZ, [matrix(R, collect(m)) for m ∈ matrices(f)]) for f ∈ basis(End(Z))]

#     # Set up equations for naturality
#     nat_eqs = [(f⊗(g⊗h))∘Q∘direct_associator - Q∘direct_associator∘((f⊗g)⊗h) for f ∈ PX_basis, g ∈ PY_basis, h ∈ PZ_basis][:]
#     nat_eqs = vcat([matrices(eq) for eq ∈ nat_eqs]...)
#     nat_eqs = vcat([collect(eq)[:] for eq ∈ nat_eqs]...)
#     unique!(filter!(e -> e != 0, nat_eqs))

#     # Add permutation matrix conditions
#     perm_eqs_row = vcat([vcat([sum(m[i,:]) - 1 for i ∈ 1:k]...) for (m,k) ∈ zip(Q_mats, P_XYZ_dims)]...) 
#     perm_eqs_col = vcat([vcat([sum(m[:,i]) - 1 for i ∈ 1:k]...) for (m,k) ∈ zip(Q_mats, P_XYZ_dims)]...) 
#     ideal([nat_eqs; perm_eqs_col; perm_eqs_row])
#     ideal([nat_eqs; perm_eqs_col]) 
# end

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
issemisimple(::RingCategory) = true

issimple(X::RingCatObject) = sum(X.components) == 1

==(X::RingCatObject, Y::RingCatObject) = parent(X) == parent(Y) && X.components == Y.components
==(f::RingCatMorphism, g::RingCatMorphism) = domain(f) == domain(g) && codomain(f) == codomain(g) && f.m == g.m

decompose(X::RingCatObject) = [(x,k) for (x,k) ∈ zip(simples(parent(X)), X.components) if k != 0]

inv(f::RingCatMorphism) = RingCatMorphism(codomain(f),domain(f), inv.(f.m))

id(X::RingCatObject) = RingCatMorphism(X,X, [one(MatrixSpace(base_ring(X),d,d)) for d ∈ X.components])

function compose(f::RingCatMorphism, g::RingCatMorphism)
    @assert codomain(f) == domain(g) "Morphisms not compatible"
    return RingCatMorphism(domain(f), codomain(g), [m*n for (m,n) ∈ zip(f.m,g.m)])
end

function +(f::RingCatMorphism, g::RingCatMorphism)
    @assert domain(f) == domain(g) && codomain(f) == codomain(g) "Not compatible"
    RingCatMorphism(domain(f), codomain(f), [m + n for (m,n) ∈ zip(f.m,g.m)])
end

"""
    dual(X::RingCatObject)

Return the dual object of ``X``. An error is thrown if ``X`` is not rigid.
"""
function dual(X::RingCatObject)
    C = parent(X)

    # Dual of simple Object
    if issimple(X)
        # Check for rigidity
        i = findfirst(e -> e == 1, X.components)
        j = []
        for k ∈ 1:C.simples 
            if C.one[k] == 1
                j = [j; findall(e -> C.tensor_product[i,e,k] >= 1, 1:C.simples)]
            end
        end
        if length(j) != 1
            throw(ErrorException("Object not rigid."))
        end
        return RingCatObject(C,[i == j[1] ? 1 : 0 for i ∈ 1:C.simples])
    end

    # Build dual from simple objects
    return dsum([dual(Y)^(X.components[i]) for (Y,i) ∈ zip(simples(C), 1:C.simples)])
end

function coev(X::RingCatObject)
    if issimple(X)
        return simple_objects_coev(X)
    end
end

function ev(X::RingCatObject)
    if issimple(X)
        return simple_objects_ev(X)
    end
end

function simple_objects_coev(X::RingCatObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    cod = X ⊗ DX

    if sum(X.components) == 0 return zero_morphism(one(C), X) end

    mats = [diagonal_matrix(F(1),n,m) for (n,m) ∈ zip(C.one, cod.components)]

    return Morphism(one(C), cod, mats)

    # m = []

    # for (x,k) ∈ zip(simples(C),X.components), y ∈ simples(C)

    #     if x == dual(y)
    #         c = [F(a==b) for a ∈ 1:k, b ∈ 1:k][:]
    #         m = [m; c]
    #     else
    #         c = [0 for _ ∈ 1:(x⊗y).components[1]]
    #         m = [m; c]
    #     end
    # end

    # mats = matrices(zero_morphism(one(C), X⊗DX))
    # M = parent(mats[1])
    # mats[1] = M(F.(m))
    # return Morphism(one(C), X⊗DX, mats)
end

function simple_objects_ev(X::RingCatObject)
    DX = dual(X)
    C = parent(X)
    F = base_ring(C)

    dom = DX ⊗ X

    if sum(X.components) == 0 return zero_morphism(X,one(C)) end

    mats = [diagonal_matrix(F(1),n,m) for (n,m) ∈ zip(dom.components, C.one)]

    unscaled_ev = Morphism(dom, one(C), mats)

    factor = F((id(X)⊗unscaled_ev)∘associator(X,DX,X)∘(coev(X)⊗id(X)))


    return inv(factor) * unscaled_ev

    # # Simple Objects
    # if issimple(X)
    #     # If X is simple
    #     e = basis(Hom(DX⊗X, one(C)))[1]
    #     # Scale ev
    #     f = (id(X)⊗e)∘associator(X,DX,X)∘(coev(X)⊗id(X))
    #     return inv(F(f))*e
    # end

    # m = elem_type(F)[]
    # #Arbitrary Objects
    # for (x,k) ∈ zip(simples(C),DX.components), y ∈ simples(C)
    #     if x == dual(y)
    #         c = F(ev(y)[1]).*([F(a==b) for a ∈ 1:k, b ∈ 1:k][:])
    #         m = [m; c]
    #     else
    #         c = [0 for _ ∈ 1:(x⊗y).components[1]]
    #         m = [m; c]
    #     end
    # end

    # mats = matrices(zero_morphism(X⊗DX, one(C)))
    # M = parent(mats[1])
    # mats[1] = M(F.(m))
    # return Morphism(X⊗DX,one(C),mats)
end

function spherical(X::RingCatObject)
    C = parent(X)
    F = base_ring(C)
    sp = C.spherical
    mats = [diagonal_matrix(θ, k) for (θ,k) ∈ zip(sp, X.components)]
    return Morphism(X,X,mats)
end


*(λ,f::RingCatMorphism) = RingCatMorphism(domain(f), codomain(f), λ .*f.m)

# function tr(f::RingCatMorphism)
#     sum(tr.(f.m))
# end

# function smatrix(C::RingCategory)
#     θ = C.twist
#     #[inv(θ(i))*inv(θ(j))*sum() i ∈ simples(C), j ∈ simples(C)]
# end

function getindex(f::RingCatMorphism, i)
    simple = simples(parent(domain(f)))
    dom = simple[i]^domain(f).components[i]
    cod = simple[i]^codomain(f).components[i]
    m = zero_morphism(dom,cod).m
    m[i] = f.m[i]
    return RingCatMorphism(dom,cod,m)
end

getindex(X::RingCatObject, i::Int64) = X.components[i]

function matrices(f::RingCatMorphism)
    f.m
end

function matrix(f::RingCatMorphism)
    M = dsum([Morphism(m) for m ∈ f.m])
    return M.m
end

function (F::Field)(f::RingCatMorphism)
    if !(domain(f) == codomain(f) && issimple(domain(f)))
        throw(ErrorException("Cannot convert Morphism to $F"))
    end
    i = findfirst(e -> e == 1, domain(f).components)
    return F(f.m[i][1,1])
end

#-------------------------------------------------------------------------------
#   Tensor Product
#-------------------------------------------------------------------------------



function tensor_product(X::RingCatObject, Y::RingCatObject)
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

    return RingCatObject(C,T)
end

function tensor_product(f::RingCatMorphism, g::RingCatMorphism)
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
                    h = h ⊕ RingCatMorphism(simpl[k]^d1,simpl[k]^d2, m)
                end
            end
        end
    end
    #dom_left = dom.components - domain(h).components
    #cod_left = cod.components - codomain(h).components
    return h #⊕ zero_morphism(RingCatObject(C,dom_left), RingCatObject(C,cod_left))
end


function one(C::RingCategory) 
    if !isdefined(C, :one) 
        throw(ErrorException("There is no unit object defined"))
    end
    RingCatObject(C,C.one)
end
#-------------------------------------------------------------------------------
#   Direct sum
#-------------------------------------------------------------------------------

function dsum(X::RingCatObject, Y::RingCatObject, morphisms::Bool = false)
    @assert parent(X) == parent(Y) "Mismatching parents"
    if morphisms return dsum_with_morphisms(X,Y) end
    return RingCatObject(parent(X), X.components .+ Y.components)
end

function dsum_with_morphisms(X::RingCatObject, Y::RingCatObject)
    S = dsum(X,Y)
    ix_mats = matrices(zero_morphism(X,S))
    iy_mats = matrices(zero_morphism(Y,S))
    px_mats = matrices(zero_morphism(S,X))
    py_mats = matrices(zero_morphism(S,Y))

    for i ∈ 1:parent(X).simples
        (x,y) = X.components[i], Y.components[i]
        for j ∈ 1:x 
            ix_mats[i][j,j] = 1
            px_mats[i][j,j] = 1
        end
        for j ∈ 1:y 
            iy_mats[i][j,j+x] = 1
            py_mats[i][j+x,j] = 1
        end
    end

    ix = Morphism(X,S, ix_mats)
    px = Morphism(S,X, px_mats)
    iy = Morphism(Y,S, iy_mats)
    py = Morphism(S,Y, py_mats)

    return S,[ix,iy],[px,py]
end

function dsum(f::RingCatMorphism, g::RingCatMorphism)
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
    return RingCatMorphism(dom,cod, m)
end


zero(C::RingCategory) = RingCatObject(C,[0 for i ∈ 1:C.simples])

function zero_morphism(X::RingCatObject, Y::RingCatObject)
    return RingCatMorphism(X,Y,[zero(MatrixSpace(base_ring(X), cX, cY)) for (cX,cY) ∈ zip(X.components, Y.components)])
end

function isisomorphic(X::RingCatObject, Y::RingCatObject)
    if X != Y
        return false, nothing
    else
        return true, id(X)
    end
end
#-------------------------------------------------------------------------------
#   Simple Objects
#-------------------------------------------------------------------------------

function simples(C::RingCategory)
    n = C.simples
    [RingCatObject(C, [i == j ? 1 : 0 for j ∈ 1:n]) for i ∈ 1:n]
end

#-------------------------------------------------------------------------------
#   Kernel and Cokernel
#-------------------------------------------------------------------------------

function kernel(f::RingCatMorphism)
    C = parent(domain(f))
    kernels = [kernel(Morphism(m)) for m ∈ f.m]
    mats = [matrix(m) for (k,m) ∈ kernels]
    ker = RingCatObject(C,[dim(k) for (k,m) ∈ kernels])

    return ker, Morphism(ker, domain(f), mats)
end


function left_inverse(f::RingCatMorphism)
    inverses = [left_inverse(Morphism(m)) for m ∈ matrices(f)]
    mats = [matrix(m) for m ∈ inverses]
    return Morphism(codomain(f), domain(f), mats)
end

#-------------------------------------------------------------------------------
#   Examples
#-------------------------------------------------------------------------------

function Ising()
    Qx,x = QQ["x"]
    F,a = NumberField(x^2-2, "√2")
    C = RingCategory(F,["𝟙", "χ", "X"])
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

    set_one!(C,[1,0,0])

    set_spherical!(C, [F(1) for s ∈ simples(C)])

    a,b,c = simples(C)

    return C
end

#-------------------------------------------------------------------------------
#   Hom Spaces
#-------------------------------------------------------------------------------

struct RingCatHomSpace<: AbstractHomSpace
    X::RingCatObject
    Y::RingCatObject
    basis::Vector{RingCatMorphism}
    parent::VectorSpaces
end

function Hom(X::RingCatObject, Y::RingCatObject)
    @assert parent(X) == parent(Y) "Mismatching parents"
    Xi, Yi = X.components, Y.components
    F = base_ring(X)

    d = sum([x*y for (x,y) ∈ zip(Xi,Yi)])

    if d == 0 return RingCatHomSpace(X,Y,RingCatMorphism[], VectorSpaces(F)) end

    basis = [zero_morphism(X,Y).m for i ∈ 1:d]
    next = 1
    for k ∈ 1:parent(X).simples

        for i ∈ 1:Xi[k], j ∈ 1:Yi[k]
            basis[next][k][i,j] = 1
            next = next + 1
        end
    end
    basis_mors = [RingCatMorphism(X,Y,m) for m ∈ basis]
    return RingCatHomSpace(X,Y,basis_mors, VectorSpaces(F))
end

function express_in_basis(f::RingCatMorphism, base::Vector)
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

function show(io::IO, X::RingCatObject)
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

function show(io::IO, f::RingCatMorphism)
    print(io, """Morphism with
Domain: $(domain(f))
Codomain: $(codomain(f))
Matrices: """)
    print(io, join(["$(m)" for m ∈ f.m], ", "))
end

#-------------------------------------------------------------------------------
#   Utility
#-------------------------------------------------------------------------------
