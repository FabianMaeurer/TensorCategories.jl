#=----------------------------------------------------------
    G-Equivariantizations of Tensor Categories 
----------------------------------------------------------=#

mutable struct Equivariantization <: Category 
    category::Category 
    group::GAPGroup 
    gens::Vector{<:GroupElem}
    action::GTensorAction 

    simples::Vector{<:Object}

    function Equivariantization(category::Category, 
        group::GAPGroup, 
        gens::Vector{<:GroupElem},
        action::GTensorAction)

        new(category,group,gens,action)
    end
end 

function ==(C::Equivariantization, D::Equivariantization)
    category(C) == category(D) && gaction(C) == gaction(D)
end

struct EquivariantObject <: Object 
    parent::Equivariantization
    object::Object 
    structure_maps::Vector{<:Morphism} # One for each element of G
end

struct EquivariantMorphism <: Morphism 
    domain::EquivariantObject 
    codomain::EquivariantObject 
    morphism::Morphism 
end

function equivariantization(C::Category, G::GAPGroup, T::GTensorAction)
    Equivariantization(C,G,gens(G),T)
end

function equivariantization(C::Category, T::GTensorAction)
    equivariantization(C,group(T),T)    
end

gaction(E::Equivariantization) = E.action
morphism(X::EquivariantObject, Y::EquivariantObject, f::Morphism) = EquivariantMorphism(X,Y,f)

id(X::EquivariantObject) = morphism(X,X,id(object(X)))

is_multifusion(E::Equivariantization) = is_multifusion(category(E))
is_fusion(E::Equivariantization) = is_fusion(category(E))

is_multitesnsor(E::Equivariantization) = is_multitensor(category(E))
is_tensor(E::Equivariantization) = is_tensor(category(E))

group(E::Equivariantization) = E.group
#=----------------------------------------------------------
    direct_sum
----------------------------------------------------------=#

function direct_sum(X::EquivariantObject...)
    if length(X) == 1
        return X[1], id(X[1]), id(X[1])
    end

    S,incl,proj = direct_sum(object.(X)...)
    E = parent(X[1])
    
    T = gaction(E)
    G = group(T)

    structure_maps = [
        sum([i ∘ x.structure_maps[l] ∘ T(g)(p) for (i,p,x) ∈ zip(incl,proj,X)])
        for l ∈ 1:order(G)
    ]

    sum_in_E = EquivariantObject(E, S, structure_maps)
    incl_in_E = [morphism(x,E,i) for (x,i) ∈ zip(X,incl)]
    proj_in_E = [morphism(E,x,p) for (x,p) ∈ zip(X,proj)]

    return sum_in_E, incl_in_E, proj_in_E
end 

function direct_sum(f::EquivariantMorphism...)
    if length(f) == 1
        return f 
    end

    dom,_,_ = direct_sum(domain.(f)...)
    cod,_,_ = direct_sum(codomain.(f)...)

    return morphism(dom,cod, direct_sum(morphism.(f)...))
end

function zero(E::Equivariantization) 
    Z = zero(category(E))
    EquivariantObject(E, Z, [id(Z) for g ∈ group(E)])
end

#=----------------------------------------------------------
    Tensor Product 
----------------------------------------------------------=#    

function tensor_product(X::EquivariantObject, Y::EquivariantObject)
    @assert parent(X) == parent(Y)

    Z = object(X) ⊗ object(Y)
    E = parent(X) 

    T = gaction(E)

    structure_maps = [(u ⊗ v) ∘ monoidal_structure(Tg,object(X),object(Y)) for (u,v,Tg) ∈ zip(X.structure_maps,Y.structure_maps, images(T))]

    return EquivariantObject(E, Z, structure_maps)
end

function tensor_product(f::EquivariantMorphism, g::EquivariantMorphism)
    dom = domain(f) ⊗ domain(g) 
    cod = codomain(f) ⊗ codomain(g)

    morphism(dom, cod, morphism(f) ⊗ morphism(g))
end


function associator(X::EquivariantObject, Y::EquivariantObject, Z::EquivariantObject)
    dom = tensor_product([X,Y,Z])
    cod = X ⊗ (Y ⊗ Z)

    return morphism(dom, cod, associator(object.([X,Y,Z])...))
end

function one(E::Equivariantization)
    O = one(category(E))
    EquivariantObject(E, O, [id(O) for g ∈ group(E)])
end


#=----------------------------------------------------------
    Kernel & Cokernel       
----------------------------------------------------------=#

function kernel(f::EquivariantMorphism)
    
    K,k = kernel(morphism(f))

    inv_k = left_inverse(k)
    T = gaction(parent(f))
    G = group(T)

    structure_maps = [inv_k ∘ u ∘ T(g)(k) for (g,u) ∈ zip(elements(G), domain(f).structure_maps)]

    EK = EquivariantObject(parent(f), K, structure_maps)

    EK, EquivariantMorphism(EK,domain(f),k)
end


function cokernel(f::EquivariantMorphism)
    
    C,c = cokernel(morphism(f))

    inv_c = right_inverse(c)
    T = gaction(parent(f))
    G = group(T)

    structure_maps = [c ∘ u ∘ T(g)(inv_c) for (g,u) ∈ zip(elements(G), codomain(f).structure_maps)]

    EC = EquivariantObject(parent(f), C, structure_maps)

    EC, EquivariantMorphism(domain(f),EC,c)
end


#=----------------------------------------------------------
    Hom Spaces 
----------------------------------------------------------=#

function Hom(X::EquivariantObject, Y::EquivariantObject)
    H = basis(Hom(object(X), object(Y)))
    N = length(H)

    K = base_ring(X) 

    T = gaction(parent(X))
    G = group(T)

    Kx, x = polynomial_ring(K, N)
    equations = elem_type(Kx)[]

    for (g,u,v) ∈ zip(elements(G), X.structure_maps, Y.structure_maps)
        base = basis(Hom(T(g)(object(X)),object(Y)))
        eq = zeros(Kx, length(base))
        Tg = T(g)
        for (f,a) ∈ zip(H,x)
            coeffs = express_in_basis(v ∘ Tg(f), base)
            eq = eq .+ (a .* coeffs)
        end
        for (f,a) ∈ zip(H,x)
            coeffs = express_in_basis(f ∘ u, base)
            eq = eq .- (a .* coeffs)
        end
        equations = [equations; eq]
    end

    M = zero(matrix_space(K,length(equations),N))

    for (i,e) ∈ zip(1:length(equations),equations)
        M[i,:] = [coeff(e, a) for a ∈ x]
    end

    N = nullspace(M)[2]

    _,cols = size(N)

    basis_coeffs = [collect(N[:,i]) for i ∈ 1:cols]

    equi_basis = [EquivariantMorphism(X,Y,sum(b .* H)) for b ∈ basis_coeffs]

    return HomSpace(X,Y,equi_basis)

end

#=----------------------------------------------------------
    simples 
----------------------------------------------------------=#

function simples(E::Equivariantization)

    isdefined(E,:simples) && return E.simples

    S = simples(category(E))
    T = gaction(E)

    inductions = [equivariant_induction(s,T,E) for s ∈ S]
    ends = [equivariant_induction_adjunction(s,i,T,i) for (s,i) ∈ zip(S,inductions)]

    simpls = vcat([simple_subobjects(i,H) for (i,H) ∈ zip(inductions, ends)]...)

    simpls = unique_simples(simpls)

    E.simples = simpls 
end

#=----------------------------------------------------------
    Pretty Printing 
----------------------------------------------------------=#

function show(io::IO, X::EquivariantObject)
    print(io, "Equivariant Object: $(object(X))")
end

function show(io::IO, E::Equivariantization)
    print(io, "Equivariantization of $(category(E)) by $(group(gaction(E)))")
end

function show(io::IO, f::EquivariantMorphism)
    print(io, "EquivariantMorphism: $(morphism(f))")
end

