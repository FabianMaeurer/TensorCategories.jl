#=----------------------------------------------------------
    Make MacLanes strictness theorem explicit for fusion 
    categories.

    Every 𝒞-module endofunctor of 𝒞 is given by the tuple
    (X⊗-, (a_{X,-,-})). If we assume 𝒞 to be fusion we can 
    reduce to a_{X,Xᵢ,Xⱼ} for simple objects Xᵢ
----------------------------------------------------------=#

struct FusionModuleEndoFunctors :> Category
    category::Category
end

struct FusionModuleEndoFunctor <: Object 
    parent::FusionModuleEndoFunctors
    object::Object
    module_structure::Array{<:Morphism,2}
end

struct FusionModuleTransformation <: Morphism
    domain::FusionModuleEndoFunctor
    codomain::FusionModuleEndoFunctor
    m::Morphism
end

is_strict(C::FusionModuleEndoFunctors) = true

object(F::FusionModuleEndoFunctor) = F.object
module_structure(F::FusionModuleEndoFunctor) = F.module_structure

morphism(f::FusionModuleTransformation) = f.m

morphism(X::FusionModuleEndoFunctor, Y::FusionModuleEndoFunctor, f::Morphism) = FusionModuleTransformation(X,Y,f)
#=----------------------------------------------------------
    Functionality 
----------------------------------------------------------=#    

function compose(f::FusionModuleTransformation...) 
    dom = domain(f[1])
    cod = codomain(f[end])
    return FusionModuleTransformation(dom, cod, compose(morphism.(f)))
end


function FusionModuleEndoFunctor(X::Object)
    C = parent(X)
    @assert is_fusion(C)

    mod_structure = [associator(X,i,j) for i ∈ simples(C), j ∈ simples(C)]

    return FusionModuleEndoFunctor(C,X,mod_structure)
end

matrix(f::FusionModuleTransformation) = matrix(morphism(f))

#=----------------------------------------------------------
    Abelian structure 
----------------------------------------------------------=#

function direct_sum(X::FusionModuleEndoFunctor...)
    S,incl,proj = direct_sum(object.(X)...)

    S = FusionModuleEndoFunctor(parent(X[1]), S,
            [direct_sum(f...) for f ∈ zip(module_structure.(X)...)])

    incl = [FusionModuleTransformation(x,S,i) for (x,i) ∈ zip(X,incl)]
    proj = [FusionModuleTransformation(S,x,p) for (x,p) ∈ zip(X,proj)]

    return S, incl, proj
end

function direct_sum(f::FusionModuleTransformation...)
    dom,_,_ = direct_sum(domain.(f)...)
    cod,_,_ = direct_sum(codomain.(f)...)

    return FusionModuleTransformation(dom, cod, direct_sum(morphism.(f)...))
end

+(f::FusionModuleTransformation...) = morphism(domain(f[1]), domain(f[1]), +(morphism.(f)...))

*(λ::FieldElem, f::FusionModuleTransformation) = morphism(domain(f), codomain(f), λ * morphism(f))

function kernel(f::FusionModuleTransformation)
    K,k = kernel(morphism(f))
    K = FusionModuleEndoFunctor(K)
    return K, morphism(K, domain(f), k)
end

function cokernel(f::FusionModuleTransformation)
    C,c = cokernel(morphism(f))
    C = FusionModuleEndoFunctor(C)
    return C, morphism(codomain(f), C, c)
end

function zero(C::FusionModuleEndoFunctors)
    FusionModuleEndoFunctor(zero(category(C)))
end

#=----------------------------------------------------------
    Tensor structure 
----------------------------------------------------------=#

function tensor_product(X::FusionModuleEndoFunctor, Y::FusionModuleEndoFunctor)
    C = parent(X)
    S = simples(category(C))
    x = object(X)
    y = object(Y)

    mors = []