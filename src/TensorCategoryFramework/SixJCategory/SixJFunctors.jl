#=----------------------------------------------------------
    Functors for 6j-Categories building on fusion ring
    morphisms 
----------------------------------------------------------=#

struct SixJFunctor <: AbstractMonoidalFunctor 
    domain::SixJCategory
    codomain::SixJCategory 
    images::Vector{SixJObject}
    monoidal_structure::Dict
end

function (F::SixJFunctor)(X::SixJObject)
    n = parent(X).simples
    direct_sum([F.images[i]^X.components[i] for i ∈ 1:n]...)[1]
end

function (F::SixJFunctor)(f::SixJMorphism)

    dom = domain(f)
    cod = codomain(f)

    C = domain(F)

    dom_dec = vcat([[C[i] for _ ∈ 1:c] for (i,c) ∈ zip(1:C.simples, dom.components)]...)
    cod_dec = vcat([[C[i] for _ ∈ 1:c] for (i,c) ∈ zip(1:C.simples, cod.components)]...)

    _, incl, _ = direct_sum(dom_dec)
    _, _, proj = direct_sum(cod_dec)

    Fdom, _, Fproj = direct_sum(F.(dom_dec))
    Fcod, Fincl, _ = direct_sum(F.(cod_dec))

    K = base_ring(C)

    f_components = [Fi ∘ (K(p ∘ f ∘ i) * id(F(domain(i)))) ∘ Fp for (i,Fp) ∈ zip(incl, Fproj), (Fi, p) ∈  zip(Fincl,proj) if domain(i) == codomain(p)][:]

    return sum(f_components[:])
end

indecomposables(F::SixJFunctor) = simples(domain(F))


function compose(F::SixJFunctor, G::SixJFunctor)
    images = G.(F.images)
    S = indecomposables(G)

    monoidal = Dict((X,Y) => G(monoidal_structure(F, X,Y)) ∘ monoidal_structure(G,F(X),F(Y)) for X ∈ S, Y ∈ S)

    SixJFunctor(domain(F), codomain(G), images, monoidal)
end