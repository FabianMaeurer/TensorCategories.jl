#=----------------------------------------------------------
    Structures for ℤ₊-Rings  
----------------------------------------------------------=#

@attributes mutable struct ZPlusRing <: AbstractAssociativeAlgebra{ZZRingElem}
    algebra::AbstractAssociativeAlgebra{ZZRingElem}
    basis_names::Vector{String}

    function ZPlusRing(mult_table::Array{ZZRingElem,3}) 
        A = new()
        A.algebra = StructureConstantAlgebra(ZZ, mult_table)
        A.basis_names = ["x$i" for i ∈ 1:size(mult_table)[1]]
        return A
    end
    function ZPlusRing(mult_table::Array{ZZRingElem,3}, one::Vector{ZZRingElem}) 
        A = new()
        A.algebra = StructureConstantAlgebra(ZZ, mult_table,one)
        A.basis_names = ["x$i" for i ∈ 1:size(mult_table)[1]]
        return A
    end
    function ZPlusRing(A::AbstractAssociativeAlgebra)
        @assert base_ring(A) == ZZ
        B = new()
        B.algebra = A
        B.basis_names = ["x$i" for i ∈ 1:dim(A)[1]]
        return B
    end
end


function ZPlusRing(names::Vector{String}, mult_table::Array{<:Number,3}, one::Vector{<:Number})
    A = ZPlusRing(ZZ.(mult_table), ZZ.(one))
    A.basis_names = names
    return A
end

@alias ℤ₊Ring ZPlusRing
@alias ℕRing ZPlusRing

struct ZPlusRingElem <: AbstractAssociativeAlgebraElem{ZZRingElem}
    parent::ZPlusRing
    elem::AbstractAssociativeAlgebraElem{ZZRingElem}
end

@alias ℤ₊RingElem ZPlusRingElem
@alias ℕRingElem ZPlusRingElem

is_semisimple(R::ℕRing) = is_semisimple(R.algebra)

parent(x::ℕRingElem) = x.parent
base_ring(::ℕRing) = ZZ
base_ring(::ℕRingElem) = ZZ

StructureConstantAlgebra(R::ℕRing) = R.algebra
AssociativeAlgebraElem(x::ℕRingElem) = x.elem

has_one(R::ℕRing) = has_one(StructureConstantAlgebra(R))
one(R::ℕRing) = ℕRingElem(R, one(StructureConstantAlgebra(R)))
zero(R::ℕRing) = ℕRingElem(R, zero(StructureConstantAlgebra(R)))

basis(R::ℕRing) = [ℕRingElem(R, b) for b ∈ basis(StructureConstantAlgebra(R))]
getindex(R::ℕRing, k::Int) = basis(R)[k]

coefficients(x::ℕRingElem) = coefficients(AssociativeAlgebraElem(x))

(A::ℕRing)(c::Vector{ZZRingElem}) = ℕRingElem(A, AssociativeAlgebraElem{ZZRingElem, typeof(A.algebra)}(A.algebra, c))

multiplication_table(R::ℕRing) = R.algebra.mult_table
print_multiplication_table(R::ℕRing) = print_multiplication_table(multiplication_table(R))

function ==(x::ℕRingElem, y::ℕRingElem) 
    parent(x) == parent(y) && coefficients(x) == coefficients(y)
end

elem_type(T::Type{ZPlusRing}) = ZPlusRingElem

#=----------------------------------------------------------
    ℤ₊RingElem Operations 
----------------------------------------------------------=#

+(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), AssociativeAlgebraElem(x) +  AssociativeAlgebraElem(y))
-(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), AssociativeAlgebraElem(x) -  AssociativeAlgebraElem(y))
*(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), AssociativeAlgebraElem(x) *  AssociativeAlgebraElem(y))

div(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), div(AssociativeAlgebraElem(x) ,AssociativeAlgebraElem(y)))

*(λ::RingElem, x::ℕRingElem) = ℕRingElem(parent(x), λ*AssociativeAlgebraElem(x))
*(λ::Int, x::ℕRingElem) = ℕRingElem(parent(x), λ*AssociativeAlgebraElem(x))

function ^(x::ℕRingElem, k::Int) 
    if k == 0  
        return one(parent(x)) 
    elseif k == 1
        return x
    end
    *([x for i ∈ 1:k]...)
end

dim(A::ℕRing) = dim(A.algebra)

#=----------------------------------------------------------
    Based Rings 
----------------------------------------------------------=#

is_based(R::ℕRing) = get_attribute(R, :involution) !== nothing



function tau(x::ℕRingElem)
    sum(coefficients(x)[coefficients(one(parent(x))) .> 0])
end

@alias τ tau


function extension_of_scalars(A::ℕRing, K::Field)
    StructureConstantAlgebra(K,K.(A.algebra.mult_table), K.(A.algebra.one))
end

⊗(A::ℕRing, K::Field) = extension_of_scalars(A,K)


function fpdim(x::ℕRingElem)
    R = parent(x)
    K = base_ring(R)

    dims = get_attribute!(R, :fpdims) do 
        m = multiplication_table(R)
        [fp_eigenvalue(matrix(K,m[i,:,:])) for i ∈ 1:length(basis(R))]
    end

    sum(coefficients(x) .* dims)
end

#=----------------------------------------------------------
    Fusion Rings 
----------------------------------------------------------=#

function involution(x::ℕRingElem)
    if has_attribute(parent(x), :involution)
        invol = get_attribute(parent(x), :involution)
        S = basis(parent(x))
        return sum(coefficients(x) .* S[invol])
    end
    error("there is no involution on $(parent(x))")
end

#=----------------------------------------------------------
    Pretty Printing 
----------------------------------------------------------=#

function show(io::IO, A::ℕRing)
    if has_attribute(A, :name)
        print(io, "$(get_attribute(A,:name))")
        return
    end
    if has_attribute(A, :is_partial) 
        if get_attribute!(A, :is_partial)
            print(io, "ℤ₊-Ring with of unknown dimension")
            return 
        end
    end
    print(io, "ℤ₊-Ring with dimension $(length(coefficients(one(A))))")
end

function show(io::IO, x::ℕRingElem)
    if iszero(x.elem)
        print(io, "0")
        return
    end
    str = "" 
    coeffs = coefficients(x)
    for (c,n) ∈ zip(coeffs, parent(x).basis_names)
        if c != 0
            if c ∉ [1,-1]
                str = str*"$c⋅$n + "
            else
                str = str*n*" + "
            end
        end
    end

    print(io, str[1:end-2])
end


#=----------------------------------------------------------
    ℤ₊-Automorphisms         
----------------------------------------------------------=#

struct ZPlusRingMorphism
    domain::ZPlusRing
    codomain::ZPlusRing 
    images::Vector{ZPlusRingElem}
end

function hom(domain::ZPlusRing,
    codomain::ZPlusRing ,
    images::Vector{ZPlusRingElem};
    check = true)

    @req is_zplus_ring_morphism(basis,images) "Not a ℤ₊-Ring morphism"

    ZPlusRingMorphism(domain,codomain,basis(domain),images)
end

function is_zplus_ring_morphism(basis::Vector{ZPlusRingElem}, images::Vector{ZPlusRingElem})
    n = length(basis)
    for i ∈ 1:n, j ∈ 1:n 
        a,b = basis[[i,j]]
        ab = a*b
       
        image_ab = sum(c*v for (v,c) ∈ zip(images, coefficients(ab)))

        image_ab != images[i] * images[j] && return false
    end
    true
end

function automorphisms(R::ℕRing)
    B = basis(R) 
    G = symmetric_group(length(B))

    fpdims = fpdim.(B)
    perms = [g for g ∈ G if permuted(fpdims, g) == fpdims]

    automorphism_perms = [p for p ∈ perms if is_zplus_ring_morphism(B,permuted(B,p))]

    [ZPlusRingMorphism(R,R,permuted(B,p)) for p ∈ automorphism_perms]
end