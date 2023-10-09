#=----------------------------------------------------------
    Structures for ℤ₊-Rings  
----------------------------------------------------------=#

@attributes mutable struct ZPlusRing <: AbsAlgAss{fmpz}
    algebra::AbsAlgAss{fmpz}
    basis_names::Vector{String}

    function ZPlusRing(mult_table::Array{fmpz,3}) 
        A = new()
        A.algebra = AlgAss(ZZ, mult_table)
        A.basis_names = ["x$i" for i ∈ 1:size(mult_table)[1]]
        return A
    end
    function ZPlusRing(mult_table::Array{fmpz,3}, one::Vector{fmpz}) 
        A = new()
        A.algebra = AlgAss(ZZ, mult_table,one)
        A.basis_names = ["x$i" for i ∈ 1:size(mult_table)[1]]
        return A
    end
    function ZPlusRing(A::AbsAlgAss)
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

struct ZPlusRingElem <: AbsAlgAssElem{fmpz}
    parent::ZPlusRing
    elem::AbsAlgAssElem{fmpz}
end

@alias ℤ₊RingElem ZPlusRingElem
@alias ℕRingElem ZPlusRingElem

is_semisimple(R::ℕRing) = is_semisimple(R.algebra)

parent(x::ℕRingElem) = x.parent
base_ring(::ℕRing) = ZZ
base_ring(::ℕRingElem) = ZZ

AlgAss(R::ℕRing) = R.algebra
AlgAssElem(x::ℕRingElem) = x.elem

has_one(R::ℕRing) = has_one(AlgAss(R))
one(R::ℕRing) = ℕRingElem(R, one(AlgAss(R)))
zero(R::ℕRing) = ℕRingElem(R, zero(AlgAss(R)))

basis(R::ℕRing) = [ℕRingElem(R, b) for b ∈ basis(AlgAss(R))]
getindex(R::ℕRing, k::Int) = basis(R)[k]

coefficients(x::ℕRingElem) = coefficients(AlgAssElem(x))

(A::ℕRing)(c::Vector{fmpz}) = ℕRingElem(A, AlgAssElem{ZZRingElem, typeof(A.algebra)}(A.algebra, c))

multiplication_table(R::ℕRing) = R.algebra.mult_table
print_multiplication_table(R::ℕRing) = print_multiplication_table(multiplication_table(R))

#=----------------------------------------------------------
    ℤ₊RingElem Operations 
----------------------------------------------------------=#

+(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), AlgAssElem(x) +  AlgAssElem(y))
-(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), AlgAssElem(x) -  AlgAssElem(y))
*(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), AlgAssElem(x) *  AlgAssElem(y))

div(x::ℕRingElem, y::ℕRingElem) = ℕRingElem(parent(x), div(AlgAssElem(x) ,AlgAssElem(y)))

*(λ::RingElem, x::ℕRingElem) = ℕRingElem(parent(x), λ*AlgAssElem(x))



#=----------------------------------------------------------
    Based Rings 
----------------------------------------------------------=#

is_based(R::ℕRing) = get_attribute(R, :involution) !== nothing



function tau(x::ℕRingElem)
    sum(coefficients(x)[coefficients(one(parent(x))) .> 0])
end

@alias τ tau


function extension_of_scalars(A::ℕRing, K::Field)
    AlgAss(K,K.(A.algebra.mult_table), K.(A.algebra.one))
end

⊗(A::ℕRing, K::Field) = extension_of_scalars(A,K)

#=----------------------------------------------------------
    Fusion Rings 
----------------------------------------------------------=#

function involution(x::ℕRingElem)
    if has_attribute(parent(x), :involution)
        invol = get_attribute(parent(x), :involution)
        return sum(coefficients(x) .* invol)
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
    str = "" 
    coeffs = coefficients(x)
    for (c,n) ∈ zip(coeffs, parent(x).basis_names)
        if c > 0
            if c > 1
                str = str*"$c⋅$n + "
            else
                str = str*n*" + "
            end
        end
    end

    print(io, str[1:end-2])
end
