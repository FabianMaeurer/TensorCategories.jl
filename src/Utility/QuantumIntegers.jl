#=----------------------------------------------------------
    Quantum Integers 
----------------------------------------------------------=#        

struct QuantumZZRing <: Ring
end

struct QuantumZZRingElem <: RingElem
    parent::QuantumZZRing
    n::Vector{Tuple{ZZRingElem,ZZRingElem}}
end

const QZZ = QuantumZZRing()

QuantumZZRingElem(n) = n == 0 ? QuantumZZRingElem(QZZ, Tuple{ZZRingElem,ZZRingElem}[]) : QuantumZZRingElem(QZZ, [(ZZ(1),ZZ(n))])

(R::QuantumZZRing)(n) = n == 0 ? QuantumZZRingElem(R, Tuple{ZZRingElem,ZZRingElem}[]) : QuantumZZRingElem(R, [(ZZ(1),n)])

zero(R::QuantumZZRing) = Quantum

#=----------------------------------------------------------
    Pretty print
----------------------------------------------------------=#

function show(io::IO, R::QuantumZZRing)
    print(io, """Quantum Integers""")
end

function show(io::IO, x::QuantumZZRingElem)
    n = length(x.n)
    if n == 0 print(io, "[0]"); return end
    s = ""
    for i ∈ 1:n
        k,m = x.n[i]
        if k == 1
            s = s*"[$m]"
        else
            s = s*"$k[$m]"
        end
        if i < n 
            s = s*" + "
        end
    end

    print(io,s)
end