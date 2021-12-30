function (F::FqNmodFiniteField)(x::GAP.FFE)
    if GAP.Globals.Characteristic(x) != Int(characteristic(F))
        throw(ErrorException("Mismatching characteristics"))
    end
    if x == GAP.Globals.Zero(x) return F(0) end
    if x == GAP.Globals.One(x) return F(1) end

    deg = GAP.Globals.DegreeFFE(x)
    char = Int(characteristic(F))
    exponent = GAP.Globals.LogFFE(x,GAP.Globals.Z(char,deg))
    a = gen(F)

    if exponent == GAP.Globals.fail throw(ErrorException("Conversion failed")) end

    return a^exponent
end
