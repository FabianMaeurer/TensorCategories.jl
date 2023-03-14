function (F::FinField)(x::GAP.FFE)
    if GAP.Globals.Characteristic(x) != Int(characteristic(F))
        throw(ErrorException("Mismatching characteristics"))
    end
    if x == GAP.Globals.Zero(x) return F(0) end
    if x == GAP.Globals.One(x) return F(1) end

    deg = degree(F)

    if deg == 1 return F(GAP.Globals.IntFFE(x)) end

    char = Int(characteristic(F))
    exponent = GAP.Globals.LogFFE(x,GAP.Globals.Z(char,deg))

    if exponent == GAP.Globals.fail throw(ErrorException("Conversion failed")) end

    a = gen(F)
    return a^exponent
end
