function *(m::MatElem{QQBarFieldElem}, n::MatElem{QQBarFieldElem})
    CC = CalciumField()
    m2 = change_base_ring(CC, m)
    n2 = change_base_ring(CC, n)
    change_base_ring(QQBarField(), m2 * n2)
end