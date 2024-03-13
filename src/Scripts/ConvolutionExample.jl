using Revise, JuCat, Oscar

G = symmetric_group(2)
F = GF(23)
Ω = [1,2,3]

function act(x,g)
    if g == one(parent(g)) return x end
    if x == 1 return x end
    if x == 2 return 3 end
    if x == 3 return 2 end
end

X = gset(G,act,Ω)

Conv = ConvolutionCategory(X,F)

e,triv,e21,e12,d,c = simples(Conv)

#str = print_multiplication_table([e,triv,e21,e12,d,c], ["ϵ", "triv", "e21", "e12", "d", "c"])
#---------------------------------------------------------------------
# Liams Example
#---------------------------------------------------------------------
using Revise,JuCat, Oscar
@elapsed begin
    n = 6
    G = symmetric_group(n)
    F = GF(23)
    Ω = [i for i in 1:n+1]
    X = gset(G,Ω)

    Conv = ConvolutionCategory(X,F)
    A,f = groethendieck_ring(Conv)
end
