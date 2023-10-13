#=----------------------------------------------------------
    The category Rep(M₁₁) over 𝔽₂ has four simple objects 
    of dimension 1,10,32 and 44. We consider the tensor
    power category of the 10-dimensional representation.
    It has seven indecomposable objects which we call
    a,b,c,d,e,f and g, a being the unit and b the
    generator. The computation is very expensive, thus 
    here a quick construction. 

    The minimal polynomials are 

        b⁵ - 12b⁴ + 19b³ + 12b² - 20b
        c⁵ - 94c⁴ - 116c³ + 1560c² - 1920c
----------------------------------------------------------=#

mult = zeros(Int,7,7,7)

# One is the unit
mult[1,:,:] = mult[:,1,:] = [i == j for i ∈ 1:7, j ∈ 1:7]

# Multiplication with b is obtained by decomposing the modules b⊗-
mult[2,2:7,:] = mult[2:7,2,:] = [0 1 1 0 0 0 0;
                                 0 2 0 2 5 1 0;
                                 0 0 0 1 2 0 0;
                                 0 0 0 3 8 2 0;
                                 0 0 0 1 5 1 1;
                                 0 0 0 2 6 2 0]

# Multiplication with c is obtained by calculations in the Ring
# E.g we have cx = (b² - b)x
mult[3,3:7,:] = mult[3:7,3,:] = [0 0 2 16 44 10 1;
                                 0 0 0  6 16  4 0;
                                 0 0 0 26 72 26 2;
                                 0 0 0 18 48 12 0;
                                 0 0 0 20 56 12 2;]

# Multiplication with d by 
# 

R = ℕRing(["a", "b", "c", "d", "e", "f", "g"], mult, [1,0,0,0,0,0,0])