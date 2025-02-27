function extended_haagerup(K::Ring)
    H = six_j_category(K, ["1", "f⁽²⁾", "f⁽⁴⁾", "f⁽⁶⁾", "P'", "Q'"])

    mult = zeros(Int,6,6,6)

    mult[1,:,:] = [i == j for i in 1:6, j in 1:6]
    mult[2,:,:] = [0 1 0 0 0 0;
        1 1 1 0 0 0;
        0 1 1 1 0 0;
        0 0 1 1 1 1;
        0 0 0 1 2 1;
        0 0 0 1 1 0]
    mult[3,:,:] = [ 0 0 1 0 0 0;
        0 1 1 1 0 0;
        1 1 1 1 1 1;
        0 1 1 2 3 1;
        0 0 1 3 3 2;
        0 0 1 1 2 1]
    mult[4,:,:] = [0 0 0 1 0 0;
        0 0 1 1 1 1;
        0 1 1 2 3 1;
        1 1 2 4 5 3;
        0 1 3 5 6 3;
        0 1 1 3 3 2]
    mult[5,:,:] = [0 0 0 0 1 0;
        0 0 0 1 2 1;
        0 0 1 3 3 2;
        0 1 3 5 6 3;
        1 2 3 6 7 4;
        0 1 2 3 4 2]
    mult[6,:,:] = [0 0 0 0 0 1;
        0 0 0 1 1 0;
        0 0 1 1 2 1;
        0 1 1 3 3 2;
        0 1 2 3 4 2;
        1 0 1 2 2 1]

    set_tensor_product!(H, mult)
    set_one!(H, [1,0,0,0,0,0])

    return H
end