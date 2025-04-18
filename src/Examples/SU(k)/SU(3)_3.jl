#=----------------------------------------------------------
    Category su(3)₃/Z₃ as computed in 
    https://iopscience.iop.org/article/10.1088/1751-8113/43/39/395205/pdf
----------------------------------------------------------=#

function su_3_3_subcategory(K::Field = cyclotomic_field(12)[1])

    r3 = sqrt(K(3))
    r12 = 2*r3

    fusion_rules = zeros(Int, 4,4,4)

    for i ∈ 1:4
        fusion_rules[1,i,:] = fusion_rules[i,1,:] = [Int(j == i) for j ∈ 1:4]
    end

    fusion_rules[2,2,:] = [1,2,1,1]
    fusion_rules[2,3,:] = fusion_rules[3,2,:] = [0,1,0,0]
    fusion_rules[2,4,:] = fusion_rules[4,2,:] = [0,1,0,0]
    fusion_rules[3,3,:] = [0,0,0,1]
    fusion_rules[3,4,:] = fusion_rules[4,3,:] = [1,0,0,0]
    fusion_rules[4,4,:] = [0,0,1,0]

    C = six_j_category(K, fusion_rules, ["1", "8", "10", "1̅0̅"])

    set_one!(C, [1,0,0,0])

    # Set the associators 

    # F-symbol -1
    minus_one = [(2,2,3,3), (2,2,3,4), (2,2,4,3), (2,2,4,4), (2,3,3,2), (2,3,4,2), (2,4,3,2), (2,4,4,2), (3,3,2,2), (3,4,2,2), (4,3,2,2), (4,4,2,2), (3,2,2,3), (3,2,2,4), (4,2,2,3), (4,2,2,4), (3,3,4,3), (3,4,4,4), (4,3,3,3), (4,4,3,4)]

    for v ∈ minus_one
        set_associator!(C, v..., matrix(K,1,1,[-1]))
    end

    # Higher dimensions 

    m1 = matrix(K, 2,2, [-1//2 -r3//2; r3//2 -1//2])
    m2 = matrix(K, 7,7, [
        1//3 1//r3 0 0 1//r3 -1//3 -1//3;
        1//r3 -1//2 0 0 1//2 1//r12 1//r12;
        0 0 1//2 1//2 0 1//2 -1//2;
        0 0 1//2 1//2 0 -1//2 1//2;
        1//r3 1//2 0 0 -1//2 1//r12 1//r12;
        -1//3 1//r12 -1//2 1//2 1//r12 1//3 1//3;
        -1//3 1//r12 1//2 -1//2 1//r12 1//3 1//3;
    ])

    m1_index = [(2,2,2,3), (2,2,3,2), (2,4,2,2), (3,2,2,2)]
    m1_transposed = [(2,2,2,4), (2,2,4,2), (2,3,2,2), (4,2,2,2)]

    for v ∈ m1_index
        set_associator!(C, v..., m1)
    end
    for v ∈ m1_transposed
        set_associator!(C, v..., transpose(m1))
    end

    set_associator!(C, 2,2,2,2, m2)

    # Braiding 

    R = Array{MatElem}(undef, 4,4,4)

    for i ∈ 1:4 
        R[1,i,:] = R[i,1,:] = [j == i ? identity_matrix(K,1) : zero_matrix(K,0,0) for j ∈ 1:4] 
    end

    i = root_of_unity(K, 4) 

    R[2,2,:] = [
        -identity_matrix(K,1), 
        matrix(K,2,2, [i 0; 0 -i]),
        -identity_matrix(K,1),
        -identity_matrix(K,1)] 

    R[2,3,:] = R[3,2,:] = [j == 2 ? identity_matrix(K,1) : zero_matrix(K,0,0) for j ∈ 1:4]
    R[3,3,:] = [j == 4 ? identity_matrix(K,1) : zero_matrix(K,0,0) for j ∈ 1:4]
    R[4,4,:] = [j == 3 ? identity_matrix(K,1) : zero_matrix(K,0,0) for j ∈ 1:4]

    R[3,4,:] = R[4,3,:] = [j == 1 ? identity_matrix(K,1) : zero_matrix(K,0,0) for j ∈ 1:4]
    R[2,4,:] = R[4,2,:] = [j == 2 ? identity_matrix(K,1) : zero_matrix(K,0,0) for j ∈ 1:4]

    set_braiding!(C,R)

    return C
end