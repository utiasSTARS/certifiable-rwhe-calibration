using LinearAlgebra

# Based on Andreff et al. 2001 but without the unknown scale (for now)
# We conduct a simultaneous (trans. + rot.) linear solve, then orthogonalize the Rotations

function linear_hand_eye_initialization(A, B, var_trans::Array, var_rot::Array)
    n = size(A, 1)
    M = zeros(12*n, 12)
    b = zeros(12*n)
    # Formulate the joint linear system
    for i in 1:n
        Mi, bi = linear_hand_eye_system(A[i], B[i], var_trans[i], var_rot[i])
        M[12*(i-1)+1:12*i, :] = Mi
        b[12*(i-1)+1:12*i] = bi
    end
    # Solve 
    x = M\b
    r = x[1:9]
    t = x[10:12]

    # Project (orthogonalize) the rotation matrix onto SO(3)
    R = nearest_rotation_matrix(reshape(r, (3, 3)))
    return R, t
end

function linear_hand_eye_initialization(A, B)
    n = size(A, 1)
    return linear_hand_eye_initialization(A, B, ones(n), ones(n))
end

function linear_hand_eye_system(A, B, var_trans, var_rot)
    Ra = A[1:3, 1:3]
    Rb = B[1:3, 1:3]
    ta = A[1:3, 4]
    tb = B[1:3, 4]
    M = [(kron(Rb, Ra)-I)/sqrt(var_rot) zeros(9, 3); 
        (-kron(tb', I(3)))/sqrt(var_trans) (Ra-I(3))/sqrt(var_trans)] 
    b = [zeros(9); -ta/sqrt(var_trans)]
    return M, b
end

# This is an old solution, but see Levinson et al. (2020) for details of special orthogonalization
function nearest_rotation_matrix(M)
    U, _, V = svd(M)
    S = diagm([1, 1, det(U*V')])
    return U*S*V'
end
