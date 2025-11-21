"""
Quadratic cost functions resulting from AX=XB problems.

For the interpretation in Giamou et al. 2019 and Wise et al. 2020:
    - A is a concatenation (along the first dimension) of R_b_i or T_b_i, which are the rotation/translation of 
    sensor b in its own frame (egomotion measurements)
    - B is a concatenation (along the first dimension) of R_a_i or T_a_i, which are the rotation/translation of 
    sensor a in its own frame (egomotion measurements)

The flip between the indices (_a and _b) and the definition of A and B is because this code uses the common AX=XB convention.
This is more than a bit confusing, we should have switched the definitions of a and b in that paper to match 
the majority of the AX=XB or hand-eye calibration literature.

The variable ordering in the combined vector is: x = [t α r s]', where r = vec(R). 

"""

# Create a cost matrix for a set of observed SO(d) rotations
function rotation_hand_eye_cost(A::AbstractArray, B::AbstractArray, κ::Array)
    n = length(κ)
    dim = size(A, 2)
    M = zeros(dim^2 + 1, dim^2 + 1)
    for i in 1:n
        M_i = kron(I(dim), A[i, :, :]) - kron(B[i, :, :]', I(dim))
        M[1:end-1, 1:end-1] += κ[i]*M_i'*M_i
    end
    return M
end

function rotation_hand_eye_cost(A::AbstractArray, B::AbstractArray)
    return rotation_hand_eye_cost(A, B, ones(size(A, 1)))
end

function transformation_hand_eye_cost(A::AbstractArray, B::AbstractArray)
    return transformation_hand_eye_cost(A, B, ones(size(A, 1)), ones(size(A, 1)))
end

function transformation_hand_eye_cost_scale(A::AbstractArray, B::AbstractArray)
    return transformation_hand_eye_cost(A, B, true, ones(size(A, 1)), ones(size(A, 1)))
end

function translation_from_rotation(R::Matrix, Q::Matrix)
    r = [reshape(R, 9); 1]
    return -inv(Q[1:3, 1:3])*Q[1:3, 4:end]*r
end

function translation_and_scale_from_rotation(R::Matrix, Q::Matrix)
    r = [reshape(R, 9); 1]
    t_scale = -inv(Q[1:4, 1:4])*Q[1:4, 5:end]*r
    return t_scale[1:3], t_scale[4]
end

function translation_cost(A::AbstractArray, B::AbstractArray, τ::Array)
    n = length(τ)
    dim = size(A, 2) - 1
    Q_t = zeros(dim^2 + dim + 1, dim^2 + dim + 1)
    for i in 1:n
        M_i = [(I(dim)-A[i,1:3, 1:3]) kron(B[i, 1:3, 4]', I(dim)) (-A[i,1:3, 4])]
        Q_t += τ[i]*M_i'*M_i
    end
    return Q_t
end

function translation_cost_scale(A::AbstractArray, B::AbstractArray, τ::Array)
    n = length(τ)
    dim = size(A, 2) - 1
    Q_t = zeros(dim^2 + dim + 2, dim^2 + dim + 2)
    for i in 1:n
        M_i = [(I(dim)-A[i,1:3, 1:3]) (-A[i,1:3, 4]) kron(B[i, 1:3, 4]', I(dim)) zeros(3)]
        Q_t += τ[i]*M_i'*M_i
    end
    return Q_t
end

# Create a cost matrix for a set of observed Special Euclidean (SE(d)) transformations, where A has unknown scale α > 0
function transformation_hand_eye_cost(A::AbstractArray, B::AbstractArray, scale::Bool, κ::Array, τ::Array)
    n = length(κ)
    dim = size(A, 2) - 1

    # Form the translational cost component
    if scale
        Q_t = translation_cost_scale(A, B, τ)
    else
        Q_t = translation_cost(A, B, τ)
    end

    # Get the rotation term 
    Q_r = rotation_hand_eye_cost(A[:, 1:dim, 1:dim], B[:, 1:dim, 1:dim], κ)

    # Combine the rotation term with the translation term
    Q = [zeros(dim+scale, size(Q_t, 1));
         zeros(10, dim+scale) Q_r] + Q_t

    # Schur complement for rotation-only component
    Q_cost = Q[4+scale:end, 4+scale:end] - Q[1:3+scale, 4+scale:end]'*inv(Q[1:3+scale, 1:3+scale])*Q[1:3+scale, 4+scale:end]
    return  Q_cost, Q
end

function transformation_hand_eye_cost(A::AbstractArray, B::AbstractArray, κ::Array, τ::Array)
    return transformation_hand_eye_cost(A, B, false, κ, τ)
end

function transformation_hand_eye_cost_scale(A::AbstractArray, B::AbstractArray, κ::Array, τ::Array)
    return transformation_hand_eye_cost(A, B, true, κ, τ)
end