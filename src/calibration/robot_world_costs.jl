using LinearAlgebra

"""
Quadratic cost functions resulting from AX=YB problems.

The variable ordering in the combined vector is: x = [t_x t_y α r_x r_y s]', where r_z = vec(R_z), z ∈ {x, y}. 
"""

# Create a cost matrix for a a set of observed SO(d) rotations
function rotation_robot_world_cost(A::AbstractArray, B::AbstractArray, κ::Array)
    n = length(κ)
    dim = size(A, 2)
    M = zeros(2*dim^2 + 1, 2*dim^2 + 1)
    #Q = diagm([sum(κ)*ones(2*dim^2); 0])
    #for i in 1:dim
    #    for j in 1:dim
    #        Q[dim^2 + 3 * i - 2:dim^2 + 3*i, 3 * j - 2:3 * j] .= -reshape(sum(κ .* B[:, i, j] .*A, dims=1), (3,3))
    #    end
    #end
    #Q[1:dim^2, dim^2 + 1:end-1] .= Q[dim^2 + 1:end-1, 1:dim^2]'
    for i in 1:n
        M_i = [kron(B[i, :, :], A[i, :, :]) -I(dim^2)]
        M[1:end-1, 1:end-1] += κ[i]*M_i'*M_i
    end
    return M
end

function rotation_robot_world_cost(A::AbstractArray, B::AbstractArray)
    return rotation_robot_world_cost(A, B, ones(size(A, 1)))
end

function transformation_robot_world_cost(A::AbstractArray, B::AbstractArray)
    return transformation_robot_world_cost(A, B, ones(size(A, 1)), ones(size(A, 1)))
end

function transformation_robot_world_cost_scale(A::AbstractArray, B::AbstractArray)
    return transformation_robot_world_cost(A, B, true, ones(size(A, 1)), ones(size(A, 1)))
end

function translations_from_rotations(R_x::Matrix, R_y::Matrix, Q::Matrix)
    r_x = reshape(R_x, 9)
    r_y = [reshape(R_y, 9); 1]
    r = [r_x; r_y]
    return -inv(Q[1:6, 1:6])*Q[1:6, 7:end]*r
end

function translation_cost(A::AbstractArray, B::AbstractArray, τ::Array)
    n = length(τ)
    dim = size(A, 2) - 1
    Q_t = zeros(2*dim^2 + 2*dim + 1, 2*dim^2 + 2*dim + 1)
    for i in 1:n
        M_i = [A[i, 1:3, 1:3] -I(3) zeros(3, 9) -kron(B[i, 1:3, 4]', I(dim)) A[i, 1:3, 4]]
        Q_t += τ[i]*M_i'*M_i
    end
    return Q_t
end

function translation_cost_scale(A::AbstractArray, B::AbstractArray, τ::Array)
    n = length(τ)
    dim = size(A, 2) - 1
    Q_t = zeros(2*dim^2 + 2*dim + 2, 2*dim^2 + 2*dim + 2)
    for i in 1:n
        M_i = [A[i, 1:3, 1:3] -I(3) A[i, 1:3, 4] zeros(3, 9) -kron(B[i, 1:3, 4]', I(dim)) zeros(3)]
        Q_t += τ[i]*M_i'*M_i
    end
    return Q_t
end

# Create a cost matrix for a set of observed Special Euclidean (SE(d)) transformations, where A has unknown scale α > 0
function transformation_robot_world_cost(A::AbstractArray, B::AbstractArray, scale::Bool, κ::Array, τ::Array)
    n = length(κ)
    dim = size(A, 2) - 1

    # Form the translational cost component
    if scale
        Q_t = translation_cost_scale(A, B, τ)
    else
        Q_t = translation_cost(A, B, τ)
    end

    # Get the rotation term 
    Q_r = rotation_robot_world_cost(A[:, 1:dim, 1:dim], B[:, 1:dim, 1:dim], κ)

    # Combine the rotation term with the translation term
    Q = [zeros(2*dim+scale, size(Q_t, 1));
         zeros(2*dim^2+1, 2*dim+scale) Q_r] + Q_t

    # Schur complement for rotation-only component
    Q_cost = Q[7+scale:end, 7+scale:end] - Q[1:6+scale, 7+scale:end]'*inv(Q[1:6+scale, 1:6+scale])*Q[1:6+scale, 7+scale:end]
    return  Q_cost, Q
end

function transformation_robot_world_cost(A::AbstractArray, B::AbstractArray, κ::Array, τ::Array)
    return transformation_robot_world_cost(A, B, false, κ, τ)
end

function transformation_robot_world_cost_scale(A::AbstractArray, B::AbstractArray, κ::Array, τ::Array)
    return transformation_robot_world_cost(A, B, true, κ, τ)
end

function sparse_robot_world_cost(A::AbstractArray, B::AbstractArray, scale::Bool, κ::Array, τ::Array)
    Qr_coo = sparse_rotation_cost(A, B, scale, κ)
    Qt_coo = sparse_translation_cost(A, B, scale, τ)
    if scale
        Qr_coo = [Qr_coo; [26 26 0]]
        Qt_coo = [Qt_coo; [26 26 0]]
    end
    Q_sparse= sparse(Qr_coo[:, 1], Qr_coo[:, 2], Qr_coo[:, 3]) + sparse(Qt_coo[:, 1], Qt_coo[:, 2], Qt_coo[:, 3])
    return Q_sparse
end

function sparse_rotation_cost(A::AbstractArray, B::AbstractArray, scale::Bool, κ::Array)
    #Compute the sparse rotation cost for the robot world problem
    dim = size(A, 2)
    #Add the Diagonal Cost Values
    Qr_coo = zeros(2*dim^2+2*dim^4, 3)
    Qr_coo[1:2*dim^2, :] = get_identity_costs(range(1, 2*dim^2, step=1)|>collect, κ)
    #Qr_coo[1:2*dim^2, :] = [range(1, 2*dim^2, step=1)|>collect range(1, 2*dim^2, step=1)|>collect sum(κ)*ones(2*dim^2)]

    #Compute the dense off diagonal values
    Q_off = zeros(dim^2, dim^2)
    for pair in zip(κ, axes(A, 1)|>collect)
        Q_off .-= pair[1] * kron(B[pair[2], 1:3, 1:3], A[pair[2], 1:3, 1:3])
    end

    #Get the row cols for each entry of Q_off
    indices = range(1, dim^2, step=1) |>collect #Assume that A[i, :, :] is  a square matrix
    row_col = [vec(repeat(indices, 1, length(indices))).+9 vec(repeat(indices', length(indices), 1))]

    #Add dense off diagonal values to COO representation
    Qr_coo[2*dim^2 + 1:end, :] = [[row_col vec(Q_off)]; [row_col[:, 2] row_col[:, 1] vec(Q_off)]]

    #Shift to correct position
    if scale
        Qr_coo[:, 1] .+= 7
        Qr_coo[:, 2] .+= 7
    else
        Qr_coo[:, 1] .+= 6
        Qr_coo[:, 2] .+= 6
    end
    return Qr_coo
end

function get_identity_costs(indices::Array, weights::Array)
    return [indices indices sum(weights)*ones(length(indices))]
end

function sparse_translation_cost(A::AbstractArray, B::AbstractArray, scale::Bool, τ::Array)
    #Compute the sparse rotation cost for the robot world problem
    dim = size(A, 2) - 1

    Qt_coo = zeros(154, 3)

    #Add the Identity Cost Values
    #Qt_coo[1:6, :] = [range(1, 2*dim, step=1)|>collect range(1, 2*dim, step=1)|>collect sum(τ)*ones(2*dim)]
    Qt_coo[1:6, :] = get_identity_costs(range(1, 2*dim, step=1)|>collect, τ)

    #Add dense and sparse entries to Qt_coo
    Qt_coo[7:109, :] = get_dense_terms_sparse_translation_cost(A, B, τ)
    Qt_coo[110:154, :] = get_sparse_terms_sparse_translation_cost(B, τ)
    
    #Convert the scale matrix to the normal one
    if ~scale
        Qt_coo = map_scale_to_scaled(Qt_coo)
    end
    return Qt_coo
end

function get_dense_terms_sparse_translation_cost(A::AbstractArray, B::AbstractArray, τ::Array)
    #Initialize temporary storage
    Q_coo = [get_row_cols_dense_terms_sparse_translation_cost() zeros(103)]

    #Compute dense entries
    for pair in zip(τ, axes(A, 1)|>collect)
        Q_coo[1:18, 3] .-= repeat(vec(pair[1] .* A[pair[2], 1:3, 1:3]), 2)
        Q_coo[19:24, 3] .+= repeat(vec(pair[1] .* A[pair[2], 1:3, 4]'* A[pair[2], 1:3, 1:3]), 2)
        Q_coo[25:30, 3] .-= repeat(vec(pair[1] .* A[pair[2], 1:3, 4]'), 2)
        Q_coo[31, 3] += pair[1] * A[pair[2], 1:3, 4]'*A[pair[2], 1:3, 4]
        Q_coo[32:85, 3] .-= repeat(pair[1] .* vec(kron(B[pair[2], 1:3, 4], A[pair[2], 1:3, 1:3])), 2)
        Q_coo[86:103, 3] .-= repeat(pair[1] .* kron(B[pair[2], 1:3, 4], A[pair[2], 1:3, 4]), 2)
    end

    #Concatenate entries
    return Q_coo
end

function get_row_cols_dense_terms_sparse_translation_cost()
    return [[repeat(range(4, 6)|>collect, 3) vec(repeat(range(1, 3)'|>collect, 3))];
    [vec(repeat(range(1, 3)'|>collect, 3)) repeat(range(4, 6)|>collect, 3)];
    [repeat([7], 3) range(1, 3)|>collect];
    [range(1, 3)|>collect repeat([7], 3)];
    [repeat([7], 3) range(4, 6)|>collect];
    [range(4, 6)|>collect repeat([7], 3)];
    [7 7];
    [repeat(range(17, 25)|>collect, 3) vec(repeat(range(1,3)'|>collect, 9))];
    [vec(repeat(range(1,3)'|>collect, 9)) repeat(range(17, 25)|>collect, 3)];
    [range(17,25)|>collect repeat([7], 9)];
    [repeat([7], 9) range(17,25)|>collect]]
end

function get_sparse_terms_sparse_translation_cost(B::AbstractArray, τ::Array)
    #Initialize temporary storage
    Q_coo =[get_row_cols_sparse_terms_sparse_translation_cost() zeros(45)]

    #Set sparse structure 1
    Q_coo[1:9, 3] = vec(repeat(sum(τ .* B[:, 1:3, 4], dims=1), 3))
    Q_coo[10:18, 3] = vec(repeat(sum(τ .* B[:, 1:3, 4], dims=1), 3))
    
    #Set sparse structure 2 diag
    Q_coo[19:27, 3] = vec(repeat(sum(τ .* B[:, 1:3, 4] .* B[:, 1:3, 4], dims=1), 3))

    #Set sparse structure 2 offdiag elements
    Q_coo[[28:30|>collect; 37:39|>collect], 3] = repeat(repeat(sum(τ .* B[:, 2, 4] .* B[:, 1, 4], dims=1), 3), 2)
    Q_coo[[31:33|>collect; 40:42|>collect], 3] = repeat(repeat(sum(τ .* B[:, 3, 4] .* B[:, 1, 4], dims=1), 3), 2)
    Q_coo[[34:36|>collect; 43:45|>collect], 3] = repeat(repeat(sum(τ .* B[:, 3, 4] .* B[:, 2, 4], dims=1), 3), 2)


    return Q_coo
end

function get_row_cols_sparse_terms_sparse_translation_cost()
    return [[range(17, 25)|>collect repeat([4; 5; 6], 3)];
    [repeat([4; 5; 6], 3) range(17, 25)|>collect];
    [range(17, 25)|>collect range(17, 25)|>collect];
    [[range(20, 25)|>collect; range(23, 25)] [range(17, 19)|>collect; range(17, 22)]];
    [[range(17, 19)|>collect; range(17, 22)] [range(20, 25)|>collect; range(23, 25)]]]
end

function map_scale_to_scaled(Qt_coo::Matrix)
    Qt_coo[Qt_coo[:, 1] .== 7, 1] .= 26
    Qt_coo[Qt_coo[:, 2] .== 7, 2] .= 26
    Qt_coo[Qt_coo[:, 1] .> 7, 1] .-= 1
    Qt_coo[Qt_coo[:, 2] .> 7, 2] .-= 1
    return Qt_coo
end

function sparse_robot_world_transformation_cost(A::AbstractArray, B::AbstractArray, κ::Array, τ::Array)
    #Get sparse rotation cost
    rot_cost_coo = [sparse_rotation_cost(A[:, 1:3, 1:3], B[:, 1:3, 1:3], false, κ); 25 25 0] #Add placeholder zero to ensure correct dimension
    Qr_sparse = sparse(rot_cost_coo[:, 1], rot_cost_coo[:, 2], rot_cost_coo[:, 3])

    #Get sparse translation cost
    trans_cost_coo = sparse_translation_cost(A, B, false, τ)
    Qt_sparse = sparse(trans_cost_coo[:, 1], trans_cost_coo[:, 2], trans_cost_coo[:, 3])
    Q = Qr_sparse + Qt_sparse
    return (Q + transpose(Q))/2
end

function sparse_robot_world_transformation_scale_cost(A::AbstractArray, B::AbstractArray, κ::Array, τ::Array)
    #Get sparse rotation cost
    rot_cost_coo = [sparse_rotation_cost(A[:, 1:3, 1:3], B[:, 1:3, 1:3], true, κ); 26 26 0] #Add homogeneous variable
    Qr_sparse = sparse(rot_cost_coo[:, 1], rot_cost_coo[:, 2], rot_cost_coo[:, 3])


    #Get sparse translation cost
    trans_cost_coo = [sparse_translation_cost(A, B, true, τ); 26 26 0] #Add homogeneous variable
    Qt_sparse = sparse(trans_cost_coo[:, 1], trans_cost_coo[:, 2], trans_cost_coo[:, 3])
    
    return Qr_sparse + Qt_sparse
end

function get_ith_eye_jth_base_cost(i::Int64, j::Int64, N::Int64, M::Int64, Q::AbstractArray, scale::Bool)
    #Unlike the constraints Q is changing for each pair, so each dense Q would need to be sparsified.
    #Would it be better to generate a sparse cost matrix from the start?
    nonzero_vals_temp = findnz(Q)
    nonzero_vals = [nonzero_vals_temp[1] nonzero_vals_temp[2] nonzero_vals_temp[3]];
    if any(nonzero_vals[:, 1] .== 26) == 0 & scale
        nonzero_vals = vcat(nonzero_vals, [26 26 0])
    end
    sparse_idxs = map_i_of_n_eyes_j_of_m_bases(nonzero_vals[:, 1:2], i, j, N, M, scale)
    return sparse(sparse_idxs[:, 1], sparse_idxs[:, 2], nonzero_vals[:, 3])
end

function map_i_of_n_eyes_j_of_m_bases(sparse_idxs::AbstractArray, i::Int64, j::Int64, N::Int64, M::Int64, scale::Bool)
    #sparse_idxs should be a Px2 matrix of non-zero indices
    discontinuities = get_parameter_ranges(scale)
    offsets = get_offsets(i, j, N, M, scale)
    for idx in eachindex(discontinuities)
        if length(discontinuities[idx]) == 1
            sparse_idxs[sparse_idxs .==discontinuities[idx]] .= offsets[idx]
        elseif length(discontinuities[idx]) ==2
            sparse_idxs[discontinuities[idx][1] .≤ sparse_idxs .&&  sparse_idxs .≤ discontinuities[idx][2]] .+= offsets[idx]
        end
    end
    return sparse_idxs
end

function get_parameter_ranges(scale::Bool)
    if scale
        return [26;
                (17, 25);
                (8,16);
                7;
                (4,6);
                (1,3)]
    else
        return [25;
                (16, 24);
                (7,15);
                (4,6);
                (1,3)]
    end
end

function get_offsets(i::Int64, j::Int64, N::Int64, M::Int64, scale::Bool)
    if scale
        return [12*(M+N)+2,
            3*(N+M)+9*N+9*(j-1)-15,
            3*(N+M)+9*(i-1)-6,
            3*(N+M)+1,
            3*N+3*(j-1)-3,
            3*(i-1)]
    else
        return [12*(M + N)+1,
                 12*N+3*M+9*(j-1)-15,
                 3*(N+M)+9*(i-1)-6,
                 3*N+3*(j-1)-3,
                 3*(i-1)]
    end
end

function get_empty_sparse_cost_matrix(N::Int64, M::Int64, scale::Bool)
    if scale
        return spzeros(12 * (M + N) + 2, 12 * (M + N) + 2)
    else
        return spzeros(12 * (M + N) + 1, 12 * (M + N) + 1)
    end
end