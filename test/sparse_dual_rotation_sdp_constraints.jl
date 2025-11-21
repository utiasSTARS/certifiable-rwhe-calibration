using Test
using Rotations
using SparseArrays

include("../src/rotation_sdp_constraints.jl")

function test_sparse_single_constraints(R1::RotMatrix3)
    P_dict = get_sparse_dual_rotation_constraints()
    r = [vec(R1); 1]
    key_list = keys(P_dict)|>collect
    display(length(key_list))
    [r'*sparse(P_dict[p_key][:, 1], P_dict[p_key][:, 2], P_dict[p_key][:, 3])*r for p_key in keys(P_dict)]
end

function test_sparse_ith_of_n_constraints(i::Int64)
    n = i + rand(1:10)
    P_dict = get_i_of_n_sparse_rotation_constraints(i, n)
    r = rand(3*n + 9*n + 1)
    r[3*n + 9 * (i - 1) + 1: 3*n + 9 * i] = vec(rand(RotMatrix{3}))
    r[end] = 1
    [r'*sparse(P_dict[p_key][:, 1], P_dict[p_key][:, 2], P_dict[p_key][:, 3])*r for p_key in keys(P_dict)]
end

function test_sparse_ith_of_n_constraints_w_scale(i::Int64)
    n = i + rand(1:10)
    P_dict = get_i_of_n_sparse_rotation_constraints_w_scale(i, n)
    r = rand(3*n + 1 + 9*n + 1)
    r[3*n + 1 + 9 * (i - 1) + 1: 3*n + 1 + 9 * i] = vec(rand(RotMatrix{3}))
    r[end] = 1
    [r'*sparse(P_dict[p_key][:, 1], P_dict[p_key][:, 2], P_dict[p_key][:, 3])*r for p_key in keys(P_dict)]
end

function test_sparse_n_constraints(n::Int64)
    P_dict = get_n_sparse_rotation_constraints(n, false)
    r = rand(3 * n + 9 * n + 1)
    for i =1:n
        r[3*n + 9*(i - 1) + 1: 3*n + 9*i] = vec(rand(RotMatrix{3}))
    end
    r[end] = 1
    [r'*sparse(P_dict[p_key][:, 1], P_dict[p_key][:, 2], P_dict[p_key][:, 3])*r for p_key in keys(P_dict)]
end

function test_sparse_dual_dense_primal_equivalency()
    #Specifically for the Robot-World Problem
    #Collect Primal Constraints
    A, b = get_double_rotation_constraints()
    A_full = zeros(size(A, 1), 25, 25)
    A_full[:, 7:25, 7:25] = A
    A_full[1:end-1, 25, 25] = -b[1:size(A, 1)-1]

    #Collect Dual Constraints
    P_dict = get_n_sparse_rotation_constraints(2, false)
    P_dict = add_sparse_homogeneous_constraint(25, P_dict)
    result = []
    for key in sort(keys(P_dict)|>collect)
        tuple_of_nz_values = findnz(sparse(A_full[key, :, :]))
        temp_A = hcat(tuple_of_nz_values[1], tuple_of_nz_values[2], tuple_of_nz_values[3])
        if size(P_dict[key], 1) > size(temp_A, 1)
            result = [result; ~iszero(P_dict[key][1:end-1, :] - temp_A)]
        else
            result = [result; ~iszero(P_dict[key] - temp_A)]
        end
    end

    result
end

@test all([≈(val, 0, atol=1e-9) for val in test_sparse_single_constraints(rand(RotMatrix{3}))])
@test all([≈(val, 0, atol=1e-9) for val in test_sparse_ith_of_n_constraints(rand(1:10))])
@test all([≈(val, 0, atol=1e-9) for val in test_sparse_ith_of_n_constraints_w_scale(rand(1:10))])
@test all([≈(val, 0, atol=1e-9) for val in test_sparse_n_constraints(rand(2:10))])
@test all([≈(val, 0, atol=1e-9) for val in test_sparse_dual_dense_primal_equivalency()])