using Test
using Rotations
using SparseArrays

include("../src/calibration/robot_world_costs.jl")

function test_rotation_cost()
    #Generate a 3 pose problem
    X = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    Y = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    B = permutedims([[rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]], [3, 1, 2])
    A = zeros(3, 4, 4)
    for i=1:3
        A[i, :, :] = Y * B[i, :, :] * inv(X)
    end

    #Get sparse rotation cost
    rot_cost_coo = [sparse_rotation_cost(A[:, 1:3, 1:3], B[:, 1:3, 1:3], false, ones(3)); 25 25 0]
    Q_sparse = sparse(rot_cost_coo[:, 1], rot_cost_coo[:, 2], rot_cost_coo[:, 3])

    #Get dense rotation cost
    Q_dense = [zeros(6, 6) zeros(6, 19);
    zeros(19, 6) rotation_robot_world_cost(A[:, 1:3, 1:3], B[:, 1:3, 1:3], ones(3))]

    return vec(Q_dense - Q_sparse)
end

function test_translation_cost()
    #Generate a 3 pose problem
    X = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    Y = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    B = permutedims([[rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]], [3, 1, 2])
    A = zeros(3, 4, 4)
    for i=1:3
       A[i, :, :] = Y * B[i, :, :] * inv(X)
    end
    τ = ones(3)

    #Get sparse cost
    trans_cost_coo = sparse_translation_cost(A, B, false, τ)
    Q_sparse = sparse(trans_cost_coo[:, 1], trans_cost_coo[:, 2], trans_cost_coo[:, 3])

    #Get dense cost
    Q_dense = translation_cost(A, B, τ)

    return vec(Q_dense - Q_sparse)
end

@test all([≈(val, 0, atol=1e-9) for val in test_rotation_cost()])
@test all([≈(val, 0, atol=1e-9) for val in test_translation_cost()])
