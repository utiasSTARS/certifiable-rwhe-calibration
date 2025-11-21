using MKL
using Test
using Rotations
using SparseArrays

include("../src/calibration/robot_world_costs.jl");
include("../src/rotation_sdp_solver_jump.jl");

function test_solver()
    #Generate a 3 pose problem
    α = rand()
    X = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    Y = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    B = permutedims([[rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]], [3, 1, 2])
    A = zeros(3, 4, 4)
    for i=1:3
       A[i, :, :] = Y * B[i, :, :] * inv(X)
    end
    τ = ones(3)
    κ = ones(3)

    #Get sparse cost
    Q = sparse_robot_world_transformation_cost(A, B, κ, τ)
    Z, model = solve_sdp_dual_jump(Q, 2, true, true)
    solution = extract_solution_from_dual(Z);
    
   return [solution[1:3] - X[1:3, 4]; 
   solution[4:6] - Y[1:3, 4]; 
   vec(reshape(solution[7:15], 3, 3) * X[1:3, 1:3]' - Matrix(I, 3, 3)); 
   vec(reshape(solution[16:24], 3, 3) * Y[1:3, 1:3]' - Matrix(I, 3, 3))]
end

@test all([≈(val, 0, atol=1e-8) for val in test_solver()])
