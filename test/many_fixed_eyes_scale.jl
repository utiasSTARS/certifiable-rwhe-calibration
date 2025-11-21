using MKL
using Test
using Profile
using Rotations
using SparseArrays

include("../src/calibration/robot_world_costs.jl");
#include("../src/rotation_sdp_solver.jl");
include("../src/rotation_sdp_solver_jump.jl");

function generate_1000_fixed_eyes()
    n = 20
    #Generate Y and A
    X = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    A = permutedims([[rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; 
                     [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; 
                     [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]], [3, 1, 2])
    
    #Initialize X and B
    Y = zeros(n, 4, 4)
    B = zeros(n, 3, 4, 4)

    #Generate X and B
    Y[:, 4, 4] .= 1
    Y[:, 1:3, 4] = rand(Float64, n, 3)
    for i=1:n
        Y[i, 1:3, 1:3] = rand(RotMatrix{3})
        for j=1:3
            B[i, j, :, :] = inv(Y[i, :, :]) * A[j, :, :] * X
        end
    end
    τ = ones(n, 3)
    κ = ones(n, 3)
    return A, B, X, Y, τ, κ, n
end


function test_multi_eye_cost()
    A, B, X, Y, τ, κ, n = generate_1000_fixed_eyes()
    display("Generated Data")
    Q = get_empty_sparse_cost_matrix(n, 1, true)
    for i=1:n
        Q_temp = sparse_robot_world_transformation_scale_cost(A, B[i, :, :, :], κ[i, :], τ[i, :])
        Q += get_ith_eye_jth_base_cost(1, i, 1, n, Q_temp, true)
    end
    display("Generated Cost Matrix")

    display("Percentage Non-Zero Entries in Full Matrix")
    display(count(!iszero, Q)/(12*(n+1) + 2)^2)
    display("Percentage Non-Zero Entries in Rotation Submatrix")
    display(count(!iszero, Q[3*(n+1)+2:end, 3*(n+1)+2:end])/(9*(n+ 1) + 1)^2)
    display("Percentage Non-Zero Entries in Schur Complement")
    display(count(!iszero, Q[3*(n+1)+2:end, 3*(n+1)+2:end] - Q[3*(n+1)+2:end, 1:3*(n+1)+1] * sparse(inv(Matrix(Q[1:3*(n+1)+1, 1:3*(n+1)+1]))) * Q[1:3*(n+1)+1, 3*(n+1)+2:end])/(9*(n+ 1) + 1)^2)

    #Older Solvers (Much Slower)
    #MOI = Convex.MOI
    #solver = MOI.OptimizerWithAttributes(COSMO.Optimizer, "eps_abs" => 1e-10, "eps_rel" =>1e-10)
    
    #Generates 21*n sparse matrices for constraints using Convex
    #Z, prob = solve_sdp_dual(Q, n + 1, true, true, solver)
    
    #Generates Constraints through concatenation using Convex
    #Z, prob = solve_sdp_dual_using_concat(Q, n + 1, true, true, solver)

    #Generates 21*n sparse matrices for constraints using JuMP
    #Z, model = solve_sdp_dual_jump(Q, n + 1, true, true)

    #Fastest Solver
    #Generates Constraints through concatenation using Convex
    Z, model = solve_sdp_dual_w_scale_using_concat(Q, n + 1, true, true)
    #Z_schur, model_schur = solve_sdp_dual_w_scale_using_schur(Q, n + 1, true, true)

    #solution = extract_solution_from_dual_schur(Z_schur, Q)
    solution = extract_solution_from_dual(Z)
    gt = vcat(collect(Iterators.flatten(X[1:3, 4])),
            collect(Iterators.flatten(permutedims(Y[:, 1:3, 4], [2, 1]))),
            1.0,
            collect(Iterators.flatten(X[1:3, 1:3])),
            collect(Iterators.flatten(permutedims(Y[:, 1:3, 1:3], [2, 3, 1]))),
            1.0)
    return gt - solution
end

@test all([≈(val, 0, atol=1e-6) for val in test_multi_eye_cost()])
#Profile.clear()
#f = open("multi_eye_profile.txt", "w")
#@profile test_multi_eye_cost()
#Profile.print(f)
#close(f)