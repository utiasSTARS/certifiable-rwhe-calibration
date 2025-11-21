using Test
using Rotations

include("../src/calibration/hand_eye_costs.jl")
include("../src/rotation_sdp_solver.jl")

function test_exact_rotation_measurements(X, n)
    dim = size(X, 1)
    A = zeros(n, dim, dim)
    B = zeros(n, dim, dim)
    for i in 1:n
        Ai = rand(RotMatrix{3})
        A[i, :, :] = Ai
        Bi = X'*Ai*X
        B[i, :, :] = Bi
    end
    cost = rotation_hand_eye_cost(A, B)
    Z, prob = solve_rotation_sdp(cost)
    extract_rotation(Z)
end

function test_exact_transformation_measurements(X, n)
    dim = size(X, 1) - 1
    A = zeros(n, dim+1, dim+1)
    B = zeros(n, dim+1, dim+1)
    for i in 1:n
        Ai = [rand(RotMatrix{3}) (rand(3)*0.2 .- 1.0);
              0. 0. 0. 1.]
        A[i, :, :] = Ai
        Bi = inv(X)*Ai*X
        B[i, :, :] = Bi
    end
    cost, Q_full = transformation_hand_eye_cost(A, B)
    Z, prob = solve_rotation_sdp(cost)
    R = extract_rotation(Z)
    t = translation_from_rotation(R, Q_full)
    return [R t; 0. 0. 0. 1.]
end

## Test case for 3D (noiseless)
X = RotX(π/2)
A1 = RotY(π/2)
A2 = RotZ(π/2)
B1 = X'*A1*X 
B2 = X'*A2*X

A = zeros(2, 3, 3)
A[1, :, :] = A1
A[2, :, :] = A2
B = zeros(2, 3, 3)
B[1, :, :] = B1
B[2, :, :] = B2

cost = rotation_hand_eye_cost(A, B)
Z, prob = solve_rotation_sdp(cost)
X_sol = extract_rotation(Z)
@test all(abs.(X - X_sol) .< 1e-9)

## Test random cases 
n_runs = 10
n_measurements_max = 10

for i in 1:n_runs
    n = rand(2:n_measurements_max)
    X_i = rand(RotMatrix{3})
    X_sol_i = test_exact_rotation_measurements(X_i, n)
    @test all(abs.(X_i - X_sol_i) .< 1e-5)
end

## Test random cases for full SE(3) problems
for i in 1:n_runs
    n = rand(2:n_measurements_max)
    R_i = rand(RotMatrix{3})
    X_i = [R_i (rand(3)*0.2 .- 1.0);
           0. 0. 0. 1.]
    X_sol_i = test_exact_transformation_measurements(X_i, n)
    @test all(abs.(X_i - X_sol_i) .< 1e-5)
end
