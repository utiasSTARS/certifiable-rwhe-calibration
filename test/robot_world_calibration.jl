using Test
using Rotations

include("../src/calibration/robot_world_costs.jl")
include("../src/rotation_sdp_solver.jl")

function test_exact_rotation_measurements(X, Y, n)
    dim = size(X, 1)
    A = zeros(n, dim, dim)
    B = zeros(n, dim, dim)
    for i in 1:n
        Ai = rand(RotMatrix{3})
        A[i, :, :] = Ai
        B[i, :, :] = Y'*Ai*X
    end
    cost = rotation_robot_world_cost(A, B)
    Z, _ = solve_double_rotation_sdp(cost)
    extract_double_rotation(Z)
end

function test_exact_transformation_measurements(X, Y, n)
    dim = size(X, 1) - 1
    A = zeros(n, dim+1, dim+1)
    B = zeros(n, dim+1, dim+1)
    for i in 1:n
        Ai = [rand(RotMatrix{3}) (rand(3)*0.2 .- 1.0);
              0. 0. 0. 1.]
        A[i, :, :] = Ai
        B[i, :, :] = inv(Y)*Ai*X
    end
    cost, Q_full = transformation_robot_world_cost(A, B)
    Z, prob = solve_double_rotation_sdp(cost)
    Rx, Ry = extract_double_rotation(Z)
    t = translations_from_rotations(Rx, Ry, Q_full)
    tx = t[1:3]
    ty = t[4:6]
    return [Rx tx; 0. 0. 0. 1.], [Ry ty; 0. 0. 0. 1.]
end

## Test case for 3D (noiseless)
n = 3
X = RotX(π/2)
Y = RotY(π/2)
A1 = RotY(π/2)
A2 = RotZ(π/2)
A3 = RotX(π/2)
B1 = Y'*A1*X 
B2 = Y'*A2*X
B3 = Y'*A3*X

A = zeros(n, 3, 3)
A[1, :, :] = A1
A[2, :, :] = A2
A[3, :, :] = A3
B = zeros(n, 3, 3)
B[1, :, :] = B1
B[2, :, :] = B2
B[3, :, :] = B3

cost = rotation_robot_world_cost(A, B)
Z, prob = solve_double_rotation_sdp(cost)
Rx_sol, Ry_sol = extract_double_rotation(Z)
@test all(abs.(X - Rx_sol) .< 1e-9)
@test all(abs.(Y - Ry_sol) .< 1e-9)


## Test random cases 
n_runs = 10
n_measurements_max = 10

for i in 1:n_runs
    n = rand(3:n_measurements_max)
    X_i = rand(RotMatrix{3})
    Y_i = rand(RotMatrix{3})
    X_sol_i, Y_sol_i = test_exact_rotation_measurements(X_i, Y_i, n)
    @test all(abs.(X_i - X_sol_i) .< 1e-5)
    @test all(abs.(Y_i - Y_sol_i) .< 1e-5)
end

## Test random cases for full SE(3) problems
for i in 1:n_runs
    n = rand(3:n_measurements_max)
    R_x_i = rand(RotMatrix{3})
    X_i = [R_x_i (rand(3)*0.2 .- 1.0);
           0. 0. 0. 1.]
    R_y_i = rand(RotMatrix{3})
    Y_i = [R_y_i (rand(3)*0.2 .- 1.0);
           0. 0. 0. 1.]
    X_sol_i, Y_sol_i = test_exact_transformation_measurements(X_i, Y_i, n)
    println("Iteration "*string(i))
    println(X_i)
    println(X_sol_i)
    @test all(abs.(X_i - X_sol_i) .< 1e-5)
    @test all(abs.(Y_i - Y_sol_i) .< 1e-5)
end

# Sometimes the solver fails (need to investigate observability)