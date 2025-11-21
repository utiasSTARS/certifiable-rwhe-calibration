using Test
using Rotations
using LinearAlgebra
include("../src/calibration/hand_eye_costs.jl")

function norm_two(x)
    return norm(x, 2)
end

## Test random cases 
dim = 3
n_runs = 10
n_tests_per_run = 10
n_measurements_max = 10

# Rotations only
for i in 1:n_runs
    n = rand(2:n_measurements_max)
    # Form the measurements
    A = zeros(n, dim, dim)
    B = zeros(n, dim, dim)
    for i in 1:n
        Ai = rand(RotMatrix{3})
        A[i, :, :] = Ai
        Bi = rand(RotMatrix{3})
        B[i, :, :] = Bi
    end
    # Form the cost matrix
    M = rotation_hand_eye_cost(A, B)

    for j in 1:n_tests_per_run
        R_j = rand(RotMatrix{3})
        r_j = [reshape(R_j, 9); 1]
        matrix_cost = r_j'*M*r_j
        true_cost = sum([norm_two(A[i, :, :]*R_j - R_j*B[i, :, :])^2 for i in 1:n])

        @test abs(matrix_cost - true_cost) < 1e-9 
    end
end

# Full SE(3) poses 
for i in 1:n_runs
    n = rand(2:n_measurements_max)
    # Form the measurements
    A = zeros(n, dim + 1, dim + 1)
    A[:, 4, 4] .= 1
    B = zeros(n, dim + 1, dim + 1)
    B[:, 4, 4] .= 1
    for i in 1:n
        Ai = rand(RotMatrix{3})
        A[i, 1:3, 1:3] = Ai
        A[i, 1:3, 4] = rand(3)
        Bi = rand(RotMatrix{3})
        B[i, 1:3, 1:3] = Bi
        B[i, 1:3, 4] = rand(3)
    end
    # Form the cost matrix
    _, M = transformation_hand_eye_cost(A, B)

    for j in 1:n_tests_per_run
        R_j = rand(RotMatrix{3})
        t_j = rand(3)
        r_j = [reshape(R_j, 9); 1]
        x_j = [t_j; r_j]
        matrix_cost = x_j'*M*x_j
        true_rot_cost = sum([norm_two(A[i, 1:3, 1:3]*R_j - R_j*B[i, 1:3, 1:3])^2 for i in 1:n])
        true_trans_cost = sum([norm_two(R_j*B[i, 1:3, 4] + t_j - A[i, 1:3, 1:3]*t_j - A[i, 1:3, 4])^2 for i in 1:n])
        true_cost = true_rot_cost + true_trans_cost
        @test abs(matrix_cost - true_cost) < 1e-9 
    end
end

# Full SE(3) poses with unknown scale for sensor a
for i in 1:n_runs
    n = rand(2:n_measurements_max)
    # Form the measurements
    A = zeros(n, dim + 1, dim + 1)
    A[:, 4, 4] .= 1
    B = zeros(n, dim + 1, dim + 1)
    B[:, 4, 4] .= 1
    for i in 1:n
        Ai = rand(RotMatrix{3})
        A[i, 1:3, 1:3] = Ai
        A[i, 1:3, 4] = rand(3)
        Bi = rand(RotMatrix{3})
        B[i, 1:3, 1:3] = Bi
        B[i, 1:3, 4] = rand(3)
    end
    # Form the cost matrix
    _, M = transformation_hand_eye_cost_scale(A, B)

    for j in 1:n_tests_per_run
        R_j = rand(RotMatrix{3})
        t_j = rand(3)
        α = rand()
        r_j = [reshape(R_j, 9); 1]
        x_j = [t_j; α; r_j]
        matrix_cost = x_j'*M*x_j
        true_rot_cost = sum([norm_two(A[i, 1:3, 1:3]*R_j - R_j*B[i, 1:3, 1:3])^2 for i in 1:n])
        true_trans_cost = sum([norm_two(R_j*B[i, 1:3, 4] + t_j - A[i, 1:3, 1:3]*t_j - α*A[i, 1:3, 4])^2 for i in 1:n])
        true_cost = true_rot_cost + true_trans_cost
        @test abs(matrix_cost - true_cost) < 1e-9
    end
end
