using Test
using Rotations
using LinearAlgebra
include("../src/calibration/robot_world_costs.jl")

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
    n = rand(3:n_measurements_max)
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
    M = rotation_robot_world_cost(A, B)

    for j in 1:n_tests_per_run
        R_x = rand(RotMatrix{3})
        r_x = reshape(R_x, 9)
        R_y = rand(RotMatrix{3})
        r_y = [reshape(R_y, 9); 1]
        r = [r_x; r_y]
        matrix_cost = r'*M*r
        true_cost = sum([norm_two(A[i, :, :]*R_x - R_y*B[i, :, :])^2 for i in 1:n])

        @test abs(matrix_cost - true_cost) < 1e-9 
    end
end

# Full SE(3) poses 
for i in 1:n_runs
    n = rand(3:n_measurements_max)
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
    _, M = transformation_robot_world_cost(A, B)

    for j in 1:n_tests_per_run
        R_x = rand(RotMatrix{3})
        r_x = reshape(R_x, 9)
        R_y = rand(RotMatrix{3})
        r_y = [reshape(R_y, 9); 1]
        r = [r_x; r_y]
        t_x = rand(3)
        t_y = rand(3)
        x = [t_x; t_y; r_x; r_y]
        matrix_cost = x'*M*x
        true_rot_cost = sum([norm_two(A[i, 1:3, 1:3]*R_x - R_y*B[i, 1:3, 1:3])^2 for i in 1:n])
        true_trans_cost = sum([norm_two(R_y*B[i, 1:3, 4] + t_y - A[i, 1:3, 1:3]*t_x - A[i, 1:3, 4])^2 for i in 1:n])
        true_cost = true_rot_cost + true_trans_cost
        @test abs(matrix_cost - true_cost) < 1e-9 
    end
end

# Full SE(3) poses with unknown scale for sensor a
for i in 1:n_runs
    n = rand(3:n_measurements_max)
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
    _, M = transformation_robot_world_cost_scale(A, B)

    for j in 1:n_tests_per_run
        R_x = rand(RotMatrix{3})
        r_x = reshape(R_x, 9)
        R_y = rand(RotMatrix{3})
        r_y = [reshape(R_y, 9); 1]
        r = [r_x; r_y]
        t_x = rand(3)
        t_y = rand(3)
        α = rand()
        x = [t_x; t_y; α; r_x; r_y]

        matrix_cost = x'*M*x
        true_rot_cost = sum([norm_two(A[i, 1:3, 1:3]*R_x - R_y*B[i, 1:3, 1:3])^2 for i in 1:n])
        true_trans_cost = sum([norm_two(R_y*B[i, 1:3, 4] + t_y - A[i, 1:3, 1:3]*t_x - α*A[i, 1:3, 4])^2 for i in 1:n])
        true_cost = true_rot_cost + true_trans_cost
        @test abs(matrix_cost - true_cost) < 1e-9
    end
end
