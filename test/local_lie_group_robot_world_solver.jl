using Test
using Rotations

include("../src/local_solver/local_lie_group_robot_world_solver.jl")

# Test single sensor AX=YB with no noise
n = 3
X = Pose3D([RotX(π/2) rand(3); 0 0 0 1])
Y = Pose3D([RotY(π/2) rand(3); 0 0 0 1])
A1 = Pose3D([RotY(π/2) rand(3); 0 0 0 1])
A2 = Pose3D([RotZ(π/2) rand(3); 0 0 0 1])
A3 = Pose3D([RotX(π/2) rand(3); 0 0 0 1])
B1 = inv_SE3(Y)*A1*X 
B2 = inv_SE3(Y)*A2*X
B3 = inv_SE3(Y)*A3*X

A = zeros(n, 4, 4)
A[1, :, :] = A1
A[2, :, :] = A2
A[3, :, :] = A3
B = zeros(n, 4, 4)
B[1, :, :] = B1
B[2, :, :] = B2
B[3, :, :] = B3
A = [A]
B = [B]

# Give a good initialization
dx = se3(0.1*ones(6))
dy = copy(dx)
dy[1:3] .= 0.0
X0 = [apply_perturbation(X, dx)]
Y0 = apply_perturbation(Y, dy)

X_sol, Y_sol = solve_SE3_linearized_calibration(A, B, X0, Y0, 2)