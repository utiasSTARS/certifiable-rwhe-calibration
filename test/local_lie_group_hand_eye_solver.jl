using Test
using Rotations

include("../src/local_solver/local_lie_group_hand_eye_solver.jl")

# Test AX=XB with no noise
X = Pose3D([RotX(π/2) rand(3); 0 0 0 1])
A1 = Pose3D([RotY(π/2) rand(3); 0 0 0 1])
A2 = Pose3D([RotZ(π/2) rand(3); 0 0 0 1])
A3 = Pose3D([RotX(π/2) rand(3); 0 0 0 1])
B1 = inv_SE3(X)*A1*X 
B2 = inv_SE3(X)*A2*X
B3 = inv_SE3(X)*A3*X

A = [A1, A2, A3]
B = [B1, B2, B3]

# Give a good initialization
dx = se3(1.5*ones(6))  # At 1.5 and above, seems to drive it to a local minimum (even in this noise free case)
                       # Check (via Jacobian) that it is in fact a local minimum?
# dx[4:6] .= 2.0
X0 = apply_perturbation(X, dx)
X0[1:3, 4] = X[1:3, 4] # Assume that only the initial rotation guess is poor

X_sol = solve_SE3_linearized_calibration(A, B, X0, 10)
