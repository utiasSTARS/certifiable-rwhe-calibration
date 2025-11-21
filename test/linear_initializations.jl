using Test
using Rotations

include("../src/local_solver/linear_initializations.jl")
include("../src/local_solver/lie_theory_helpers.jl")

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

R, t = linear_hand_eye_initialization(A, B)
X_sol = [R t; 0 0 0 1]