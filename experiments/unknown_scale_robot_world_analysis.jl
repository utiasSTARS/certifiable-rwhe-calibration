using Rotations
include("../src/calibration/robot_world_costs.jl")
include("../src/rotation_sdp_solver.jl")

# Setup a minimal (3 measurement) noiseless base case
n = 3
R_x = RotZ(π/2)  # Ground truth extrinsic rotation
t_x = [0; 1; 0]  # Ground truth extrinsic translation
T_x = [R_x t_x; 0 0 0 1]
R_y = RotX(π/2)  # Ground truth extrinsic rotation
t_y = [1; 0; 0]  # Ground truth extrinsic translation
T_y = [R_y t_y; 0 0 0 1]
α = 0.5

Ta = zeros(n, 4, 4)
Ta[:, 4, 4] .= 1
Ta[1, 1:3, 1:3] = RotX(π/2)
Ta[1, 1:3, 4]   = [0; 0; 1]
Ta[2, 1:3, 1:3] = RotY(π/2)
Ta[2, 1:3, 4] = [0; 0; 1]
Ta[3, 1:3, 1:3] = RotZ(π/2)
Ta[3, 1:3, 4] = [1; 0; 0]

Tb = zeros(size(Ta))
Tb[1, :, :] = inv(T_y)*Ta[1, :, :]*T_x
Tb[2, :, :] = inv(T_y)*Ta[2, :, :]*T_x
Tb[3, :, :] = inv(T_y)*Ta[3, :, :]*T_x

# Apply the scale 
Ta_known_α = copy(Ta)
Ta[:, 1:3, 4] = Ta[:, 1:3, 4]/α

# Check the inverted element for the Schur complement:
M_tα1 = [Ta[1, 1:3, 1:3] -I(3) Ta[1, 1:3, 4]]
Q_tα  = M_tα1'*M_tα1
M_tα2 = [Ta[2, 1:3, 1:3] -I(3) Ta[2, 1:3, 4]]
M_tα3 = [Ta[3, 1:3, 1:3] -I(3) Ta[3, 1:3, 4]]
Q_tα  += M_tα2'*M_tα2 + M_tα3'*M_tα3


# # See if using the nullspace of Q_tα leads to the same cost function:
# δ = 1.5
# t_perturbed = t + δ*eigvecs(Q_tα)[1:3, 1]
# α_perturbed = α + δ*eigvecs(Q_tα)[4, 1]

# true_trans_cost = sum([norm2(R*Tb[i, 1:3, 4] + t - Ta[i, 1:3, 1:3]*t - α*Ta[i, 1:3, 4])^2 for i in 1:2])

# perturbed_trans_cost = sum([norm2(R*Tb[i, 1:3, 4] + t_perturbed - Ta[i, 1:3, 1:3]*t_perturbed - α_perturbed*Ta[i, 1:3, 4])^2 for i in 1:2])


# Gives inversion error for now
# M_rot, M_full = transformation_robot_world_cost_scale(Ta, Tb)

# sol, prob = solve_rotation_sdp(M_rot)

