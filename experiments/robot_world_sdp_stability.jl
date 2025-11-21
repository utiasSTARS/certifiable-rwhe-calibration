""" Setup an experiment to show the zero-duality-gap (ZDG) region of robot-world calibration.
"""

include("../src/calibration/robot_world_costs.jl")
include("../src/rotation_sdp_solver.jl")
using Plots
using GR
using DataFrames
import GR.heatmap

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
Ra = Ta[:, 1:3, 1:3]

Ta_1_saved = copy(Ta[1, :, :])

Tb = zeros(size(Ta))
Tb[1, :, :] = inv(T_y)*Ta[1, :, :]*T_x
Tb[2, :, :] = inv(T_y)*Ta[2, :, :]*T_x
Tb[3, :, :] = inv(T_y)*Ta[3, :, :]*T_x
Rb = Tb[:, 1:3, 1:3]

# Perturb 2 translation elements in the first rotation to start
n_grid = 21 # 101
ϵ_max = 10.0 # π
ϵ_x = range(-ϵ_max, stop=ϵ_max, length=n_grid)
ϵ_y = ϵ_x

# Solver constraint settings
col_orthog = false
handed = true

# Base cost function
cost_gt = rotation_robot_world_cost(Ra, Rb)

eig_ratio_results = zeros(n_grid, n_grid)
determinant = zeros(size(eig_ratio_results))
orthog_resid = zeros(size(eig_ratio_results))
solver_status = zeros(size(eig_ratio_results))

# Solve for each perturbation and record the duality gap (and other stats)
for (i, dx) in enumerate(ϵ_x)
    for (j, dy) in enumerate(ϵ_y)
        Ta[1, 1, 4] = dx
        Ta[1, 2, 4] = dy
        cost, Q_full = transformation_robot_world_cost(Ta, Tb)

        # Ra[1, 1:3, 1:3] = Ta_1_saved[1:3, 1:3]*RotX(dx)*RotY(dy)
        # cost = rotation_robot_world_cost(Ra, Rb)

        # cost = copy(cost_gt)
        # cost[1, 1] += dx 
        # cost[2, 2] += dy

        Z, prob = solve_double_rotation_sdp(cost, col_orthog, handed, Mosek.Optimizer)
        R1, R2 = extract_double_rotation(Z)
        # global determinant[i, j] = 0.5*(real(det(R1)) + real(det(R2)))
        Z_eigs_sdp = eigvals(Z)
        global eig_ratio_results[i, j] = abs(real(Z_eigs_sdp[end]/Z_eigs_sdp[end-1]))
        global orthog_resid[i, j] = max(norm(I(3) - R1*R1'), norm(I(3) - R2*R2'))
        global solver_status[i, j] = Int(prob.status)

        # t = translations_from_rotations(Rx, Ry, Q_full)
        # tx = t[1:3]
        # ty = t[4:6]
        
    end
end

# Solve the "base case" (no noise)
# Ra[1, 1:3, 1:3] = Ta_1_saved[1:3, 1:3]
# cost_noise_free = rotation_robot_world_cost(Ra, Rb)
# Z_noise_free, prob_noise_free = solve_double_rotation_sdp(cost_noise_free, col_orthog, handed)
# R1_noise_free, R2_noise_free = extract_double_rotation(Z_noise_free)

# Visualize as a heatmap centered around the noise-free case
heatmap(orthog_resid)