using DataFrames
# using COSMO
# using COSMOAccelerators
using MosekTools
using LinearAlgebra
using Statistics
using Rotations

include("../src/local_solver/local_lie_group_hand_eye_solver.jl")
include("../src/rotation_sdp_solver.jl")
include("../src/calibration/hand_eye_costs.jl")
include("../src/utils/experiment_setup.jl")
include("../src/local_solver/linear_initializations.jl")


# Experiment parameters
n_runs = 10
n_motions = 20
# These parameters describe the distribution of the platform's groundtruth motion
mean_rot = 20*π/180
rot_std_dev = 5*π/180
mean_trans = 0.1
trans_std_dev = 0.025

# These parameters describe the distribution of zero-mean measurement errors
var_r_vec = ([0.05, 0.2, 0.5]*mean_rot).^2  # Lie-algebraic (gets converted to Langevins)
var_t_vec = ([0.05, 0.2, 0.5]*mean_trans).^2

use_column_orthog = true
use_handedness = true

# Store data for analysis
n_data = n_runs*length(var_r_vec)*length(var_t_vec)

trans_std_dev_vec = zeros(n_data)
rot_std_dev_vec = zeros(n_data) 
time_sdp = zeros(n_data)
eig_ratio_sdp = zeros(n_data)
rot_err_sdp = zeros(n_data)
trans_err_sdp = zeros(n_data)

time_local = zeros(n_data)
rot_err_local = zeros(n_data)
trans_err_local = zeros(n_data)
local_iters_vec = zeros(n_data)
rot_err_init = zeros(n_data)
trans_err_init = zeros(n_data)


ind = 0
for vr in var_r_vec
    κ = convert_lie_variance_to_langevins_concentration(vr)
    for vt in var_t_vec
        for i in 1:n_runs
            global ind += 1
            T_gt = [rand(RotMatrix{3}) 2*rand(3).-1; 0 0 0 1]
            # Last boolean flag is Langevins (uses Lie algebra for noise otherwise)
            A, B = random_test_instance(κ, vt, mean_rot, rot_std_dev, mean_trans, trans_std_dev, T_gt, n_motions, true)
            A_mat = vec_of_matrices_to_tensor(A)
            B_mat = vec_of_matrices_to_tensor(B)
            # Form cost for SDP
            cost_matrix, Q = transformation_hand_eye_cost(A_mat, B_mat, κ*ones(n_motions), (1/vt)*ones(n_motions))

            # COSMO.Optimizer or Mosek.Optimizer
            global time_sdp[ind] = @elapsed Z_i_sdp, prob_sdp = solve_rotation_sdp(cost_matrix, use_column_orthog, use_handedness, Mosek.Optimizer)  
            Z_eigs_sdp = eigvals(Z_i_sdp)
            global eig_ratio_sdp[ind] = Z_eigs_sdp[end]/Z_eigs_sdp[end-1]
            R_sdp = extract_rotation(Z_i_sdp)
            t_sdp = translation_from_rotation(R_sdp, Q)
            global rot_err_sdp[ind] = abs(rotation_angle(RotMatrix{3}(R_sdp'*T_gt[1:3, 1:3])))
            global trans_err_sdp[ind] = norm(t_sdp - T_gt[1:3, 4])

            # Linearized init. + local solver
            R_init, t_init = linear_hand_eye_initialization(A, B, vt*ones(n_motions), vr*ones(n_motions))
            T_init = [R_init t_init; 0 0 0 1]
            global rot_err_init[ind] = abs(rotation_angle(RotMatrix{3}(T_init[1:3, 1:3]'*T_gt[1:3, 1:3])))
            global trans_err_init[ind] = norm(T_init[1:3, 4] - T_gt[1:3, 4])
            T_init[1:3, 1:3] = T_init[1:3, 1:3]*random_rotation_sample_normal_magnitude(0., 0.1)
            T_init[1:3, 4] = T_init[1:3, 4] + randn(3)*0.2
            global time_local[ind] = @elapsed T_lin, local_iters = solve_SE3_linearized_calibration(A, B, Pose3D(T_init), vt, vr, 20, 1e-6)
            # Local solver is TERRIBLE with identity as the initializer
            # global time_local[ind] = @elapsed T_lin, local_iters = solve_SE3_linearized_calibration(A, B, Pose3D(I(4)), vt, vr, 100, 1e-6)
            global rot_err_local[ind] = abs(rotation_angle(RotMatrix{3}(T_lin[1:3, 1:3]'*T_gt[1:3, 1:3])))
            global trans_err_local[ind] = norm(T_lin[1:3, 4] - T_gt[1:3, 4])
            global rot_std_dev_vec[ind] = sqrt(vr)
            global trans_std_dev_vec[ind] = sqrt(vt)
            global local_iters_vec[ind] = local_iters
        end
    end
end

# Create dataframe
df = DataFrame(rot_std_dev=rot_std_dev_vec, trans_std_dev=trans_std_dev_vec, runtime_sdp=time_sdp, runtime_local=time_local, eig_ratio=eig_ratio_sdp, 
               rot_err_sdp=rot_err_sdp, rot_err_local=rot_err_local, rot_err_init=rot_err_init, local_iters=local_iters_vec,
               trans_err_sdp=trans_err_sdp, trans_err_local=trans_err_local, trans_err_init=trans_err_init)

# Example of how to get summaries:
describe(df[(df.rot_std_dev .≈ 0.017453292519943295) .& (df.trans_std_dev .≈ 0.005), :].rot_err_init)
describe(df[(df.rot_std_dev .≈ 0.017453292519943295) .& (df.trans_std_dev .≈ 0.005), :].rot_err_local)
describe(df[(df.rot_std_dev .≈ 0.017453292519943295) .& (df.trans_std_dev .≈ 0.005), :].rot_err_sdp)

describe(df[(df.rot_std_dev .≈ 0.0174532925) .& (df.trans_std_dev .≈ 0.005), :].trans_err_init)
describe(df[(df.rot_std_dev .≈ 0.0174532925) .& (df.trans_std_dev .≈ 0.005), :].trans_err_local)
describe(df[(df.rot_std_dev .≈ 0.0174532925) .& (df.trans_std_dev .≈ 0.005), :].trans_err_sdp)
