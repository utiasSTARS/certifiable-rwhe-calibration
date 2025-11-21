using DataFrames
using COSMO
using COSMOAccelerators
using MosekTools
using LinearAlgebra
using Statistics
include("../src/rotation_sdp_solver.jl")
include("../src/low_rank_solver.jl")

# Experiment parameters
# TODO: add duality gap
n_runs = 10
use_column_orthog = true
use_handedness = true

# BM parameters
bm_rank = 9

eig_ratio_sdp = zeros(n_runs)
runtime_sdp = zeros(n_runs)
orthog_resid_sdp = zeros(n_runs)
determinant_sdp = zeros(n_runs)
cost_sdp = zeros(n_runs)

eig_ratio_bm = zeros(n_runs)
runtime_bm = zeros(n_runs)
orthog_resid_bm = zeros(n_runs)
determinant_bm = zeros(n_runs)
cost_bm = zeros(n_runs)

for i in 1:n_runs
    cost_matrix = rand(19, 19)
    cost_matrix = cost_matrix + cost_matrix'

    time_i_sdp = @elapsed Z_i_sdp, prob_sdp = solve_double_rotation_sdp(cost_matrix, use_column_orthog, use_handedness, COSMO.Optimizer)  # Mosek.Optimizer)
    Z_eigs_sdp = eigvals(Z_i_sdp)
    global eig_ratio_sdp[i] = Z_eigs_sdp[end]/Z_eigs_sdp[end-1]
    global runtime_sdp[i] = time_i_sdp
    R1_sdp, R2_sdp = extract_double_rotation(Z_i_sdp)
    global orthog_resid_sdp[i] = max(norm(I(3) - R1_sdp*R1_sdp'), norm(I(3) - R2_sdp*R2_sdp'))
    global determinant_sdp[i] = 0.5*(det(R1_sdp) + det(R2_sdp))
    global cost_sdp[i] = prob_sdp.optval

    time_i_bm = @elapsed Y_i_bm, bm_model = low_rank_double_rotation_qcqp(cost_matrix, bm_rank, use_column_orthog, use_handedness)
    Z_i_bm = Y_i_bm*Y_i_bm'
    Z_eigs_bm = eigvals(Z_i_bm)
    global eig_ratio_bm[i] = Z_eigs_bm[end]/Z_eigs_bm[end-1]
    global runtime_bm[i] = time_i_bm
    R1_bm, R2_bm = extract_double_rotation(Z_i_bm)
    global orthog_resid_bm[i] = max(norm(I(3) - R1_bm*R1_bm'), norm(I(3) - R2_bm*R2_bm'))
    global determinant_bm[i] = 0.5*(det(R1_bm) + det(R2_bm))
    global cost_bm[i] = objective_value(bm_model)
end

# TODO: remove verbose, get number of iterations, compare with IPOPT solver of low-rank Burer-Monteiro relaxation
df_sdp = DataFrame(runtime=runtime_sdp, eig_ratio=eig_ratio_sdp, orthog_resid=orthog_resid_sdp, determinant=determinant_sdp, cost=cost_sdp, method="SDP")
df_bm = DataFrame(runtime=runtime_bm, eig_ratio=eig_ratio_bm, orthog_resid=orthog_resid_bm, determinant=determinant_bm, cost=cost_bm, method="BM"*string(bm_rank))

