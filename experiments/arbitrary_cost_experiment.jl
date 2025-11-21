using DataFrames
using COSMO
using COSMOAccelerators
using MosekTools
using LinearAlgebra
using Statistics
include("../src/rotation_sdp_solver.jl")
include("../src/low_rank_solver.jl")

# TODO: test COSMOAccelerators
# settings = COSMO.Settings(accelerator = with_options(AndersonAccelerator, mem = 15))
# settings = COSMO.Settings(safeguard = false)

# TODO: test parallelization 
# TODO: think about warm start in the DL setting (incremental updates to parameters)

# Experiment parameters
# TODO: add duality gap
n_runs = 100
use_column_orthog = true
use_handedness = true

# BM parameters
bm_rank = 5

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
    cost_matrix = rand(10, 10)
    cost_matrix = cost_matrix + cost_matrix'

    time_i_sdp = @elapsed Z_i_sdp, prob_sdp = solve_rotation_sdp(cost_matrix, use_column_orthog, use_handedness, COSMO.Optimizer)  # Mosek.Optimizer)
    Z_eigs_sdp = eigvals(Z_i_sdp)
    global eig_ratio_sdp[i] = Z_eigs_sdp[end]/Z_eigs_sdp[end-1]
    global runtime_sdp[i] = time_i_sdp
    R_sdp = extract_rotation(Z_i_sdp)
    global orthog_resid_sdp[i] = norm(I(3) - R_sdp*R_sdp')
    global determinant_sdp[i] = det(R_sdp)
    global cost_sdp[i] = prob_sdp.optval

    time_i_bm = @elapsed Y_i_bm, bm_model = low_rank_rotation_qcqp(cost_matrix, bm_rank, use_column_orthog, use_handedness)
    Z_i_bm = Y_i_bm*Y_i_bm'
    Z_eigs_bm = eigvals(Z_i_bm)
    global eig_ratio_bm[i] = Z_eigs_bm[end]/Z_eigs_bm[end-1]
    global runtime_bm[i] = time_i_bm
    R_bm = extract_rotation(Z_i_bm)
    global orthog_resid_bm[i] = norm(I(3) - R_bm*R_bm')
    global determinant_bm[i] = det(R_bm)
    global cost_bm[i] = objective_value(bm_model)
end

# TODO: remove verbose, get number of iterations, compare with IPOPT solver of low-rank Burer-Monteiro relaxation
df_sdp = DataFrame(runtime=runtime_sdp, eig_ratio=eig_ratio_sdp, orthog_resid=orthog_resid_sdp, determinant=determinant_sdp, cost=cost_sdp, method="SDP")
df_bm = DataFrame(runtime=runtime_bm, eig_ratio=eig_ratio_bm, orthog_resid=orthog_resid_bm, determinant=determinant_bm, cost=cost_bm, method="BM"*string(bm_rank))

