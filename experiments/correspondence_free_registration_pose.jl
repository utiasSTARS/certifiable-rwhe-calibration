include("../src/registration/correspondence_free_cost.jl")
include("../src/permutation_sdp_constraints.jl")
include("../src/rotation_sdp_solver.jl")
using DataFrames, Statistics

function solve_known_correspondence(X_gt, P, Q)
    Qc = Q*X_gt
    A = P*Qc'
    U, _, V = svd(A)
    σ = ones(size(P, 1))
    Σ = diagm(σ)
    Σ[end] = det(U*V')
    Rt = U*Σ*V'
    return Rt'
end

# Used to test the expression used in solve_known_correspondence() above
function registration_cost_matrix(P, Q)
    n = size(P, 2)
    dim = size(P, 1)
    A = zeros(dim, dim)
    for i in 1:n
        A = A + P[:, i]*Q[:, i]'
    end
    return A
end

# High-level parameters
n = 5
dim = 3
sd_noise = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05]
n_noise = length(sd_noise)
n_runs = 3
rank_tol = 1e-6
sdp_cost_base = zeros(n_runs, n_noise)
# sdp_cost_ortho = zeros(n_runs)
# sdp_cost_ds = zeros(n_runs)
sdp_cost_full = zeros(n_runs, n_noise)
primal_cost = zeros(n_runs, n_noise)

rank_base = zeros(n_runs, n_noise)
rank_full = zeros(n_runs, n_noise)

# Ground-truth permutation and rotation
X_gt = [0 1 0 0 0;
        0 0 0 0 1;
        1 0 0 0 0;
        0 0 1 0 0;
        0 0 0 1 0]
# X_gt = I(n)

if dim == 2
    # R_gt = [0 1;
    #        -1 0]
    R_gt = [sqrt(2)/2 sqrt(2)/2;
            -sqrt(2)/2 sqrt(2)/2]
else
    # R_gt = [0 1 0;
    #        -1 0 0;
    #         0 0 1]
    R_gt = [sqrt(2)/2 sqrt(2)/2 0;
           -sqrt(2)/2 sqrt(2)/2 0;
            0         0         1]
end

t_gt = rand(dim)

# Noise-free P and Q 
# P_gt = [0 1 0 0;
#         0 0 1 0;
#         0 0 0 1]
# P_gt = [0 1 0 0 -1;
#         0 0 2 0 -1;
#         0 0 0 3 -1]
P_gt = rand(dim, n)
Q_gt = R_gt*P_gt*X_gt' .+ t_gt


# Setup different formulations
A_eq_base, b_eq_base, A_ineq_base = permutation_with_rotation_constraints(n, dim, true, false, false)
b_ineq_base = zeros(size(A_ineq_base, 1))

A_eq_ortho, b_eq_ortho, A_ineq_ortho = permutation_with_rotation_constraints(n, dim, true, true, false)
b_ineq_ortho = zeros(size(A_ineq_ortho, 1))

A_eq_ds, b_eq_ds, A_ineq_ds = permutation_with_rotation_constraints(n, dim, true, false, true)
b_ineq_ds = zeros(size(A_ineq_ds, 1))

A_eq_full, b_eq_full, A_ineq_full = permutation_with_rotation_constraints(n, dim, true, true, true)
b_ineq_full = zeros(size(A_ineq_full, 1))

for i in 1:n_runs
    for j in 1:n_noise
        Q_noise = randn(size(P_gt))*sd_noise[j]
        Q = Q_gt + Q_noise 

        # Get cost
        cost = correspondence_free_registration_cost_pose(P_gt, Q)

        # Solve for each constraint set
        Z_sol_base, _ = solve_constrained_sdp(cost, A_eq_base, b_eq_base, A_ineq_base, b_ineq_base, Mosek.Optimizer)
        # Z_sol_ortho, _ = solve_constrained_sdp(cost, A_eq_ortho, b_eq_ortho, A_ineq_ortho, b_ineq_ortho, Mosek.Optimizer)
        # Z_sol_ds, _ = solve_constrained_sdp(cost, A_eq_ds, b_eq_ds, A_ineq_ds, b_ineq_ds, Mosek.Optimizer)
        Z_sol_full, _ = solve_constrained_sdp(cost, A_eq_full, b_eq_full, A_ineq_full, b_ineq_full, Mosek.Optimizer)

        r_gt = reshape(R_gt, (dim^2,))
        x_gt = reshape(X_gt, (n^2,))
        z_gt = [t_gt; x_gt; r_gt; 1]
        Z_gt = z_gt*z_gt'

        sdp_cost_base[i, j] = tr(Z_sol_base*cost)
        # sdp_cost_ortho[i] = tr(Z_sol_ortho*cost)
        # sdp_cost_ds[i] = tr(Z_sol_ds*cost)
        sdp_cost_full[i, j] = tr(Z_sol_full*cost)

        rank_base[i, j] = rank(Z_sol_base, rank_tol)
        rank_full[i, j] = rank(Z_sol_full, rank_tol)

        # TODO: implement translation as well for proper duality gap
        # R_opt = solve_known_correspondence(X_gt, P_gt, Q)
        # r_opt = reshape(R_opt, (dim^2,))
        # z_opt = [x_gt; r_opt; 1]
        # Z_opt = z_opt*z_opt'

        # primal_cost[i, j] = tr(Z_opt*cost)
    end
end



# Plot the results
# using Plots, GR
# gr()

# gap_base = primal_cost .- sdp_cost_base
# gap_full = primal_cost .- sdp_cost_full

# # Plots.scatter(sd_noise, gap_base, color=:red, label="Base SDP")
# # Plots.scatter!(sd_noise, gap_full, color=:blue, label="Redundant SDP")
# Plots.plot(sd_noise, maximum(gap_base', dims=2), color=:red, label="Base SDP")
# Plots.plot!(sd_noise, maximum(gap_full', dims=2), color=:blue, label="Redundant SDP")
# xlabel!("Std. Dev. of Gaussian Noise [m]")
# ylabel!("Relaxation Gap [m^2]")
# title!("Base and Redundant SDP Relaxation")
