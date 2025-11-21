include("../src/registration/correspondence_free_cost.jl")
include("../src/permutation_sdp_constraints.jl")
include("../src/rotation_sdp_solver.jl")
include("../src/low_rank_solver.jl")
include("../src/registration/permutation_utilities.jl")
using DataFrames, Statistics, Rotations

# High-level parameters
m = 6
n = 4
dim = 3
sd_noise = 0.00
sd_noises = [0.00, 0.01, 0.05]
n_noise = length(sd_noises)
n_runs = 2
rank_tol = 1e-6
stiefel = true

sdp_cost_stief = zeros(n_runs, n_noise)
primal_cost = zeros(n_runs, n_noise)
rank_stief = zeros(n_runs, n_noise)

# Ground-truth permutation and rotation
X_gt = random_injection(m, n)
R_gt = rand(RotMatrix{dim})
t_gt = rand(dim)

P_gt = rand(dim, n)
Q_gt = R_gt*P_gt*X_gt' .+ t_gt

# TODO: figure out the "false" indices later
# Q_gt[:, false_ind] .= Q_gt[:, false_ind] + rand(dim)
Q_noise = randn(size(Q_gt))*sd_noise
Q = Q_gt + Q_noise 

# Initialize the problem
A_eq_base, b_eq_base, A_ineq_base, b_ineq_base = injective_matching_with_rotation_constraints(m, n, dim, false, false)
A_eq_stief, b_eq_stief, A_ineq_stief, b_ineq_stief = injective_matching_with_rotation_constraints(m, n, dim, false, true)
cost, cost_full = correspondence_free_registration_cost_pose_schur(P_gt, Q)

num_constr_eq = size(A_eq_stief, 1)
num_constr_ineq = size(A_ineq_stief, 1)
num_constr = num_constr_eq + num_constr_ineq

rank_cifuentes = cifuentes_rank_bound(num_constr)

# sol_base, model_base = low_rank_qcqp(cost, A_eq_base, b_eq_base, A_ineq_base, b_ineq_base, rank_cifuentes)
# sol_stief, model_stief = low_rank_qcqp(cost, A_eq_stief, b_eq_stief, A_ineq_stief, b_ineq_stief, rank_cifuentes)

# sol9, _ = low_rank_qcqp(cost_stief, A_eq_stief, b_eq_stief, A_ineq_stief, b_ineq_stief, 9)


prob, cost_var, Z_var = initialize_constrained_sdp(cost, A_eq_stief, b_eq_stief, A_ineq_stief, b_ineq_stief)
# Z_sol_stief, _ = solve_constrained_sdp(cost_stief, A_eq_stief, b_eq_stief, A_ineq_stief, b_ineq_stief, Mosek.Optimizer)

# # Solve
Z_sol_stief = solve_problem(prob, Z_var, cost_var,  cost, Mosek.Optimizer)

# z_sol = Z_sol_stief[:, end]
# # z_sol = z_sol/z_sol[end]
# t_sol = get_translation_solution(cost_full, z_sol, dim)
# x_sol = z_sol[1:m*n]
# X_sol = reshape(x_sol, (m, n))
# r_sol = z_sol[m*n+1:end-1]
# R_sol = reshape(r_sol, (3, 3))

# eigvals(Z_sol_stief)

# for i in 1:n_runs
#     for j in 1:n_noise
#         Q_noise = randn(size(P_gt))*sd_noise[j]
#         Q = Q_gt + Q_noise 

#         # Get cost
#         cost, cost_stief = correspondence_free_registration_cost_pose_schur(P_gt, Q)
#         Z_sol_stief, _ = solve_constrained_sdp(cost, A_eq_stief, b_eq_stief, A_ineq_stief, b_ineq_stief, COSMO.Optimizer)

#         sdp_cost_base[i, j] = tr(Z_sol_base*cost)
#         # sdp_cost_ortho[i] = tr(Z_sol_ortho*cost)
#         # sdp_cost_ds[i] = tr(Z_sol_ds*cost)
#         sdp_cost_stief[i, j] = tr(Z_sol_stief*cost)

#         rank_base[i, j] = rank(Z_sol_base, rank_tol)
#         rank_stief[i, j] = rank(Z_sol_stief, rank_tol)

#         # TODO: implement translation as well for proper duality gap
#         # R_opt = solve_known_correspondence(X_gt, P_gt, Q)
#         # r_opt = reshape(R_opt, (dim^2,))
#         # z_opt = [x_gt; r_opt; 1]
#         # Z_opt = z_opt*z_opt'

#         # primal_cost[i, j] = tr(Z_opt*cost)
#     end
# end



# Plot the results
# using Plots, GR
# gr()

# gap_base = primal_cost .- sdp_cost_base
# gap_stief = primal_cost .- sdp_cost_stief

# # Plots.scatter(sd_noise, gap_base, color=:red, label="Base SDP")
# # Plots.scatter!(sd_noise, gap_stief, color=:blue, label="Redundant SDP")
# Plots.plot(sd_noise, maximum(gap_base', dims=2), color=:red, label="Base SDP")
# Plots.plot!(sd_noise, maximum(gap_stief', dims=2), color=:blue, label="Redundant SDP")
# xlabel!("Std. Dev. of Gaussian Noise [m]")
# ylabel!("Relaxation Gap [m^2]")
# title!("Base and Redundant SDP Relaxation")
