using Convex, LinearAlgebra, MosekTools, COSMO, DataFrames
using Plots
using LaTeXStrings
using GR
import GR.heatmap
gr()

begin
	d = 2
	n = 3
	x0 = [0.; 0.]
	goal = [1.; 2.];
end

Z = Semidefinite((n-1)*d + 1);

begin
	A0 = zeros(size(Z))
	A0[end, end] = 1
	c0 = tr(A0*Z) == 1;
end

begin
	A1 = zeros(size(Z))
	A1[1, 1] = 1
	A1[2, 2] = 1
	c1 = tr(A1*Z) == 1

	A2 = zeros(size(Z))
	A2[1, 1] = 1
	A2[2, 2] = 1
	A2[3, 3] = 1
	A2[4, 4] = 1
	A2[1, 3] = -1
	A2[3, 1] = -1
	A2[2, 4] = -1
	A2[4, 2] = -1
	c2 = tr(A2*Z) == 1

	A3 = zeros(size(Z))
	A3[3, 3] = 1
	A3[4, 4] = 1
	A3[3:4, end] = -goal
	A3[end, 3:4] = -goal'
	c3 = tr(A3*Z) == 1 - norm(goal)^2;
end

function nearest_neighbour_cost(x1, x2)
	C = Matrix(1.0I(5))
	C[end, end] = norm(x1)^2 + norm(x2)^2
	C[1:2, end] = -x1
	C[end, 1:2] = -x1
	C[3:4, end] = -x2
	C[end, 3:4] = -x2
	return C
end

begin
	x1_gt = [0; 1]
	x2_gt = [1; 1]
	x_gt = [x1_gt; x2_gt; 1]
	Z_gt = x_gt*x_gt'
	C0 = nearest_neighbour_cost(x1_gt, x2_gt)
end

begin
	eq_constraints = [c0, c1, c2, c3]
	prob = minimize(tr(C0*Z), eq_constraints)
	solve!(prob, Mosek.Optimizer)
	evaluate(Z)
end


z_rank_1 = eigvecs(evaluate(Z))[:, end]

z_sol = z_rank_1/z_rank_1[end]

function extract_solution(Z)
	z_rank_1 = eigvecs(evaluate(Z))[:, end]
	return z_rank_1/z_rank_1[end]
end

function total_constraint_violation(x, A, b)
	return sum([abs(x'*A[:, :, i]*x - b[i]) for i in 1:size(A, 3)])
end

begin
	perturb_max = 2.
	n_grid = 100
	dx1 = range(-perturb_max, stop=perturb_max, length=n_grid)
	dx2 = dx1
	
	eig_ratios = zeros(n_grid, n_grid)
	total_constraint_violations = zeros(size(eig_ratios))
	# solver_statuses = zeros(Int, size(eig_ratios))

	A = cat(A0, A1, A2, A3, dims=3)
	b = [1; 1; 1; 1 - norm(goal)^2];
	
end

begin
	for i in 1:n_grid
		for j in 1:n_grid
			x1_perturbed = x1_gt + [dx1[i]; dx2[j]]
            C = nearest_neighbour_cost(x1_perturbed, x2_gt)
			# x2_perturbed = x2_gt + [dx1[i]; dx2[j]]
			# C = nearest_neighbour_cost(x1_gt, x2_perturbed)
			prob = minimize(tr(C*Z), eq_constraints)
			solve!(prob, Mosek.Optimizer)
			x_sol = extract_solution(evaluate(Z))

			eig_ratios[i, j] = eigvals(evaluate(Z))[end]/eigvals(evaluate(Z))[end-1]
			total_constraint_violations[i, j] = total_constraint_violation(x_sol, A, b);
			# solver_statuses[i, j] = prob.status;
		end
	end
end

# Post-process to get rank
eig_ratio_thresh = 1e6
rank_results = eig_ratios .< eig_ratio_thresh

# Visualize case without inequality constraint
# Plot contour map of eig_ratios
font_family = "Computer Modern"
guide_font_size = 12
title_font_size = 13
lw = 2.5
x_extent = [-perturb_max, perturb_max]
y_extent = [-perturb_max, perturb_max]

# Base case
Plots.heatmap(dx1, dx2, rank_results', aspect_ratio=:equal, levels=2, colorbar=false, framestyle=:box, widen=false)
plot!(fontfamily=font_family)
plot!(titlefontsize=title_font_size)    
plot!(legendfontsize=guide_font_size)
title!("SDP Solution Rank (Planar IK)")
xlabel!(L"\delta x")
ylabel!(L"\delta y")
xaxis!(x_extent)
yaxis!(y_extent)
Plots.savefig("planar_ik_stability.pdf")

begin
	A_ineq = zeros(size(A0))
	A_ineq[3, 3] = 1
	A_ineq[4, 4] = 1
	c_ineq = tr(A_ineq*Z) >= 2.;
end

begin
	constraints = [c0, c1, c2, c3, c_ineq]
	prob_ineq = minimize(tr(C0*Z), constraints)
	solve!(prob_ineq, Mosek.Optimizer)
	evaluate(Z)
end

extract_solution(evaluate(Z))

eigvals(evaluate(Z))

# Experiment for inequalities
begin
	eig_ratios_ineq = zeros(n_grid, n_grid)
	total_constraint_violations_ineq = zeros(size(eig_ratios))
	for i in 1:n_grid
		for j in 1:n_grid
			x1_perturbed = x1_gt + [dx1[i]; dx2[j]]
            C = nearest_neighbour_cost(x1_perturbed, x2_gt)
			# x2_perturbed = x2_gt + [dx1[i]; dx2[j]]
			# C = nearest_neighbour_cost(x1_gt, x2_perturbed)
			prob_ineq = minimize(tr(C*Z), constraints)
			solve!(prob_ineq, Mosek.Optimizer, verbose=false)
			x_sol = extract_solution(evaluate(Z))

			eig_ratios_ineq[i, j] = eigvals(evaluate(Z))[end]/eigvals(evaluate(Z))[end-1]
			total_constraint_violations_ineq[i, j] = total_constraint_violation(x_sol, A, b) + max(0., 2.0 - x_sol'*A_ineq*x_sol);
			# solver_statuses[i, j] = prob.status;
		end
	end
end

rank_results_ineq = eig_ratios_ineq .< eig_ratio_thresh

# Plot case with inequalities
Plots.heatmap(dx1, dx2, rank_results_ineq', aspect_ratio=:equal, colorbar=false, framestyle=:box, widen=false)
plot!(fontfamily=font_family)
plot!(titlefontsize=title_font_size)    
plot!(legendfontsize=guide_font_size)
title!("SDP Solution Rank (Planar IK With Inequality)")
xlabel!(L"\delta x")
ylabel!(L"\delta y")
xaxis!(x_extent)
yaxis!(y_extent)
Plots.savefig("planar_ik_inequality_stability.pdf")
