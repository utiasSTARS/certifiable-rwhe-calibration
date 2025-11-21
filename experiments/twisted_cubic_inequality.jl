using LinearAlgebra
using Convex 
using MosekTools
using Optim
using Plots 
using LaTeXStrings
using GR
import GR.contour
gr()

function twisted_cubic(t)
    return [t; t^2; t^3]
end


function solve_primal_twisted_cubic(θ, t_min=-5.0, t_max=5.0)
    f_univariate(t) = norm(twisted_cubic(t) - θ)^2
    res = optimize(f_univariate, t_min, t_max)
    return res.minimum, res.minimizer
end

function solve_twisted_cubic_sdp_relaxation(θ, inequality=false)
    Z = Semidefinite(4)
    # Nearest-point cost function 
    C = [I(3)    -θ; 
         -θ'   θ'*θ]
    # Constraints
    A = zeros(4, 4, 4)
    # y2 = y1^2
    A[1, 1, 1] = -1 
    A[1, 2, 4] = 0.5
    A[1, 4, 2] = 0.5
    # y3 = y1*y2
    A[2, 4, 3] = 0.5
    A[2, 3, 4] = 0.5
    A[2, 1, 2] = -0.5
    A[2, 2, 1] = -0.5
    # s^2 = 1
    A[3, 4, 4] = 1

    if inequality
        # y1 >= 0
        A[4, 1, 4] = 0.5
        A[4, 4, 1] = 0.5
        constraints = [tr(A[1, :, :]'*Z) == 0, tr(A[2, :, :]'*Z) == 0, tr(A[3, :, :]'*Z) == 1, tr(A[4, :, :]'*Z) >= 0]
    else
        constraints = [tr(A[1, :, :]'*Z) == 0, tr(A[2, :, :]'*Z) == 0, tr(A[3, :, :]'*Z) == 1]
    end
        # Solve with Convex library (MOSEK as backend)
    prob = minimize(tr(C'*Z), constraints)
    solve!(prob, Mosek.Optimizer, verbose=false)
    return evaluate(Z), prob
end

# Experiment parameters
θ1_range = [-1, 1]
θ3_range = [-3, 3]
n_samples = 100
dθ1 = range(θ1_range[1], stop=θ1_range[2], length=n_samples)
dθ3 = range(θ3_range[1], stop=θ3_range[2], length=n_samples)
eig_ratios = zeros(n_samples, n_samples)
eig_ratios_ineq = zeros(n_samples, n_samples)
gap = zeros(n_samples, n_samples)
gap_ineq = zeros(n_samples, n_samples)

# Run experiment
for i in 1:length(dθ1) 
    for j in 1:length(dθ3)
        θ = [dθ1[i]; dθ1[i]^2; dθ3[j]]

        # Solve primal first
        min_val, t_min =  solve_primal_twisted_cubic(θ)
        min_val_ineq, t_min_ineq = solve_primal_twisted_cubic(θ, 0.0, 5.0)

        # Solve without inequalities
        Z, prob = solve_twisted_cubic_sdp_relaxation(θ)
        Z_eigs = eigvals(Z)
        eig_ratios[i, j] = abs(real(Z_eigs[end])/real(Z_eigs[end-1]))
        gap[i, j] = min_val - prob.optval
        
        Z_ineq, prob_ineq = solve_twisted_cubic_sdp_relaxation(θ, true)
        Z_ineq_eigs = eigvals(Z_ineq)
        eig_ratios_ineq[i, j] = abs(real(Z_ineq_eigs[end])/real(Z_ineq_eigs[end-1]))
        gap_ineq[i, j] = min_val_ineq - prob_ineq.optval

    end
end

# Plot contour map of eig_ratios
font_family = "Computer Modern"
guide_font_size = 12
title_font_size = 13
lw = 2.5
x_extent = [-1., 1.]
y_extent = [-3., 3.]
# Base case
Plots.contour(dθ1, dθ3, gap', levels=[-1e-6, 1e-3, 0.2, 0.4, 0.6, 0.8, 1.0], fill=true, legend=false, colorbar=true, framestyle=:box, widen=false)
plot!(dθ1, dθ1.^3, linecolor=:white, linewidth=lw, label="Twisted Cubic")
plot!(fontfamily=font_family)
plot!(titlefontsize=title_font_size)    
plot!(legendfontsize=guide_font_size)
title!("SDP Relaxation Gap (Twisted Cubic)")
xlabel!(L"\gamma_1")
ylabel!(L"\gamma_3")
xaxis!(x_extent)
yaxis!(y_extent)
Plots.savefig("twisted_cubic.pdf")

# Inequality case 
Plots.contour(dθ1, dθ3, gap_ineq', levels=[-1e-6, 1e-3, 0.5, 1, 2, 3, 4], fill=true, legend=false, colorbar=true, framestyle=:box, widen=false)
plot!(dθ1[Int(n_samples/2):end], dθ1[Int(n_samples/2):end].^3, linecolor=:white, linewidth=lw, label="Truncated Twisted Cubic")
plot!(fontfamily=font_family)
plot!(titlefontsize=title_font_size)    
plot!(legendfontsize=guide_font_size)
title!("SDP Relaxation Gap (Truncated Twisted Cubic)")
xlabel!(L"\gamma_1")
ylabel!(L"\gamma_3")
xaxis!(x_extent)
yaxis!(y_extent)
Plots.savefig("truncated_twisted_cubic.pdf")