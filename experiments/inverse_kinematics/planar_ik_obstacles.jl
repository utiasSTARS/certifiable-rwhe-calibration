using LinearAlgebra
using Plots
using Convex, MosekTools, COSMO


# Try the "tighter" rank-d convex relaxation
function equality_constraints(link_lengths, goal)
    d = length(goal)
    n = length(link_lengths)
    A = zeros(n-1+d, n-1+d, n)
    b = zeros(n)
    A[1, 1, 1] = 1
    b[1] = link_lengths[1]^2
    for i in 2:n-1
        A[i, i, i] = 1
        A[i-1, i-1, i] = 1
        A[i, i-1, i] = -1
        A[i-1, i, i] = -1
        b[i] = link_lengths[i]^2 
    end
    A[n-1, n-1, end] = 1
    A[n-1, end-d+1:end, end] = -goal'
    A[end-d+1:end, n-1, end] = -goal
    b[end] = link_lengths[end]^2 - norm(goal)^2
    return A, b
end

function obstacle_constraints(P, radii, n, d)
    m = size(P, 2)
    A = zeros(n-1+d, n-1+d, (n-1)*m)
    b = zeros((n-1)*m)
    k = 1
    for i = 1:m
        p = P[:, i]
        r = radii[i]
        for j = 1:n-1
            A[j, j, k] = 1
            A[j, end-d+1:end, k] = -p'
            A[end-d+1:end, j, k] = -p
            b[k] = r^2 - norm(p)^2
            k += 1
        end
    end
    return A, b
end

function nearest_neighbour_cost(P)
    d = size(P, 1)
    n = size(P, 2)
    C = zeros(n + d, n + d)
    for i in 1:n
        p = P[:, i]
        C[i, i] = 1
        C[i, end-d+1:end] = -p
        C[end-d+1:end, i] = -p
        C[end, end] = norm(p)^2
    end
    return C
end

function solve_sdp(A_eq, b_eq, A_ineq, b_ineq, cost, d, solver=Mosek.Optimizer)
    Z = Semidefinite(size(A_eq, 1))
    constraints_ineq = [tr(A_ineq[:, :, i]'*Z) >= b_ineq[i] for i in 1:size(A_ineq)[3]]
    constraints_eq = [tr(A_eq[:, :, i]'*Z) == b_eq[i] for i in 1:size(A_eq)[3]]
    constraint_homog = Z[end-d+1:end, end-d+1:end] == I(d)
    constraints = constraints_eq + constraints_ineq
    push!(constraints, constraint_homog)
    prob = minimize(tr(cost'*Z), constraints)
    solve!(prob, solver)
    return evaluate(Z), prob
end

function iterative_nearest_point(A_eq, b_eq, A_ineq, b_ineq, d, n, solver=Mosek.Optimizer, max_iter=10, eig_ratio_thresh=1e8)
    
    P_init = rand(d, n-1)
    Z = zeros(n-1+d, n-1+d, max_iter)
    C = nearest_neighbour_cost(P_init)
    for i in 1:max_iter
        Z_i, _ = solve_sdp(A_eq, b_eq, A_ineq, b_ineq, C, d, solver)
        eig_Z_i = eigvals(Z_i)
        Z[:, :, i] = Z_i
        if abs(eig_Z_i[d+1]/eig_Z_i[d]) > eig_ratio_thresh
            return Z[:, :, 1:i]
        end
        P_i = Z_i[end-d+1:end, 1:n-1]
        C = nearest_neighbour_cost(P_i)
    end
    return Z
end

function plot_planar_manipulator(P_rob, goal, P_obs, obs_radii)

    P = [zeros(2, 1) P_rob goal]
    p = plot(P[1, :], P[2, :], linewidth=5, linemarker=:circle, aspect_ratio=:equal)
    scatter!(P_rob[1, :], P_rob[2, :], markersize=10)

end

n = 3
d = 2
# P_init = [1 1;
#           0 1]
P_init = [sqrt(2)/2 sqrt(2);
          sqrt(2)/2 sqrt(2)]

goal = [1; 2]
P_obs = [1 0;
         0 1]
A, b = equality_constraints([1; 1; 1], goal)
radii_obs = [0.5; 0.5]
A_obs, b_obs = obstacle_constraints(P_obs, radii_obs, n, d)

Z_iterative = iterative_nearest_point(A, b, A_obs, b_obs, d, n)

P_final = Z_iterative[end-d+1:end, 1:end-d, end]

plot_planar_manipulator(P_final, goal, P_obs, radii_obs)

# C = nearest_neighbour_cost(P_init)
# Z, prob = solve_sdp(A, b, A_obs, b_obs, C, d)

# # Investigate the solution
# eigvals(Z)
# P_sol = Z[end-d+1:end, 1:n-1]
