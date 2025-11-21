using Optim

include("inertial_point_mass_placement_analysis.jl")

# Define the P_cube base
p0 = zeros(3)
p1 = [1; 0; 0]
p2 = [0; 1; 0]
p3 = [0; 0; 1]
P_base = [p0 p1 p2 p3]
p4 = [1; 1; 0]
p5 = [1; 0; 1]
p6 = [0; 1; 1]
p7 = [1; 1; 1]
P_cube = [P_base p4 p5 p6 p7]

# Define our model
n = 10
dim = 3
n_vars = n*dim
# P_init = rand(n_vars) # zeros(n_vars)
P_init_extra = [0.5 0.5;
                0.0 1.0;
                0.5 0.5]
P_init_matrix = [P_cube P_init_extra]

# Need to randomly wiggle the cube away from its degenerate (rank-9) starting point!
P_init_matrix = P_init_matrix*0.5 .- 0.2 + rand(dim, n)*0.1
# P_init_matrix = rand(dim, n)
P_init = reshape(P_init_matrix, :, 1)
lower = -ones(size(P_init))
upper =  ones(size(P_init))

function optim_point_placement_cost(P_vec)
    P = reshape(P_vec, 3, :)
    M = points_to_data_matrix(P)
    return -logdet(M)  # Assumes optimizer is minimizing
end

function optim_point_placement_cost_svd(P_vec)
    P = reshape(P_vec, 3, :)
    M = points_to_data_matrix(P)
    return -svdvals(M)[end]  # Assumes optimizer is minimizing
end

inner_optimizer = BFGS()
# res = optimize(optim_point_placement_cost, lower, upper, P_init, Fminbox(inner_optimizer))
res = optimize(optim_point_placement_cost_svd, lower, upper, P_init, Fminbox(inner_optimizer))

P_res = reshape(res.minimizer, 3, :)
M_res = points_to_data_matrix(P_res)

using Plots, GR
# backend(:gr)
plotlyjs()  # Gives an interactive plot!
plot_points = scatter3d(P_res[1, :], P_res[2, :], P_res[3, :])
scatter3d!(P_init_matrix[1, :], P_init_matrix[2, :], P_init_matrix[3, :], color=:red)

