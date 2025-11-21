using LinearAlgebra

function point_to_data_matrix_column(p::Vector)
    p_affine = [1; p]
    inertia_vec = [p[2]^2 + p[3]^2;
                   p[1]^2 + p[3]^2;
                   p[1]^2 + p[2]^2;
                   -p[1]*p[2];
                   -p[1]*p[3];
                   -p[2]*p[3]]
    return [p_affine; inertia_vec]
end

function points_to_data_matrix(P::Matrix)
    return hcat([point_to_data_matrix_column(P[:, i]) for i=1:size(P, 2)]...)
end

function point_to_data_matrix_column_jacobian(p::Vector)
    x, y, z = p
    d_inertia = [0   2*y 2*z;
                 2*x 0   2*z;
                 2*x 2*y 0;
                 -y  -x  0;
                 -z  0   -x;
                 0   -z  -y]
    return [I(3); d_inertia]
end

function construct_full_rank_data_matrix_from_base(P::Matrix)
    n = size(P, 2)
    M = points_to_data_matrix(P)
    M = [M zeros(10, 10-n)]
    P = [P zeros(3, 10-n)]
    p_base = P[:, 2]  # Always use the same base point for now
    for i = n+1:10
        Φ_i = M[2:end, i-1]
        Mi = M[2:end, 1:i-1]
        jac = point_to_data_matrix_column_jacobian(p_base)
        dp = -pinv(Mi'*jac)*Mi'*Φ_i
        p_new = p_base + dp 
        M_test = point_to_data_matrix_column(p_new)
        # println("Iteration "*string(i))
        # println("Delta p: ")
        # println(dp)
        # println("Linearly approximated nullspace element: ")
        # println(M_test)
        # println("Actual nullspace: ")
        # println(nullspace(Mi'))
        # println("Residual check: ")
        # println(Mi'*M_test[2:end])
        p_new = p_new/norm(p_new)  # Prevents blow-up 
        P[:, i] = p_new
        M[:, i] = point_to_data_matrix_column(p_new)
    end
    return P, M
end

# p0 = zeros(3)
# p1 = [1; 0; 0]
# p2 = [0; 1; 0]
# p3 = [0; 0; 1]
# P_base = [p0 p1 p2 p3]
# M_base = points_to_data_matrix(P_base)

# # See what happens when you just add a linear scaling of one of the base points
# n_lin = 7
# P_lin = hcat([i*p3 for i=2:(n_lin+1)]...)
# M_lin = points_to_data_matrix([P_base P_lin])

# # Test a random placement
# n_rand = 10
# P_rand = rand(3, n_rand)
# M_rand = points_to_data_matrix(P_rand)

# # Test a cube - only has a rank of 7! That's a good clue.
# p4 = [1; 1; 0]
# p5 = [1; 0; 1]
# p6 = [0; 1; 1]
# p7 = [1; 1; 1]

# P_cube = [P_base p4 p5 p6 p7]
# M_cube = points_to_data_matrix(P_cube)

# # Try out my perturbation algorithm
# P_alg, M_alg = construct_full_rank_data_matrix_from_base(P_base)


# # # Plot the perturbations results 
# using Plots, GR
# backend(:gr)

# plot_points = scatter3d(P_alg[1, :], P_alg[2, :], P_alg[3, :])
# scatter3d!(P_cube[1, :], P_cube[2, :], P_cube[3, :], color=:red)

