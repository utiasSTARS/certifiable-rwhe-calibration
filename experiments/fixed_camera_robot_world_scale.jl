using MKL
using LinearAlgebra
using Rotations
using SparseArrays
using DelimitedFiles

include("../src/calibration/robot_world_costs.jl");
include("../src/rotation_sdp_solver_jump.jl");
include("../src/utils/rotation_noise.jl")

function generate_camera_poses(T_bi)

    #Create camera poses
    m = 4
    T_ic = zeros(m, 4, 4)
    X = zeros(m, 4, 4)
    T_ic[:, 4, 4] .= 1.0

    T_ic[:, 1:3, 4] = transpose([-1.0 0.0 -1.0 -1.0
        0.0 0.0 1.0 -1.0
        0.5 1.0 1.0 1.0])

    for j in range(1, 4)
        T_ic[j, 1:3, 3] = -normalize(T_ic[j, 1:3, 4])
        if transpose(T_ic[j, 1:3, 3]) * [0; 1; 0] == 0
            T_ic[j, 1:3, 2] = [0 1 0]
        else
            T_ic[j, 1:3, 2] = normalize([-1 0 -1])
        end
        T_ic[j, 1:3, 1] = cross(T_ic[j, 1:3, 2], T_ic[j, 1:3, 3])
        X[j, :, :] = T_bi * T_ic[j, :, :]
    end
    return T_ic, X
end

function generate_arm_camera_measurements(Y, T_bi, T_ic)

    #Create target translation
    n_grid = 3;
    m = size(T_ic, 1)
    num_points_on_line = 4 #Number of points to trace a line
    max_dist = 0.6 #Maximum distance from frame i 
    t_i_pi = hcat([[x;y;z] for x in range(-max_dist, max_dist, n_grid) for y in range(-max_dist, max_dist, n_grid) for z in range(-max_dist, max_dist, n_grid)]...)
    n= size(t_i_pi, 2)
    line_to_trace = zeros(3, num_points_on_line)
    line_to_trace[1, :] .= range(-10.0, -1.0, num_points_on_line)|>collect
    line_to_trace[2, :] = range(-10.0, 10.0, num_points_on_line) |> collect
    line_to_trace[3, :] = range(0.0, 2.0, num_points_on_line) |> collect

    T_ip = zeros(n * num_points_on_line, 4, 4)
    A = zeros(n * num_points_on_line, 4, 4)
    T_ip[:, 4, 4] .= 1.0
    for j in range(1, n)
        for k in range(1, num_points_on_line)
            T_ip[(j-1)*num_points_on_line+k, 1:3, 4] = t_i_pi[:, j]
            T_ip[(j-1)*num_points_on_line+k, 1:3, 3] = normalize(line_to_trace[:, k] - t_i_pi[:, j])
            T_ip[(j-1)*num_points_on_line+k, 1:3, 2] = normalize(cross(T_ip[(j-1)*num_points_on_line+k, 1:3, 3], diff(line_to_trace, dims=2)[:, 1]))
            T_ip[(j-1)*num_points_on_line+k, 1:3, 1] = cross(T_ip[(j-1)*num_points_on_line+k, 1:3, 2], T_ip[(j-1)*num_points_on_line+k, 1:3, 3])
            T_ip[(j-1)*num_points_on_line+k, 1:3, 1:3] *= Matrix(RotZ(2*pi*(k-1.0)/(num_points_on_line-1.0)))
            A[(j-1)*num_points_on_line+k, :, :] = Y * inv(T_bi * T_ip[(j-1)*num_points_on_line+k, :, :])
        end
    end

    B = zeros(m, n * num_points_on_line, 4, 4)
    for i in range(1, m)
        for j in range(1, n * num_points_on_line)
            B[i, j, :, :] = inv(T_ip[j, :, :]) * T_ic[i, :, :]
            #B[i, j, 1:3, 4] .*=0.5
            #display(A[j, :, :]*X[i, :, :] - Y*B[i, j, :, :])
        end
    end

    return A, B
end

function format_data_for_csv(poses::AbstractArray)
    #Format single pose
    if length(size(poses))==2
        return vcat(Rotations.params(QuatRotation(poses[1:3, 1:3])), poses[1:3, 4])'
    end

    #Format (N, 4, 4) pose matrix
    num_poses = size(poses, 1);
    pose_vector= zeros(num_poses, 7);
    for i in 1:num_poses
        pose_vector[i, :] = vcat(Rotations.params(QuatRotation(poses[i, 1:3, 1:3])), poses[i, 1:3, 4])';
    end
    return pose_vector
end

function perform_experiment()
    #Set Variables
    num_experiments = 100
    scale_factor = 0.5
    half_langevin_k = 12 / 2
    trans_std_dev = 0.01

    #Create constants
    Y = [1 0 0 0
        0 1 0 0
        0 0 1 0.1
        0 0 0 1]

    T_bi = [-1 0 0 0.5
        0 -1 0 0
        0 0 1 0
        0 0 0 1]

    T_ic, X = generate_camera_poses(T_bi)

    test_case_result = zeros(num_experiments, 11)

    for l in range(1, num_experiments)
        A, B = generate_arm_camera_measurements(Y, T_bi, T_ic)

        #Save Data
        folder_structure = string("data/fixed_camera_robot_world_scale/kappa_", trunc(Int, 2 * half_langevin_k), "_sd_trans_", trunc(Int, trans_std_dev * 100), "/exp_", l, "/")
        if ~isdir(folder_structure)
            mkpath(folder_structure)
        end
        A_save_format = format_data_for_csv(A)
        writedlm(string(folder_structure, "A.csv"), A_save_format, ',')
        Y_save_format = format_data_for_csv(Y)
        writedlm(string(folder_structure, "Y.csv"), Y_save_format, ',')

        n = size(B, 1)
        num_meas = size(B, 2)

        for i in range(1, n)
            for j in range(1, num_meas)
                B[i, j, 1:3, 1:3] *= langevin_sample(half_langevin_k)
                B[i, j, 1:3, 4] += trans_std_dev * randn() .* normalize(rand(Float64, 3))
                B[i, j, 1:3, 4] .*= scale_factor;
                B_save_format = format_data_for_csv(B[i, :, :, :])
                writedlm(string(folder_structure, "B", i, ".csv"), B_save_format, ',')
                X_save_format = format_data_for_csv(X[i, :, :])
                writedlm(string(folder_structure, "X", i, ".csv"), X_save_format, ',')
            end
        end

        display("Generated Data")
        Q = get_empty_sparse_cost_matrix(n, 1, true)
        for i = 1:n
            Q_temp = sparse_robot_world_transformation_scale_cost(A, B[i, :, :, :], half_langevin_k .* ones(num_meas), (0.5/(trans_std_dev^2)) .* ones(num_meas))
            Q += get_ith_eye_jth_base_cost(i, 1, n, 1, Q_temp, true)
        end

        Z_schur, model = solve_sdp_dual_w_scale_using_schur(Q, n + 1, true, true)
        solution = extract_solution_from_dual_schur(Z_schur, Q)
        error = zeros(11)
        for i in range(1, n + 1)
            if i < n + 1
                error[i] = norm(solution[3*i-2:3*i]/solution[3*(n + 1) + 1] - X[i, 1:3, 4])
                error[i+n+2] = norm(inv(Rotations.ExponentialMap())(nearest_rotation(reshape(solution[9*(i-1)+3*(n+1)+2:9*i+3*(n+1)+1], 3, 3) * X[i, 1:3, 1:3]')))
            else
                error[i] = norm(solution[3*i-2:3*i]/solution[3*(n + 1) + 1] - Y[1:3, 4])
                error[i+n+2] = norm(inv(Rotations.ExponentialMap())(nearest_rotation(reshape(solution[9*(i-1)+3*(n+1)+2:9*i+3*(n+1)+1], 3, 3) * Y[1:3, 1:3]')))
            end
        end
        error[n+2] = abs(solution[3*(n + 1) + 1] - scale_factor)
        display(solution[3*(n + 1) + 1])
        display(transpose(error))
        test_case_result[l, :] = error

    end
    writedlm(string("data/fixed_camera_robot_world_scale/kappa_", trunc(Int, 2 * half_langevin_k), "_sd_trans_", trunc(Int, trans_std_dev * 100), "/convex_results.csv"), test_case_result, ',')
end

perform_experiment();