using LinearAlgebra
using Plots
using Rotations
using NPZ
using DelimitedFiles

include("../src/calibration/robot_world_costs.jl");
include("../src/rotation_sdp_solver_jump.jl");
include("../src/utils/rotation_noise.jl")

#B-Spline Functions
function mixing_matrix(k::Int64)
    #Computes the kth order mixing matrix for a Uniform B-Spline
    M = zeros(k, k);
    for i in 1:k
        for j in 1:k
            i_temp = i-1;
            j_temp = j-1;
            Cin = binomial(k-1, k-1-i_temp);
            for s=j_temp:k-1
                M[j,i] += (-1)^(s-j_temp)*binomial(k, s-j_temp)*(k-s-1)^(k-1-i_temp);
            end
            M[j, i] *= 1/factorial(k-1)*Cin;
        end
    end
    return M    
end

function u2vec(u::Float64, k::Int64, n::Int64)
    #Computes the nth derivative u vector for a Uniform B-Spline of order k 
    uvec = zeros(k);
    uvec[n + 1] = 1;
    for i in n+2:k
        uvec[i] = uvec[i-1]*u;
    end
    return uvec
end

function get_spline_knot_ids(dt::Float64, t::Float64, k::Int64)
    ids = floor(t/dt) + 1:k + floor(t/dt);
    u = (t - floor(t/dt)*dt)/dt;
    return ids, u
end

function bsplineeval(knots, knot_interval, t, k, n)
    #k is the spline order
    #n is the order of the derivative
    #t is the time to evaluate the spline
    ids, u = get_spline_knot_ids(knot_interval, t, k);
    ids = convert(Vector{Int64}, ids)
    M = mixing_matrix(k);
    uvec = u2vec(u, k, n);
    p = (1.0/knot_interval)^n*knots[:, ids]*M*uvec;
    return p
end

function angles_to_T_ic(theta::Float64, phi::Float64, r::Float64)
return [-cos(theta)*cos(phi) sin(theta) -cos(theta)*sin(phi) r*cos(theta)*sin(phi);
        sin(phi) 0 -cos(phi) r*cos(phi);
        -sin(theta)*cos(phi) -cos(theta) -sin(theta)*sin(phi) r*sin(theta)*sin(phi);
        0 0 0 1];
end

function angles_to_rot_vel_c_ci(phi::Float64, theta_dot::Float64, phi_dot::Float64)
    return [-sin(phi) 0;
            0 1;
            cos(phi) 0]*[theta_dot;phi_dot]
end

function angles_to_rot_accel_c_ci(phi::Float64, theta_dot::Float64, phi_dot::Float64, theta_2dot::Float64, phi_2dot::Float64)
    return [-cos(phi) -sin(phi) 0;
            0 0 1;
            -sin(phi) cos(phi) 0]*[theta_dot*phi_dot;theta_2dot;phi_2dot]
end

function generate_camera_data(num_meas::Int64,
    knot_interval::Float64,
    k::Int64,
    max_theta_val::Float64,
    min_theta_val::Float64,
    max_phi_val::Float64,
    min_phi_val::Float64,
    radius_of_sphere::Float64,
    visualize::Bool)

    #Settings
    sampling_interval = 0.5*knot_interval;
    max_time = num_meas*sampling_interval;
    num_knots = convert(Int64, max_time÷knot_interval) + k;
    eval_times = 0:sampling_interval:max_time-sampling_interval; 
    theta_knots = rand(min_theta_val:(max_theta_val-min_theta_val)/1000:max_theta_val, num_knots);
    phi_knots = rand(min_phi_val:(max_phi_val-min_phi_val)/1000:max_phi_val, num_knots);

    #Compute the splines and associated derivatives
    theta_spline = [bsplineeval(theta_knots', knot_interval, time, k, 0) for time in eval_times];
    phi_spline = [bsplineeval(phi_knots', knot_interval, time, k, 0) for time in eval_times];
    theta_dot_spline = [bsplineeval(theta_knots', knot_interval, time, k, 1) for time in eval_times];
    phi_dot_spline = [bsplineeval(phi_knots', knot_interval, time, k, 1) for time in eval_times];
    theta_2dot_spline = [bsplineeval(theta_knots', knot_interval, time, k, 2) for time in eval_times];
    phi_2dot_spline = [bsplineeval(phi_knots', knot_interval, time, k, 2) for time in eval_times];

    #Poses, Velocities, and accelerations
    T_ic = [angles_to_T_ic(θ, ϕ, radius_of_sphere) for (θ,ϕ) in zip(theta_spline,phi_spline)];
    ω_c_ci = hcat([angles_to_rot_vel_c_ci(ϕ, θ_dot, ϕ_dot) for (ϕ,θ_dot,ϕ_dot) in zip(phi_spline, theta_dot_spline, phi_dot_spline)]...);
    α_c_ci = hcat([angles_to_rot_accel_c_ci(ϕ, θ_dot, ϕ_dot, θ_2dot, ϕ_2dot) for (ϕ,θ_dot,ϕ_dot,θ_2dot,ϕ_2dot) in zip(phi_spline, theta_dot_spline, phi_dot_spline, theta_2dot_spline, phi_2dot_spline)]...);

    if visualize
        plotting_interval = sampling_interval/10;
        plotting_times = 0:plotting_interval:max_time-plotting_interval;
        theta_plotting_spline = [bsplineeval(theta_knots', knot_interval, time, k, 0) for time in plotting_times];
        phi_plotting_spline = [bsplineeval(phi_knots', knot_interval, time, k, 0) for time in plotting_times];
        T_ic_plotting = [angles_to_T_ic(θ, ϕ, radius_of_sphere) for (θ,ϕ) in zip(theta_plotting_spline,phi_plotting_spline)];
        r_i_ci = hcat([pose[1:3, end] for pose in T_ic_plotting]...);
        display(plot(r_i_ci[1, :], r_i_ci[2, :], r_i_ci[3, :]))
    end
    return T_ic, ω_c_ci, α_c_ci
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
    #Camera Trajectory Generator Settings
    num_experiments = 100
    num_meas=100;
    knot_interval=1.0;
    k=4;
    max_theta_val=3*pi/4;
    min_theta_val=pi/4;
    max_phi_val=3*pi/4;
    min_phi_val=pi/4;
    radius_of_sphere = 1.0;
    visualize = false;
    half_langevin_k = 12/2;
    trans_std_dev = 0.01;
    scale_factor = 0.5; 

    #Results storage
    test_case_result = zeros(num_experiments, 5)

    for i in 1:num_experiments
        #Generate X and Y

        T_bi = [-1.0 0 0 0.5;
                 0 -1.0 0 0;
                 0 0 1 0;
                 zeros(3)' 1]; #Y
        T_hc = [1.0 0 0 0;
                0 1.0 0 0;
                0 0 1 0.1;
                zeros(3)' 1]; #X

        T_ic, ω_c_ci, α_c_ci = generate_camera_data(Int64(num_meas/2),
                                                    knot_interval,
                                                    k,
                                                    max_theta_val,
                                                    min_theta_val,
                                                    max_phi_val,
                                                    min_phi_val,
                                                    radius_of_sphere,
                                                    visualize);
        
        T_ic_2, ω_c_ci_2, α_c_ci_2 = generate_camera_data(Int64(num_meas/2),
                                                    knot_interval,
                                                    k,
                                                    max_theta_val,
                                                    min_theta_val,
                                                    max_phi_val,
                                                    min_phi_val,
                                                    0.3* radius_of_sphere,
                                                    visualize);
        
        T_ic = vcat(T_ic, T_ic_2)
        #T_ic = [[RotXYZ(2*pi*rand(), 2*pi*rand(),2*pi*rand()) [1;1;0.5].*rand(3, 1)-[0.5;0.5;0]; 0 0 0 1] for i in 1:num_meas]
        T_bh = [T_bi*pose*inv(T_hc) for pose in T_ic];

        #Format Data
        A = zeros(num_meas, 4, 4);
        B = zeros(num_meas, 4, 4);
        for j in 1:num_meas
            A[j, :,:] = T_bh[j];
            B_temp = T_ic[j];
            B_temp[1:3, 1:3] *=langevin_sample(half_langevin_k);
            noise_vec = rand(Float64, 3);
            noise_vec .*= trans_std_dev * randn()/(transpose(noise_vec)*noise_vec); 
            B_temp[1:3, 4] += noise_vec;
            B_temp[1:3, 4] .*= scale_factor;
            B[j, :,:] = B_temp;
        end

        #Save Data
        folder_structure = string("data/robot_world_poses_on_sphere_scale/kappa_", trunc(Int, 2*half_langevin_k), "_sd_trans_", trunc(Int, trans_std_dev*100) ,"/exp_", i, "/")
        if ~isdir(folder_structure)
            mkpath(folder_structure)
        end
        A_save_format = format_data_for_csv(A)
        writedlm(string(folder_structure,"A.csv"),  A_save_format, ',')
        B_save_format = format_data_for_csv(B)
        writedlm(string(folder_structure,"B.csv"),  B_save_format, ',')
        X_save_format = format_data_for_csv(T_hc)
        writedlm(string(folder_structure,"X.csv"),  X_save_format, ',')
        Y_save_format = format_data_for_csv(T_bi)
        writedlm(string(folder_structure,"Y.csv"),  Y_save_format, ',')


        #Get Cost
        Q = sparse_robot_world_transformation_scale_cost(A, B, half_langevin_k .* ones(num_meas), (0.5/(trans_std_dev^2)) .* ones(num_meas))
        #Q = sparse_robot_world_transformation_cost(A, B, ones(num_meas),  ones(num_meas))

        #Solve the problem
        #Z, model = solve_sdp_dual_w_scale_using_concat(Q, 2, true, true)
        Z, model = solve_sdp_dual_w_scale_using_schur(Q, 2, true, true)
        #Z, model = solve_sdp_dual_using_concat(Q, 2, true, true)
        solution = extract_solution_from_dual_schur(Z, Q);
        #solution = extract_solution_from_dual(Z);

        error = [norm(solution[1:3]/solution[7] - T_hc[1:3, 4]); 
        norm(solution[4:6]/solution[7] - T_bi[1:3, 4]);
        abs(solution[7] - scale_factor); 
        norm(inv(Rotations.ExponentialMap())(nearest_rotation(reshape(solution[8:16], 3, 3) * T_hc[1:3, 1:3]'))); 
        norm(inv(Rotations.ExponentialMap())(nearest_rotation(reshape(solution[17:25], 3, 3) * T_bi[1:3, 1:3]')))];

        #error = [solution[1:3] - T_hc[1:3, 4]; 
        #solution[4:6] - T_bi[1:3, 4];
        #inv(Rotations.ExponentialMap())(nearest_rotation(reshape(solution[7:15], 3, 3) * T_hc[1:3, 1:3]')); 
        #inv(Rotations.ExponentialMap())(nearest_rotation(reshape(solution[16:24], 3, 3) * T_bi[1:3, 1:3]'))]

        display(transpose(error))
        test_case_result[i, :] = error;
    end

    writedlm(string("data/robot_world_poses_on_sphere_scale/kappa_", trunc(Int, 2*half_langevin_k), "_sd_trans_", trunc(Int, trans_std_dev*100), "/convex_results.csv"),  test_case_result, ',');
end

perform_experiment()