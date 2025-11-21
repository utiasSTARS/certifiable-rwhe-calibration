using Plots
gr()
Plots.GRBackend()

include("../src/utils/experiment_setup.jl")
include("../src/local_solver/local_lie_group_hand_eye_solver.jl")
include("../src/local_solver/linear_initializations.jl")


# Plot cost function vs. step size
# Also check that Lie algebra cost function is close enough to QCQP's

n_motions = 10
# These parameters describe the distribution of the platform's groundtruth motion
mean_rot = 40*π/180
rot_std_dev = 10*π/180
mean_trans = 0.25
trans_std_dev = 0.0625

# These parameters describe the distribution of zero-mean measurement errors
vr = (0.001*mean_rot)^2  # Lie-algebraic (gets converted to Langevins)
vt = (0.001*mean_trans)^2
κ = convert_lie_variance_to_langevins_concentration(vr)


function investigate_SE3_linearized_calibration(A::Vector, B::Vector, X0::Pose3D, var_t=1, var_r=1, max_iters=5, abstol=1e-6, β=0.5, τ=0.5)
    X = X0
    iter_count = 0
    delta_sizes = []
    for _ in 1:max_iters
        
        # Previous iteration's solution
        X_iter = copy(X)
        # Solve for the SE(3) perturbation
        dz, G, e0 = iterate_SE3_linearized_calibration(A, B, X, var_t, var_r)
        # Perform line search for stepsize with Armijo's rule
        # X, γ = armijo_line_search(A, B, X_iter, dz, G, e0, τ, β, var_t, var_r)
        X, γ = discretized_line_search(A, B, X_iter, dz, var_t, var_r, 20)
        iter_count += 1
        # Plot line search
        step_increments = range(0, 1, step=0.01)
        cost_vs_step = []
        for step in step_increments
            X_stepped = apply_perturbation(X_iter, se3(step*dz .+ eps()))
            push!(cost_vs_step, hand_eye_cost(A, B, X_stepped, var_t, var_r))
        end
        p_iter = Plots.plot(step_increments, cost_vs_step, linewidth=3, shape=:star5, framestyle=:box, color=:blue)
        title!("Cost vs. Step Size for Iteration "*string(iter_count))
        xlabel!("Step Size Portion")
        ylabel!("Cost Function")
        scatter!([γ], [hand_eye_cost(A, B, X, var_t, var_r)], markersize=8, color=:black, legend=false)
        display(p_iter)
        
        # Check convergence criterion
        push!(delta_sizes, norm(dz))
        if norm(dz) <= abstol
            break
        end
    end
    return X, iter_count, delta_sizes
end

T_gt = [rand(RotMatrix{3}) 2*rand(3).-1; 0 0 0 1]
A, B = random_test_instance(κ, vt, mean_rot, rot_std_dev, mean_trans, trans_std_dev, T_gt, n_motions)

R_init, t_init = linear_hand_eye_initialization(A, B, vt*ones(n_motions), vr*ones(n_motions))

R_init_perturbed = R_init*random_rotation_sample_normal_magnitude(0., 0.5)
t_init_perturbed = t_init + randn(3)*0.3

T_init = [R_init t_init; 0 0 0 1]
T_init_perturbed = [R_init_perturbed t_init_perturbed; 0 0 0 1]
T_lin, iter_count, delta_sizes = investigate_SE3_linearized_calibration(A, B, Pose3D(T_init_perturbed), vt, vr, 5, 1e-6)

# Check errors:
rot_err_init = abs(rotation_angle(RotMatrix{3}(T_init[1:3, 1:3]'*T_gt[1:3, 1:3])))
trans_err_init = norm(T_init[1:3, 4] - T_gt[1:3, 4])
rot_err_local = abs(rotation_angle(RotMatrix{3}(T_lin[1:3, 1:3]'*T_gt[1:3, 1:3])))
trans_err_local = norm(T_lin[1:3, 4] - T_gt[1:3, 4])
