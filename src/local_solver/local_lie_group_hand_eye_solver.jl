include("lie_theory_helpers.jl")
"""
Locally solve single-sensor AX=XB problem with Barfoot's SE(3) perturbation 
strategy. Compare this custom solver to Optim (and maybe Ipopt in JuMP).

Main reference: pose graph relaxation example (pg. 329) in Barfoot 2017.
"""

using Optim

# Use Gauss-Newton with SE(3)-perturbations as unconstrained variables to solve AX=XB 
function solve_SE3_linearized_calibration(A::Vector, B::Vector, X0::Pose3D, var_t=1, var_r=1, max_iters=10, abstol=1e-6, β=0.2, τ=0.5)
    X = X0
    iter_count = 0
    for _ in 1:max_iters
        X_iter = copy(X)
        # Solve for the SE(3) perturbation
        dz, G, e0 = iterate_SE3_linearized_calibration(A, B, X, var_t, var_r)
        # Perform line search for stepsize with Armijo's rule
        # X, γ = armijo_line_search(A, B, X_iter, dz, G, e0, τ, β, var_t, var_r)
        J_inv = inv_jacobian(so3(e0[4:6]))
        X, γ = discretized_line_search(A, B, X_iter, J_inv, dz, var_t, var_r, 20)
        iter_count += 1
        # Check convergence criterion
        if norm(dz) <= abstol
            break
        end
    end
    return X, iter_count
end

function hand_eye_cost(A::Vector, B::Vector, X::Pose3D, J_inv, var_t=1, var_r=1)
    cost = 0.0
    for i in 1:length(A)
        Ai = Pose3D(A[i])
        Bi = Pose3D(B[i])
        e = vee(ln(Ai*X*inv_SE3(Bi)*inv_SE3(X)))
        cost += (1/var_t) * norm(J_inv*e[1:3]*J_inv')^2 + (1/var_r) * norm(e[4:6])^2
    end
    return cost
end

# Compute a residual using the SE(3) perturbation in Barfoot 2017.
function linearized_residual(A::Pose3D, B::Pose3D, X0::Pose3D)

    # TODO: Talked to Emmett, try inv_SE3(X0)*inv_SE3(A)*X0*B instead! Will change e_0 and G, derive it!
    e_0 = vee(ln(A*X0*inv_SE3(B)*inv_SE3(X0)))

    TB = adjoint(B)
    TX0 = adjoint(X0)

    G = -inv(jacobian(-e_0))*(TX0*TB*inv(TX0) - I)

    return e_0, G
end

# Construct the SE(3) linearized cost matrix and vector for a multi-sensor AX=XB problem. TODO: add covariances
function iterate_SE3_linearized_calibration(A::Vector, B::Vector, X0::Pose3D, var_t, var_r)
    N = 6
    G = zeros(N, N)
    e0 = zeros(6)
    for i in 1:size(A, 1)
        e0_i, G_i = linearized_residual(Pose3D(A[i]), Pose3D(B[i]), X0)
        G += G_i
        e0 += e0_i
    end
    Σ_inv = diagm([(1/var_t)*ones(3); (1/var_r)*ones(3)])

    # Experimental: se(3) is NOT just translation, it is J*ρ, see SE(3) identities in Barfoot 2017. 
    e_J_inv = inv_jacobian(so3(e0[4:6]))
    Σ_inv[1:3, 1:3] = e_J_inv*Σ_inv[1:3, 1:3]*e_J_inv'

    M = G'*Σ_inv*G
    b = G'*Σ_inv*e0
    # println("b (gradient): ")
    # println(b)
    return M\b, G, e0
end

function armijo_line_search(A, B, X_iter, dz, G, e0, τ, β, var_t, var_r)
    J_inv = inv_jacobian(so3(e0[4:6]))
    J0 = hand_eye_cost(A, B, X_iter, J_inv, var_t, var_r)
    X = apply_perturbation(X_iter, se3(dz))
    γ = 1.0
    k = 1
    while hand_eye_cost(A, B, X, J_inv, var_t, var_r) - 2*β*γ*dz'*G'*e0 > J0
        γ = τ*γ
        X = apply_perturbation(X_iter, se3(γ*dz))
        k += 1
    end
    return X, γ
end

function discretized_line_search(A, B, X_iter, J_inv, dz, var_t, var_r, n_grid=10)
    steps = range(0, 1, length=n_grid)
    min_step = NaN
    min_step_cost = Inf
    for step in steps
        X = apply_perturbation(X_iter, se3(step*dz .+ eps()))
        cost = hand_eye_cost(A, B, X, J_inv, var_t, var_r)
        if cost < min_step_cost
            min_step = step
            min_step_cost = cost
        end
    end
    return apply_perturbation(X_iter, se3(min_step*dz .+ eps())), min_step
end
