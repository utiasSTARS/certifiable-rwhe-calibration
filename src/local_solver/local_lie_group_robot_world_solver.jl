include("lie_theory_helpers.jl")
"""
Locally solve the multi-sensor AX=YB problem with Barfoot's SE(3) perturbation 
strategy. Compare this custom solver to Optim (and maybe Ipopt in JuMP).

Main reference: pose graph relaxation example (pg. 329) in Barfoot 2017.
"""

using Optim

# Use Gauss-Newton with SE(3)-perturbations as unconstrained variables to solve AX=YB 
function solve_SE3_linearized_calibration(A::Vector, B::Vector, X0::Vector, Y0::Pose3D, max_iters=5)
    X = X0
    Y = Y0
    for iter in 1:max_iters
        # Solve for the SE(3) perturbation
        dz = iterate_SE3_linearized_calibration(A, B, X, Y)
        # Apply the perturbed solution
        X, Y = apply_perturbation(X, Y, dz)
    end
    return X, Y
end

function apply_perturbation(X::Vector, Y::Pose3D, dz)
    X_out = []
    for i in 1:length(X)
        push!(X_out, apply_perturbation(X[i], se3(dz[6*i+1:6*i+6])))
    end
    return X_out, apply_perturbation(Y, se3(dz[1:6]))
end

# Compute a residual using the SE(3) perturbation in Barfoot 2017.
function linearized_residual(A::Pose3D, B::Pose3D, X0::Pose3D, Y0::Pose3D)
    e_0 = vee(ln(A*X0*inv_SE3(B)*inv_SE3(Y0)))

    TY = adjoint(Y0)
    TB = adjoint(B)
    TX0 = adjoint(X0)

    G = [-inv(jacobian(-e_0))*TY*TB*inv(TX0) inv(jacobian(-e_0))]

    return e_0, G
end

# Construct the SE(3) linearized cost matrix and vector for a multi-sensor AXᵢ=YB problem. TODO: add covariances
function iterate_SE3_linearized_calibration(A::Vector, B::Vector, X0::Vector, Y0::Pose3D)
    n_sensors = length(A)
    N = 6*(n_sensors + 1)
    M = zeros(N, N)
    b = zeros(N)
    for i in 1:n_sensors
        Ai = A[i]
        Bi = B[i]
        m = size(Ai, 1)
        G = zeros(6, 12)
        e0 = zeros(6)
        X0i = X0[i]
        for j in 1:m
            e0j, Gj = linearized_residual(Pose3D(Ai[j, :, :]), Pose3D(Bi[j, :, :]), X0i, Y0)
            G += Gj
            e0 += e0j
        end
        # Place the sub-blocks in the batch variables A, b
        Pi = variable_projection(i+1, n_sensors)
        M += Pi'*G'*G*Pi
        b += Pi'*G'*e0
    end
    # Solve the resulting linear system M*dz = b
    # println(M)
    # println(b)
    return M\b
end

function variable_projection(i::Int, n_sensors::Int)
    p = zeros(2, n_sensors+1)
    p[2, 1] = 1
    p[1, i] = 1
    return kron(p, I(6))
end


# Simple test - can't use Identity (breaks the inverse Jacobian)
A = Pose3D([rand(RotMatrix{3}) rand(3); 0 0 0 1]) 
B = Pose3D([rand(RotMatrix{3}) rand(3); 0 0 0 1])
X0 = Pose3D([rand(RotMatrix{3}) rand(3); 0 0 0 1])
Y0 = Pose3D([rand(RotMatrix{3}) rand(3); 0 0 0 1])

e0, G = linearized_residual(A, B, X0, Y0)
