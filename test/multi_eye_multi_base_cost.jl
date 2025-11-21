using MKL
using Test
using Rotations
using SparseArrays

include("../src/calibration/robot_world_costs.jl");

function generate_n_eyes()
    n = rand(1:256)
    #Generate Y and A
    Y = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    A = permutedims([[rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; 
                     [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; 
                     [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]], [3, 1, 2])
    
    #Initialize X and B
    X = zeros(n, 4, 4)
    B = zeros(n, 3, 4, 4)

    #Generate X and B
    X[:, 4, 4] .= 1
    X[:, 1:3, 4] = rand(Float64, n, 3)
    for i=1:n
        X[i, 1:3, 1:3] = rand(RotMatrix{3})
        for j=1:3
            B[i, j, :, :] = inv(Y) * A[j, :, :] * X[i, :, :]
        end
    end
    τ = ones(n, 3)
    κ = ones(n, 3)
    return A, B, X, Y, τ, κ, n
end

function generate_n_bases()
    n = rand(1:256)
    #Generate Y and A
    X = [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]
    B = permutedims([[rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; 
                     [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1];;; 
                     [rand(RotMatrix{3}) rand(Float64, 3); zeros(3)' 1]], [3, 1, 2])
    
    #Initialize X and B
    Y = zeros(n, 4, 4)
    A = zeros(n, 3, 4, 4)

    #Generate X and B
    Y[:, 4, 4] .= 1
    Y[:, 1:3, 4] = rand(Float64, n, 3)
    for i=1:n
        Y[i, 1:3, 1:3] = rand(RotMatrix{3})
        for j=1:3
            A[i, j, :, :] = Y[i, :, :] * B[j, :, :] * inv(X)
        end
    end
    τ = ones(n, 3)
    κ = ones(n, 3)
    return A, B, X, Y, τ, κ, n
end

function construct_n_eyes_1_base_r(X::AbstractArray, Y::AbstractArray, n::Int64)
    r = ones(12*(1+n)+1)
    for i=1:n
        r[3*(i-1)+1:3*i] = X[i, 1:3, 4]
        r[3*(1+n)+9*(i-1)+1:3*(1+n)+9*i] = vec(X[i, 1:3, 1:3])
    end
    r[3*n+1:3*n+3] = Y[1:3, 4]
    r[12*n+3+1:12*n+12] = Y[1:3, 1:3]
    r[end] = 1
    return r
end

function construct_1_eye_n_bases_r(X::AbstractArray, Y::AbstractArray, m::Int64)
    r = ones(12*(1+m)+1)
    r[1:3] = X[1:3, 4]
    r[3*(m+1)+1:3*(m+1)+9] = vec(X[1:3, 1:3])
    for i=1:m
        r[3*i+1:3*(i+1)] = Y[i, 1:3, 4]
        r[3*(m+1)+9*i+1:3*(1+m)+9*(i+1)] = vec(Y[i, 1:3, 1:3])
    end
    r[end] = 1
    return r
end

function test_multi_eye_cost()
    A, B, X, Y, τ, κ, n = generate_n_eyes()
    r = construct_n_eyes_1_base_r(X,Y,n)
    Q = get_empty_sparse_cost_matrix(n, 1, false)
    for i=1:n
        Q_temp = sparse_robot_world_transformation_cost(A, B[i, :, :, :], κ[i, :], τ[i, :])
        Q += get_ith_eye_jth_base_cost(i, 1, n, 1, Q_temp, false)
    end
    return [r'*Q*r] 
end

function test_multi_base_cost()
    A, B, X, Y, τ, κ, n = generate_n_bases()
    r = construct_1_eye_n_bases_r(X,Y,n)
    Q = get_empty_sparse_cost_matrix(1, n, false)
    for i=1:n
        Q_temp = sparse_robot_world_transformation_cost(A[i, :, :, :], B, κ[i, :], τ[i, :])
        Q += get_ith_eye_jth_base_cost(1, i, 1, n, Q_temp, false)
    end
    return [r'*Q*r] 
end

@test all([≈(val, 0, atol=1e-9) for val in test_multi_eye_cost()])
@test all([≈(val, 0, atol=1e-9) for val in test_multi_base_cost()])