using Rotations
using Convex, MosekTools, COSMO
using SparseArrays
include("rotation_sdp_constraints.jl")
include("sparse_rotation_sdp_dual_constraints.jl")

# TODO: add the dual solver, compare speed with COSMO and others; also add dual solution extractor from primal (and vice-versa) 
# TODO: refactor this to be a bit more modular (experimental needs will guide)

function solve_equality_constrained_sdp(cost, A, b, solver=COSMO.Optimizer)
    n = size(cost)[1]
    Z = Semidefinite(n)
    constraints = [tr(A[i, :, :]'*Z) == b[i] for i in 1:size(A)[1]]
    prob = minimize(tr(cost'*Z), constraints)
    solve!(prob, solver)
    return evaluate(Z), prob
end

function initialize_constrained_sdp(cost, A_eq, b_eq, A_ineq, b_ineq)
    n = size(cost)[1]
    C = Variable(n, n)
    Z = Semidefinite(n)
    constraints = [tr(A_eq[i, :, :]'*Z) == b_eq[i] for i in 1:size(A_eq)[1]]
    constraints += [tr(A_ineq[i, :, :]'*Z) <= b_ineq[i] for i in 1:size(A_ineq)[1]]
    prob = minimize(tr(cost'*Z), constraints)

    return prob, C, Z
end

function solve_problem(prob, Z, C_var, cost, solver=MOSEK.Optimizer)
    fix!(C_var, cost)
    solve!(prob, solver)
    return evaluate(Z)
end

function solve_constrained_sdp(cost, A_eq, b_eq, A_ineq, b_ineq, solver=COSMO.Optimizer)
    n = size(cost)[1]
    Z = Semidefinite(n)
    constraints = [tr(A_eq[i, :, :]'*Z) == b_eq[i] for i in 1:size(A_eq)[1]]
    constraints += [tr(A_ineq[i, :, :]'*Z) <= b_ineq[i] for i in 1:size(A_ineq)[1]]
    prob = minimize(tr(cost'*Z), constraints)
    solve!(prob, solver)
    return evaluate(Z), prob
end

function solve_rotation_sdp(cost::Matrix, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    dim = Int(sqrt(size(cost, 1) - 1))
    A, b = rotation_constraints(usecolumnorthogonality, usehandedness, dim)
    return solve_equality_constrained_sdp(cost, A, b, solver)
end

function solve_SE3_pose_sdp(cost::Matrix, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    dim = 3
    A, b = padded_rotation_constraints(dim, usecolumnorthogonality, usehandedness, dim)
    return solve_equality_constrained_sdp(cost, A, b, solver)
end

function extract_rotation(Z::Matrix)
    r = eigvecs(Z)[:, end]
    r = r[1:end-1]/r[end]  # Assumes that there's a homogenizing scalar variable 
    return reshape(r, 3, 3)
end

# Function for solving a rotation SDP with an additional inequality-constrained scale variable


# Functions for SDPs over two rotations (i.e., over SO(3)×SO(3))
function solve_double_rotation_sdp(cost::Matrix, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    dim = Int(sqrt((size(cost, 1) - 1)/2))
    A, b = get_double_rotation_constraints(usecolumnorthogonality, usehandedness, dim)
    return solve_equality_constrained_sdp(cost, A, b, solver)
end

function extract_double_rotation(Z::Matrix)
    R1 = extract_rotation(Z[[1:9; 19], [1:9; 19]])
    R2 = extract_rotation(Z[10:19, 10:19])
    return R1, R2
end

function solve_sdp_dual_base(Q::AbstractArray, dict_of_constraints::Dict, solver=COSMO.Optimizer)
    ν = Variable(length(keys(dict_of_constraints)))
    Z = Q
    for key in keys(dict_of_constraints)
        Z += ν[key] * sparse(dict_of_constraints[key][:, 1], dict_of_constraints[key][:, 2], dict_of_constraints[key][:, 3])
    end
    display("Constraints Generated")
    constraints = [Z in :SDP]
    prob = maximize(-ν[end], constraints) #note the final constraint in the dict must be the homogeneous constraint
    result = solve!(prob, solver) #really slow before actually solving with COSMO idk why
    return evaluate(Z), prob
end

#Functions for solving the dual of the robot world problem
function solve_sdp_dual(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    dict_of_constraints = get_n_sparse_rotation_constraints(n, false, usecolumnorthogonality, usehandedness, 3)
    last_row_col = 3 * n + 9 * n + 1
    dict_of_constraints = add_sparse_homogeneous_constraint(last_row_col, dict_of_constraints)
    return solve_sdp_dual_base(Q, dict_of_constraints, solver)
end

#Functions for solving the dual of the robot world problem
function solve_sdp_dual_scale(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    dict_of_constraints = get_n_sparse_rotation_constraints(n, true, usecolumnorthogonality, usehandedness, 3)
    last_row_col = 3 * n + 1 + 9 * n + 1
    dict_of_constraints = add_sparse_homogeneous_constraint(last_row_col, dict_of_constraints)
    return solve_sdp_dual_base(Q, dict_of_constraints, solver)
end

function extract_solution_from_dual(Z::SparseMatrixCSC)
    data = eigen(Matrix(Z))
    solution_index = .≈(data.values, 0, atol=1e-5)
    solution = data.vectors[:, solution_index]
    solution ./= solution[end]
    return solution
end

function solve_sdp_dual_using_concat(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    P, ν = generate_n_rotation_dual_constraints_using_concat_at_i(n, false, usecolumnorthogonality, usehandedness)
    Z = Q + P
    display("Constraints Generated")
    constraints = [Z in :SDP]
    prob = maximize(ν, constraints)
    result = solve!(prob, solver)
    return evaluate(Z), prob
end