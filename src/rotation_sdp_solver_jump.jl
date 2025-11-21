using Rotations
using JuMP, MosekTools, COSMO
using SparseArrays

include("rotation_sdp_constraints.jl")
include("sparse_rotation_sdp_dual_constraints_jump.jl")

function solve_sdp_dual_jump(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true)
    dict_of_constraints = get_n_sparse_rotation_constraints(n, false, usecolumnorthogonality, usehandedness, 3)
    last_row_col = 3 * n + 9 * n + 1
    dict_of_constraints = add_sparse_homogeneous_constraint(last_row_col, dict_of_constraints)
    model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer,"eps_abs" => 5e-11, "eps_rel" =>5e-11));
    @variable(model, x[1:length(keys(dict_of_constraints))])
    Z = Q
    for key in keys(dict_of_constraints)
        Z += x[key] * sparse(dict_of_constraints[key][:, 1], dict_of_constraints[key][:, 2], dict_of_constraints[key][:, 3])
    end
    @constraint(model, Z >= 0, PSDCone());
    display("Constraints Generated")
    @objective(model, Max, -x[end]);
    JuMP.optimize!(model);
    return value.(Z), model
end

function solve_sdp_dual_using_concat(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer,"eps_abs" => 5e-11, "eps_rel" =>5e-11, "max_iter"=>1000000, "rho"=>1e-2));
    #P = generate_n_rotation_dual_constraints_using_concat_at_i(n, false, model, usecolumnorthogonality, usehandedness)
    P = generate_n_rotation_dual_constraints_using_block_diag(n, false, model, usecolumnorthogonality, usehandedness)
    Z = Q + P
    display("Constraints Generated")
    @constraint(model, Z >= 0, PSDCone());
    JuMP.optimize!(model);
    return value.(Z), model
end

function solve_sdp_dual_w_scale_using_concat(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer,"eps_abs" => 5e-11, "eps_rel" =>5e-11));
    #P = generate_n_rotation_dual_constraints_using_concat_at_i(n, false, model, usecolumnorthogonality, usehandedness)
    P = generate_n_rotation_dual_constraints_using_block_diag(n, true, model, usecolumnorthogonality, usehandedness)
    Z = Q + P
    display("Constraints Generated")
    @constraint(model, Z >= 0, PSDCone());
    JuMP.optimize!(model);
    return value.(Z), model
end

function solve_sdp_dual_using_concat_schur(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    #model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer,"eps_abs" => 5e-11, "eps_rel" =>5e-11, "max_iter"=>1000000, "rho"=>1e-4));
    #model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer,"eps_abs" => 1e-8, "eps_rel" =>1e-8, "rho"=>1e-4)); #Tuned for thousand_eyes_restrict
    model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer,"eps_abs" => 5e-11, "eps_rel" =>5e-11, "max_iter"=>1000000, "rho"=>1e-4, "alpha"=>1.0));
    #model = JuMP.Model(Mosek.Optimizer);

    #P = generate_n_rotation_dual_constraints_using_concat_at_i(n, false, model, usecolumnorthogonality, usehandedness)
    #set_silent(model)
    P = generate_n_rotation_dual_constraints_using_block_diag(n, false, model, usecolumnorthogonality, usehandedness)
    Q_schur = Q[3*n+1:end, 3*n+1:end] - Q[3*n+1:end, 1:3*n] * sparse(inv(Matrix(Q[1:3*n,1:3*n]))) * Q[1:3*n, 3*n+1:end]
    Q_schur = (Q_schur + transpose(Q_schur))/2
    Z = Q_schur + P[3*n+1:end, 3*n+1:end];
    #Z = (Z + transpose(Z))/2
    display("Constraints Generated")
    @constraint(model, Z >= 0, PSDCone());
    JuMP.optimize!(model);
    return value.(Z), model
end

function solve_sdp_dual_w_scale_using_schur(Q::AbstractArray, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, solver=COSMO.Optimizer)
    model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer,"eps_abs" => 5e-11, "eps_rel" =>5e-11, "max_iter"=>1000000, "rho"=>1e-4));
    #P = generate_n_rotation_dual_constraints_using_concat_at_i(n, false, model, usecolumnorthogonality, usehandedness)
    P = generate_n_rotation_dual_constraints_using_block_diag(n, true, model, usecolumnorthogonality, usehandedness)
    Q_schur = Q[3*n+2:end, 3*n+2:end] - Q[3*n+2:end, 1:3*n+1] * sparse(inv(Matrix(Q[1:3*n+1,1:3*n+1]))) * Q[1:3*n+1, 3*n+2:end]
    Z = Q_schur + P[3*n+2:end, 3*n+2:end];
    display("Constraints Generated")
    @constraint(model, Z >= 0, PSDCone());
    JuMP.optimize!(model);
    return value.(Z), model
end

function extract_solution_from_dual(Z::SparseMatrixCSC)
    data = eigen(Matrix(Z))
    solution_index = .≈(data.values, 0, atol=1e-4)
    solution = data.vectors[:, solution_index]
    solution ./= solution[end]
    return solution
end

function extract_solution_from_dual_schur(Z::SparseMatrixCSC, Q::SparseMatrixCSC)
    num_params = size(Q, 2)
    data = eigen(Matrix(Z))
    #solution_index = .≈(data.values, 0, atol=1e-5)
    solution_index = argmin(data.values)
    solution = data.vectors[:, solution_index]

    if length(solution) < num_params
        temp = -Q[1:end-length(solution), 1:end-length(solution)] \ (Q[1:end-length(solution), end-length(solution)+1:end]*solution)
        solution = vcat(temp, solution)
    end
    solution ./= solution[end]
    return solution
end