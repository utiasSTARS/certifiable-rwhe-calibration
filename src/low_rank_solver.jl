using JuMP  # TODO: may need to change to import because of clash with Convex.jl
import Ipopt

include("rotation_sdp_constraints.jl")

function low_rank_rotation_qcqp(cost, r=1, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    A, b = rotation_constraints(usecolumnorthogonality, usehandedness, dim)
    return low_rank_qcqp(cost, A, b, r)
end

function low_rank_double_rotation_qcqp(cost, r=1, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    A, b = get_double_rotation_constraints(usecolumnorthogonality, usehandedness, dim)
    return low_rank_qcqp(cost, A, b, r)
end

# TODO: consider sparse matrices!
# TODO: do we need redundant constraints for the local QCQP solver? I think so
function low_rank_qcqp(cost, A, b, r=1)

    n = size(cost, 1)

    model = Model(Ipopt.Optimizer)
    # set_silent(model)
    # Y_init = zeros(n, r)
    # Y_init[1, 1] = 1
    # Y_init[5, 1] = 1
    # Y_init[9, 1] = 1
    @variable(model, Y[1:n, 1:r])
    @objective(model, Min, tr((Y*Y')*cost))
    # set_start_value(Y[1, 1], 1)
    # set_start_value(Y[5, 1], 1)
    # set_start_value(Y[9, 1], 1)
    for i in 1:n
        for j in 1:r
            set_start_value(Y[i, j], rand()*2.0 - 1.0)
        end
    end
    for i in 1:size(A, 1)
        @constraint(model, tr((Y*Y')*A[i, :, :]) == b[i])
    end

    optimize!(model)
    
    return value.(Y), model
end

# TODO: add Y_init
function low_rank_qcqp(cost, A_eq, b_eq, A_ineq, b_ineq, r=1)

    n = size(cost, 1)

    model = Model(Ipopt.Optimizer)
    # set_silent(model)
    # Y_init = zeros(n, r)
    @variable(model, Y[1:n, 1:r])
    @objective(model, Min, tr((Y*Y')*cost))
    for i in 1:n
        for j in 1:r
            set_start_value(Y[i, j], rand())
        end
    end
    for i in 1:size(A_eq, 1)
        @constraint(model, tr((Y*Y')*A_eq[i, :, :]) == b_eq[i])
    end

    for i in 1:size(A_ineq, 1)
        @constraint(model, tr((Y*Y')*A_ineq[i, :, :]) <= b_ineq[i])
    end

    optimize!(model)
    
    return value.(Y), model
end