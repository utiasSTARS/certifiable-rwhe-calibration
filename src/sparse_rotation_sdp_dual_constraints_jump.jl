using LinearAlgebra
using JuMP
using SparseArrays
using BlockDiagonals

function symmetric_variable_matrix(x::Vector{VariableRef})
    X = [x[1] 0.5*x[2] 0.5*x[3];
        0.5*x[2] x[4] 0.5*x[5];
        0.5*x[3] 0.5*x[5] x[6]]
    return X
end

function skew_variable_matrix(x::Vector{VariableRef})
    X = sparse([0 -0.5 0; 0.5 0 0; 0 0 0])*x[3] + sparse([0 0 0.5; 0 0 0; -0.5 0 0])*x[2] + sparse([0 0 0; 0 0 -0.5; 0 0.5 0])*x[1]
    return X
end

function pad_P(P::Matrix{AffExpr}, n::Int64, scale::Bool)
    if scale
        P = vcat(hcat(spzeros(3*n + 1, 3*n + 1), spzeros(3*n + 1, 9*n+1)),
            hcat(spzeros(9*n+1, 3*n + 1), P))
    else
        P = vcat(hcat(spzeros(3*n, 3*n), spzeros(3*n, 9*n+1)),
            hcat(spzeros(9*n+1, 3*n), P))
    end
    return P
end

function generate_n_rotation_dual_constraints_using_block_diag(n::Int64, scale::Bool, model, usecolumnorthogonality::Bool=true, usehandedness::Bool=true)
    @variable(model, ν_r[1:6*n])
    
    if usecolumnorthogonality & usehandedness
        @variable(model, ν_c[1:6*n])
        @variable(model, ν_h[1:9*n])
        #Generate base P
        P = BlockDiagonal([generate_base_ortho_rows(ν_r[6*(i-1)+1:6*i]) + generate_base_ortho_cols(ν_c[6*(i-1)+1:6*i]) + generate_base_handedness_no_vec(ν_h[9*(i-1)+1:9*i]) for i=1:n])
        
        #Add homogeneous entries
        homog_entry = repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_r + repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_c
        homog_vec = -0.5*ν_h[BlockDiagonal([kron([0 1 0; 0 0 1; 1 0 0], Matrix(I, 3, 3)) for i=1:n])*(1:9*n|>collect)]
        P = vcat(hcat(P, homog_vec), hcat(homog_vec', homog_entry))
        
        #Add homogeneous constraint
        @variable(model, ν_homog)
        @objective(model, Max, ν_homog);
        P = P - vcat(hcat(spzeros(9*n, 9*n), spzeros(9*n, 1)), hcat(spzeros(1, 9*n), ν_homog))
        
        P = pad_P(P, n, scale)
        return P
    end
    if usecolumnorthogonality
        @variable(model, ν_c[1:6*n])
        #Generate base P
        P = BlockDiagonal([generate_base_ortho_rows(ν_r[6*(i-1)+1:6*i]) + generate_base_ortho_cols(ν_c[6*(i-1)+1:6*i]) for i=1:n])
        
        #Add homogeneous entries
        homog_entry = repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_r + repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_c
        P = vcat(hcat(P, zeros(size(P,1), 1)), hcat(zeros(1, size(P,1)), homog_entry))
        
        #Add homogeneous constraint
        @variable(model, ν_homog)
        @objective(model, Max, ν_homog);
        P = P - vcat(hcat(spzeros(9*n, 9*n), spzeros(9*n, 1)), hcat(spzeros(1, 9*n), ν_homog))
        
        P = pad_P(P, n, scale)
        return P
    end
    if usehandedness
        @variable(model, ν_h[1:9*n])

        #Generate base P
        P = BlockDiagonal([generate_base_ortho_rows(ν_r[6*(i-1)+1:6*i]) + generate_base_handedness_no_vec(ν_h[9*(i-1)+1:9*i]) for i=1:n])
        
        #Add homogeneous entries
        homog_entry = repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_r
        homog_vec = -0.5*ν_h[BlockDiagonal([kron([0 1 0; 0 0 1; 1 0 0], Matrix(I, 3, 3)) for i=1:n])*(1:9*n|>collect)]
        P = vcat(hcat(P, homog_vec), hcat(homog_vec', homog_entry))
        
        #Add homogeneous constraint
        @variable(model, ν_homog)
        @objective(model, Max, ν_homog);
        P = P - vcat(hcat(spzeros(9*n, 9*n), spzeros(9*n, 1)), hcat(spzeros(1, 9*n), ν_homog))
        
        P = pad_P(P, n, scale)
        return P
    end

    #Generate base P
    P = BlockDiagonal([generate_base_ortho_rows(ν_r[6*(i-1)+1:6*i]) for i=1:n])
        
    #Add homogeneous entries
    homog_entry = repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_r
    P = vcat(hcat(P, zeros(size(P,1), 1)), hcat(zeros(1, size(P,1)), homog_entry))
        
    #Add homogeneous constraint
    @variable(model, ν_homog)
    @objective(model, Max, ν_homog);
    P = P - vcat(hcat(spzeros(9*n, 9*n), spzeros(9*n, 1)), hcat(spzeros(1, 9*n), ν_homog))
        
    P = pad_P(P, n, scale)
    return P
end


function generate_n_rotation_dual_constraints_using_concat_at_i(n::Int64, scale::Bool, model, usecolumnorthogonality::Bool=true, usehandedness::Bool=true)
    @variable(model, ν_r[1:6*n])
    
    #Create base P
    P = generate_base_ortho_rows(ν_r[1:6])
    if usecolumnorthogonality
        @variable(model, ν_c[1:6*n])
        P += generate_base_ortho_cols(ν_c[1:6]) 
    end 
    if usehandedness
        @variable(model, ν_h[1:9*n])
        P_temp, homog_vec = generate_base_handedness(ν_h[1:9])
        P = P + P_temp
    end
    
    #Concatenate additional rotations
    for i=2:n
        P_sum_temp = generate_base_ortho_rows(ν_r[6*(i-1)+1:6*i])
        if usecolumnorthogonality
        P_sum_temp += generate_base_ortho_cols(ν_c[6*(i-1)+1:6*i])
        end
        if usehandedness
        P_temp, h_temp = generate_base_handedness(ν_h[9*(i-1)+1:9*i])
        P_sum_temp += P_temp
        homog_vec = vcat(homog_vec, h_temp)
        end
        P = vcat(hcat(P, spzeros(9*(i-1), 9)), hcat(spzeros(9, 9*(i-1)), P_sum_temp))
    end

    #Generate homogenous component
    homog_entry = repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_r + repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν_c
    P = vcat(hcat(P, homog_vec), hcat(homog_vec', homog_entry))

    #Add homogeneous constraint
    @variable(model, ν_homog)
    @objective(model, Max, ν_homog);
    P = P - vcat(hcat(spzeros(9*n, 9*n), spzeros(9*n, 1)), hcat(spzeros(1, 9*n), ν_homog))

    if scale
        P = vcat(hcat(spzeros(3*n + 1, 3*n + 1), spzeros(3*n + 1, 9*n+1)),
            hcat(spzeros(9*n+1, 3*n + 1), P))
    else
        P = vcat(hcat(spzeros(3*n, 3*n), spzeros(3*n, 9*n+1)),
            hcat(spzeros(9*n+1, 3*n), P))
    end

    return P

end

function generate_n_rotation_dual_constraints_using_concat(n::Int64, scale::Bool, usecolumnorthogonality::Bool=true, usehandedness::Bool=true)
    P, ν_temp = generate_n_ortho_rows(n)
    
    if usecolumnorthogonality
        P_temp, ν_temp = generate_n_ortho_cols(n)
        P = P + P_temp
    end

    if usehandedness
        P_temp, ν_temp = generate_n_handedness_constraints(n)
        P = P + P_temp
    end

    #Add homogeneous constraint
    ν_homog = Variable(1)
    P = P - vcat(hcat(spzeros(9*n, 9*n), spzeros(9*n, 1)), hcat(spzeros(1, 9*n), ν_homog))
    if scale
        P = vcat(hcat(spzeros(3*n + 1, 3*n + 1), spzeros(3*n + 1, 9*n+1)),
            hcat(spzeros(9*n+1, 3*n + 1), P))
    else
        P = vcat(hcat(spzeros(3*n, 3*n), spzeros(3*n, 9*n+1)),
            hcat(spzeros(9*n+1, 3*n), P))
    end
    return P, ν_homog
end

function generate_n_ortho_rows(n::Int64)
    return generate_n_ortho_constraints(n, generate_base_ortho_rows)
end

function generate_n_ortho_cols(n::Int64)
    return generate_n_ortho_constraints(n, generate_base_ortho_cols)
end

function generate_n_ortho_constraints(n::Int64, f::Function)
    ν =Variable(6 * n)
    
    #Create base P
    P = f(ν[1:6])
    
    #Concatenate additional rotations
    for i=2:n
        P_temp = f(ν[6*(i-1)+1:6*i])
        P = vcat(hcat(P, spzeros(9*(i-1), 9)), hcat(spzeros(9, 9*(i-1)), P_temp))
    end

    #Generate homogenous component
    homog_entry = repeat(sparsevec([1, 4, 6], [1, 1, 1]), n)'*ν
    P = vcat(hcat(P, spzeros(9*n)), hcat(spzeros(1, 9*n), homog_entry))

    return P, ν
end

function generate_n_handedness_constraints(n::Int64)
    ν = Variable(9 * n)

    #Create base P and homog_vec
    P, homog_vec = generate_base_handedness(ν[1:9])

    #Generate the values
    for i=2:n
        P_temp, h_temp = generate_base_handedness(ν[9*(i-1)+1:9*i])
        P = vcat(hcat(P, spzeros(9*(i-1), 9)), hcat(spzeros(9, 9*(i-1)), P_temp))
        homog_vec = vcat(homog_vec, h_temp)
    end

    #Place homogenous components
    P = vcat(hcat(P, homog_vec), hcat(homog_vec', 0))
    return P, ν
end

function generate_base_ortho_rows(x::Vector{VariableRef})
    X = symmetric_variable_matrix(x)
    return kron(sparse(Matrix(-I, 3, 3)), X)
end

function generate_base_ortho_cols(x::Vector{VariableRef})
    X = symmetric_variable_matrix(x)
    return kron(X, sparse(Matrix(-I, 3, 3)))
end

function generate_base_handedness(x::Vector{VariableRef})
    h_123 = x[1:3]
    h_231 = x[4:6]
    h_312 = x[7:9]
    P_base = vcat(hcat(spzeros(3, 3), -skew_variable_matrix(h_123), skew_variable_matrix(h_312)),
                    hcat(skew_variable_matrix(h_123), spzeros(3, 3), -skew_variable_matrix(h_231)),
                    -hcat(skew_variable_matrix(h_312), skew_variable_matrix(h_231), spzeros(3, 3)))
    homog_base = vcat(-0.5*h_231, -0.5*h_312, -0.5*h_123)
    return P_base, homog_base
end

function generate_base_handedness_no_vec(x::Vector{VariableRef})
    h_123 = x[1:3]
    h_231 = x[4:6]
    h_312 = x[7:9]
    P_base = vcat(hcat(spzeros(3, 3), -skew_variable_matrix(h_123), skew_variable_matrix(h_312)),
                    hcat(skew_variable_matrix(h_123), spzeros(3, 3), -skew_variable_matrix(h_231)),
                    -hcat(skew_variable_matrix(h_312), skew_variable_matrix(h_231), spzeros(3, 3)))
    return P_base
end