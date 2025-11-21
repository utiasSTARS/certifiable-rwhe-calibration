using LinearAlgebra

# TODO: save these as a big hard-coded file or data for speed for d = 2, 3, 4?
function column_pattern(i::Int, j::Int, dim::Int)
    A = zeros(dim, dim)
    if i == j
        A[i, j] = 1
    else
        A[i, j] = 0.5
        A[j, i] = 0.5
    end
    return kron(A, I(dim))
end

function row_pattern(i::Int, j::Int, dim::Int)
    A = zeros(dim, dim)
    if i == j
        A[i, j] = 1
    else
        A[i, j] = 0.5
        A[j, i] = 0.5
    end
    return kron(I(dim), A)
end

function orthogonality(pattern::Function, dim::Int=3)
    n_constraints = Int(dim*(dim+1)/2)
    A = zeros(n_constraints, dim^2 + 1, dim^2 + 1)
    b = zeros(n_constraints)
    constraint_ind = 1
    for i in 1:dim
        for j in i:dim
            A[constraint_ind, 1:end-1, 1:end-1] = pattern(i, j, dim)
            if i == j
                b[constraint_ind] = 1
            end
            constraint_ind += 1
        end
    end
    return A, b
end

# Computes the homogenous SDP matrix constraints enforcing the lifted R*R' = I  
function row_orthogonality(dim::Int=3)
    return orthogonality(row_pattern, dim)
end

# Computes the homogenous SDP matrix constraints enforcing the lifted R'*R = I 
function column_orthogonality(dim::Int=3)
    return orthogonality(column_pattern, dim)
end

# Returns the simple scalar homogeneity constraint y^2 = 1 (last variable is y)
function homogeneity(dim::Int=3)
    A = zeros(dim^2 + 1, dim^2 + 1)
    A[end, end] = 1
    return A, 1
end

# Computes the cross-product or 'handedness' constraints (only valid for dim=3)
function handedness()
    dim = 3
    ids = [1 2 3; 
           2 3 1;
           3 1 2]
    A = zeros(9, dim^2 + 1, dim^2 + 1)
    b = zeros(9)

    for id in 1:dim
        i, j, k = ids[id, :]
        Eij = zeros(dim, dim)
        Eij[i, j] = 1
        Eij[j, i] = -1
        for l in 1:dim
            Cijl = zeros(dim, dim)
            if l == 1
                Cijl[2, 3] = -1
                Cijl[3, 2] =  1 
            elseif l == 2
                Cijl[1, 3] =  1
                Cijl[3, 1] = -1 
            elseif l == 3
                Cijl[1, 2] = -1
                Cijl[2, 1] =  1 
            end
            A[dim*(id-1) + l, 1:end-1, 1:end-1] = 0.5*kron(Eij, Cijl)
            A[dim*(id-1) + l, end, dim*(k-1) + l] = 0.5
            A[dim*(id-1) + l, dim*(k-1) + l, end] = 0.5
        end
    end

    return A, b
end

function rotation_constraints(usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    A, b = row_orthogonality(dim)
    A_h, b_h = homogeneity(dim)
    A = [A; reshape(A_h, 1, size(A_h, 1), size(A_h, 2))]
    b = [b; [b_h]]

    if usecolumnorthogonality
        A_col, b_col = column_orthogonality(dim)
        A = [A; A_col]
        b = [b; b_col]
    end

    if usehandedness && dim == 3
        A_hand, b_hand = handedness()
        A = [A; A_hand]
        b = [b; b_hand]
    end
    return A, b
end

function padded_rotation_constraints(front_padding::Int=3, back_padding::Int=0, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    A, b = rotation_constraints(usecolumnorthogonality, usehandedness, dim)
    A_pad = front_pad_constraints(front_padding, A)
    if back_padding > 0
        A_pad = back_pad_constraints(back_padding, A_pad)
    end
    return A_pad, b
end

function front_pad_constraints(front_padding::Int, A)
    A_out = zeros(size(A, 1), size(A, 2) + front_padding, size(A, 3) + front_padding)
    for i in 1:size(A, 1)
        A_out[i, front_padding + 1:front_padding+size(A, 2), front_padding + 1:front_padding+size(A, 3)] = A[i, :, :]
    end
    return A_out
end

function back_pad_constraints(back_padding::Int, A)
    A_out = zeros(size(A, 1), size(A, 2) + back_padding, size(A, 3) + back_padding)
    n = size(A, 2) - 1
    for i in 1:size(A, 1)
        A_out[i, 1:n, 1:n] = A[i, 1:n, 1:n]
        A_out[i, end, 1:n] = A[i, end, 1:n]
        A_out[i, 1:n, end] = A[i, 1:n, end]
        A_out[i, end, end] = A[i, end, end]
    end
    return A_out
end

# Constraints for positive scale variable
function positive_scale_constraint(n::Int, index::Int=4)
    A = zeros(n, n)
    A[index, index] = -1
    b = 0
    return A, b
end

function get_double_rotation_constraints(usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    n = 2*dim^2 + 1
    single_constraints, b = rotation_constraints(usecolumnorthogonality, usehandedness, dim)
    # Need to remove the homogeneity constraint
    single_constraints = single_constraints[[1:6; 8:size(single_constraints, 1)], :, :]
    b = b[[1:6; 8:size(b, 1)]]
    m = size(single_constraints, 1)
    constraints = zeros(2*m+1, n, n)
    constraints[1:m, 1:(dim^2), 1:(dim^2)] = single_constraints[:, 1:end-1, 1:end-1]
    constraints[1:m, 1:(dim^2 + 1), end] = single_constraints[:, 1:end, end]
    constraints[1:m, end, 1:(dim^2 + 1)] = single_constraints[:, end, 1:end]
    constraints[m+1:2*m, (dim^2 + 1):end, (dim^2+1):end] = single_constraints
    constraints[end, end, end] = 1 # Homogeneity
    return constraints, [b; b; 1]
end

function get_sparse_dual_rotation_constraints(usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    #Convert the primal rotation constraints into the sparse dual constraints
    
    #Collect the primal constraints
    A, b = rotation_constraints(usecolumnorthogonality, usehandedness, dim) 
    
    #Remove the homogeneity constraint
    A = A[[1:6; 8:size(A, 1)], :, :]
    b = b[[1:6; 8:size(b, 1)]]

    #Make a dict of A matrices
    A_sparse_dict = Dict{Int64, Matrix{Float64}}() 

    #Get the row cols for each entry of A[i, :, :]
    indices = range(1, size(A, 2), step=1) |>collect #Assume that A[i, :, :] is  a square matrix
    row_col = [vec(repeat(indices, 1, length(indices))) vec(repeat(indices', length(indices), 1))]
    
    #Fill the sparse dict
    for i in axes(A, 1)
        A_sparse_temp = [row_col vec(A[i, :, :])]
        A_sparse_temp = A_sparse_temp[A_sparse_temp[:, 3].!=0, :]
        
        #Include b
        A_sparse_temp = [A_sparse_temp; [size(A, 2) size(A, 2) -b[i]]]
        
        #Merge to dictionary
        merge!(A_sparse_dict, Dict(i=>A_sparse_temp))
    end
    return A_sparse_dict
end

function apply_offsets_to_sparse_constraints(i::Int64, base_A_sparse_dict::Dict{Int64, Matrix{Float64}}, base_offset::Int64, last_row_col::Int64)
    #Map the sparse A constraints to the ith of n rotations
    i_of_n_A_sparse_dict = Dict{Int64, Matrix{Float64}}()
    num_constraints = length(keys(base_A_sparse_dict))

    #Loop through base set of keys
    for key in keys(base_A_sparse_dict)
        #Update key and collect entry
        new_key = key + num_constraints * (i - 1)
        temp_entry = base_A_sparse_dict[key]

        #Move homogeneous row col
        temp_entry[temp_entry[:, 1] .== 10, 1] .= last_row_col
        temp_entry[temp_entry[:, 2] .== 10, 2] .= last_row_col

        #Apply base offset to constraints
        temp_entry[temp_entry[:, 1] .< 10, 1] .+= base_offset
        temp_entry[temp_entry[:, 2] .< 10, 2] .+= base_offset

        #Merge into Dict
        merge!(i_of_n_A_sparse_dict, Dict(new_key=>temp_entry))
    end
    return i_of_n_A_sparse_dict
end

function get_i_of_n_sparse_rotation_constraints_w_scale(i::Int64, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    #Create a dict of sparse matrices for the ith rotation matrix of N
    
    #Compute the constants that map the ith of n rotations to the correct location
    base_offset = 3 * n + 1 + 9 * (i - 1) #translations + scale + rotations
    last_row_col = 3 * n + 1 + 9 * n + 1 #translations + scale + rotations + homogeneous
    
    #Get the constraints
    base_A_sparse_dict = get_sparse_dual_rotation_constraints(usecolumnorthogonality, usehandedness, dim)

    i_of_n_A_sparse_dict = apply_offsets_to_sparse_constraints(i, base_A_sparse_dict, base_offset, last_row_col)
    return i_of_n_A_sparse_dict
end

function get_i_of_n_sparse_rotation_constraints(i::Int64, n::Int64, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    #Create a dict of sparse matrices for the ith rotation matrix of N
    
    #Compute the constants that map the ith of n rotations to the correct location
    base_offset = 3 * n + 9 * (i - 1) #translations + rotations
    last_row_col = 3 * n + 9 * n + 1 #translations + rotations + homogeneous
    
    #Get the constraints
    base_A_sparse_dict = get_sparse_dual_rotation_constraints(usecolumnorthogonality, usehandedness, dim)

    i_of_n_A_sparse_dict = apply_offsets_to_sparse_constraints(i, base_A_sparse_dict, base_offset, last_row_col)
    return i_of_n_A_sparse_dict
end

function get_n_sparse_rotation_constraints(n::Int64, scale::Bool, usecolumnorthogonality::Bool=true, usehandedness::Bool=true, dim=3)
    #Create N rotation constraints

    #Initialize the sparse dict
    n_A_sparse_dict = Dict{Int64, Matrix{Float64}}()
    
    #Loop through the n dicts
    for i=1:n
        if scale
            temp_dict = get_i_of_n_sparse_rotation_constraints_w_scale(i, n, usecolumnorthogonality, usehandedness, dim)
        else
            temp_dict = get_i_of_n_sparse_rotation_constraints(i, n, usecolumnorthogonality, usehandedness, dim)
        end
        merge!(n_A_sparse_dict, temp_dict)
    end
    return n_A_sparse_dict
end

function add_sparse_homogeneous_constraint(last_row_col::Int64, A_sparse_dict::Dict{Int64, Matrix{Float64}})
    #Adds a component of the homogeneous constraint to the sparse dictionary for the dual problem
    key = maximum(keys(A_sparse_dict)|>collect) + 1
    entry = [last_row_col last_row_col 1]
    merge!(A_sparse_dict, Dict(key=>entry))
    return A_sparse_dict
end