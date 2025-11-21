include("rotation_sdp_constraints.jl")

# Define constants matrices related to QCQP constraints
const A_row, b_row = roworthogonality()
const A_hand, b_hand = handedness()
const A_homog, b_homog = homogeneity()
const symvec_index_matrix_10 = make_symvec_index_matrix(10)

# TODO: do I need the factors of 2? They cancel out, but better to include them, I think.
# I don't believe the handedness constraints are needed because the gradient is local and SO(3) and O(3)\SO(3) are disconnected
function rotation_qcqp_gradient(cost::Matrix, r::Vector)
    A_dx = solution_gradient(cost, r)
    # return -inv(A_dx)*parameter_gradient(r)
    return -A_dx\parameter_gradient(r)
end

function solution_gradient(cost::Matrix, r::Vector)
    # A_dx = hcat([A_h[i, :, :]*r for i = 1:3]...)
    A_dx = hcat(hcat([A_row[i, :, :]*r for i = 1:size(A_row, 1)]...), A_homog*r)
    return 2.0*[cost A_dx; 
            A_dx' zeros(size(A_dx, 2), size(A_dx, 2))]
end

function parameter_gradient(r::Vector, m::Int = 7)
    n = size(r, 1)
    N = Int((n + 1)*n/2)
    grad = zeros(n + m, N)

    for i = 1:n
        for j = 1:n
            grad[i, symvec_index_matrix_10[i, j]] = r[j]
        end
    end
    return 2.0*grad
end


function make_symvec_index_matrix(n::Int = 10)
    N = Int((n+1)*n/2)
    symvec_index_matrix = zeros(Int, n, n)
    count = 1
    for idx = 1:n 
        for jdx = idx:n 
            symvec_index_matrix[idx, jdx] = count
            symvec_index_matrix[jdx, idx] = count
            count += 1
        end
    end
    return symvec_index_matrix
end

# Test it out - move this to the test folder later
r = zeros(10)
r[1], r[5], r[9], r[10] = ones(4)
cost = rand(10, 10)
cost = 0.5*(cost + cost')
A_z = solution_gradient(cost, r)

println(make_symvec_index_matrix(5))

println(parameter_gradient(collect(1:10)))
