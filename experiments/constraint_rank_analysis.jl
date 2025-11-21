include("../src/rotation_sdp_constraints.jl")

function sym_to_vec(A)::Vector
    n = size(A, 1)
    v = zeros(Int(n*(n+1)/2))
    k = 0
    for i in 1:n
        for j in i:n
            k += 1
            v[k] = A[i, j]
        end
    end
    return v
end


A_row, _ = roworthogonality()
n = size(A_row, 2)
N = Int(n*(n+1)/2)
V = zeros(N, 0)
for i in 1:size(A_row, 1)
    global V = hcat(V, sym_to_vec(A_row[i, :, :]))
end
A_col, _ = columnorthogonality()
for i in 1:size(A_col, 1)
    global V = hcat(V, sym_to_vec(A_col[i, :, :]))
end
A_h, _ = handedness()
for i in 1:size(A_h, 1)
    global V = hcat(V, sym_to_vec(A_h[i, :, :]))
end

println(rank(V))