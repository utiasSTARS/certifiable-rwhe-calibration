using Test
using Rotations

include("../src/rotation_sdp_constraints.jl")

# See https://docs.julialang.org/en/v1/stdlib/Test/#Unit-Testing for help

function test_row_orthogonality(R::RotMatrix3)
    A_row, b_row = row_orthogonality()
    r = [vec(R); 1]
    [r'*A_row[i,:,:]*r - b_row[i] for i in 1:size(A_row)[1]]
end

function test_column_orthogonality(R::RotMatrix3)
    A_col, b_col = column_orthogonality()
    r = [vec(R); 1]
    [r'*A_col[i,:,:]*r - b_col[i] for i in 1:size(A_col)[1]]
end

function test_handedness(R::RotMatrix3)
    A_hand, b_hand = handedness()
    r = [vec(R); 1]
    [r'*A_hand[i,:,:]*r - b_hand[i] for i in 1:size(A_hand)[1]]
end


@test all([≈(val, 0, atol=1e-9) for val in test_row_orthogonality(rand(RotMatrix{3}))])
@test all([≈(val, 0, atol=1e-9) for val in test_column_orthogonality(rand(RotMatrix{3}))])
@test all([≈(val, 0, atol=1e-9) for val in test_handedness(rand(RotMatrix{3}))])