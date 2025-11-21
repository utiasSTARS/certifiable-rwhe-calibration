using Test
using Rotations

include("../src/rotation_sdp_constraints.jl")

# See https://docs.julialang.org/en/v1/stdlib/Test/#Unit-Testing for help

function test_constraints(R1::RotMatrix3, R2::RotMatrix3)
    A, b = get_double_rotation_constraints()
    r = [vec(R1); vec(R2); 1]
    [r'*A[i,:,:]*r - b[i] for i in 1:size(A)[1]]
end

@test all([≈(val, 0, atol=1e-9) for val in test_constraints(rand(RotMatrix{3}), rand(RotMatrix{3}))])