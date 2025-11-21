include("../src/rotation_sdp_solver.jl")

## Simple test
cost = [1 0; 0 2]
A = zeros(2, 2, 2)
A[1, :, :] = I(2)
A[2, :, :] = [0 0; 0 1]
b = [3, 2]
sol = solve_equality_constrained_sdp(cost, A, b)
##

## Random rotation matrix test
cost = rand(10, 10)
cost = cost + cost'
sol = solve_rotation_sdp(cost)
R = extract_rotation(sol)

@test det(R) ≈ 1 atol=1e-8
@test all(abs.(R*R' - I(3)) .< 1e-6)
##