using Rotations
include("../src/calibration/hand_eye_costs.jl")
include("../src/rotation_sdp_solver.jl")

# Setup a minimal (2 measurement) noiseless base case
n = 2
R = RotZ(π/2)  # Ground truth extrinsic rotation
t = [0; 1; 0]  # Ground truth extrinsic translation
T = [R t; 0 0 0 1]

# Test the MFI paper's sufficient condition for strict convexity
v = [0; 0; 1]
Ta = zeros(2, 4, 4)
Ta[:, 4, 4] .= 1
Ta[1, 1:3, 1:3] = RotX(π/2)
Ta[2, 1:3, 1:3] = RotY(π/2)
Tb = zeros(size(Ta))
# Next set rotations and translations according to v
Tb[1, 1:3, 1:3] = T[1:3, 1:3]'*Ta[1, 1:3, 1:3]*T[1:3, 1:3]
Tb[2, 1:3, 1:3] = T[1:3, 1:3]'*Ta[2, 1:3, 1:3]*T[1:3, 1:3]

Tb[1, 1:3, 4] = R'*(v - Ta[1, 1:3, 1:3]*v) # + rand(3)  # Adding the rand brings back the eigenvalue!
Tb[2, 1:3, 4] = R'*(v - Ta[2, 1:3, 1:3]*v) # + rand(3)

# Finally set the "a" translations
Ta[1, 1:3, 4] = R*Tb[1, 1:3, 4] + t - Ta[1, 1:3, 1:3]*t
Ta[2, 1:3, 4] = R*Tb[2, 1:3, 4] + t - Ta[2, 1:3, 1:3]*t

# M_rot, M_full = transformation_hand_eye_cost(Ta, Tb)
M_rot, M_full = transformation_hand_eye_cost_scale(Ta, Tb)

sol, prob = solve_rotation_sdp(M_rot)

# It appears that it can be solved even when not strictly convex! 
# That was a SUFFICIENT but not NECESSARY condition - need to use the full, new Cifuentes conditions!
# Or maybe it is necessary because there are multiple TRANSLATION solutions now? Check!
# OR, maybe co-rank-one Hessian is met but the dual optimal is not zero - check! THIS IS PROBABLY THE PROBLEM!
R_sol = real(extract_rotation(sol))
# t_sol = translation_from_rotation(R_sol, M_full)
t_sol, scale_sol = translation_and_scale_from_rotation(R_sol, M_full)
