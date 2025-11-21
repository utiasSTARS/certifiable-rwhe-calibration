using Rotations
include("../src/calibration/hand_eye_costs.jl")
include("../src/rotation_sdp_solver.jl")

# Setup a minimal (2 measurement) noiseless base case
n = 2
R = RotZ(π/2)  # Ground truth extrinsic rotation
t = [0; 1; 0]  # Ground truth extrinsic translation
T = [R t; 0 0 0 1]
α = 0.5

Ta = zeros(2, 4, 4)
Ta[:, 4, 4] .= 1
Ta[1, 1:3, 1:3] = RotX(π/2)
Ta[1, 1:3, 4]   = [0; 0; 1]
Ta[2, 1:3, 1:3] = RotY(π/2)
Ta[2, 1:3, 4] = [0; 0; 1]

Tb = zeros(size(Ta))
Tb[1, :, :] = inv(T)*Ta[1, :, :]*T
Tb[2, :, :] = inv(T)*Ta[2, :, :]*T

# Apply the scale 
Ta_known_α = copy(Ta)
Ta[:, 1:3, 4] = Ta[:, 1:3, 4]/α


# Check the rank
M_rank_check = [I(3)-Ta[1, 1:3, 1:3] R*Tb[1, 1:3, 4]; 
                I(3)-Ta[2, 1:3, 1:3] R*Tb[2, 1:3, 4]]
# Degenerate Schur inverse for these param (saved for later):
# Ta[1, 1:3, 1:3] = RotX(π/2)
# Ta[1, 1:3, 4]   = [0; 0; 1]
# Ta[2, 1:3, 1:3] = RotY(π/2)
# Ta[2, 1:3, 4] = [0; 0; 1]
# Does this affect the known scale case as well? Not clear yet. Maybe singularities are a problem?

# Next check the inverted element for the Schur complement:
M_tα1 = [I(3)-Ta[1,1:3,1:3] -Ta[1,1:3,4]]
Q_tα  = M_tα1'*M_tα1
M_tα2 = [I(3)-Ta[2,1:3,1:3] -Ta[2,1:3,4]]
Q_tα  += M_tα2'*M_tα2
# This is also related to the matrix in Eq. (42) of Andreff et al.! 
# In that paper, they claim that "two general motions with nonparallel axes" are enough...
# Is the word "general" here key? I understand that this motion is rare and degenerate, but 
# it's not very thorough to leave out these interesting edge cases!
M_andreff = [I(3)-Ta[1,1:3,1:3] -Ta[1,1:3,4];
             I(3)-Ta[2,1:3,1:3] -Ta[2,1:3,4]]

# See if using the nullspace of Q_tα leads to the same cost function:
δ = 1.5
t_perturbed = t + δ*eigvecs(Q_tα)[1:3, 1]
α_perturbed = α + δ*eigvecs(Q_tα)[4, 1]

true_trans_cost = sum([norm2(R*Tb[i, 1:3, 4] + t - Ta[i, 1:3, 1:3]*t - α*Ta[i, 1:3, 4])^2 for i in 1:2])

perturbed_trans_cost = sum([norm2(R*Tb[i, 1:3, 4] + t_perturbed - Ta[i, 1:3, 1:3]*t_perturbed - α_perturbed*Ta[i, 1:3, 4])^2 for i in 1:2])


# Gives inversion error for now
# M_rot, M_full = transformation_hand_eye_cost_scale(Ta, Tb)

# sol, prob = solve_rotation_sdp(M_rot)

