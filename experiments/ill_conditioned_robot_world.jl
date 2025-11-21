using NPZ

include("../src/calibration/robot_world_costs.jl");
include("../src/rotation_sdp_solver.jl");

A = npzread("data/robot_world_ill_conditioned/A.npy");
B = npzread("data/robot_world_ill_conditioned/B.npy");
X = npzread("data/robot_world_ill_conditioned/X.npy");
Y = npzread("data/robot_world_ill_conditioned/Y.npy");

cost, Q_full = transformation_robot_world_cost(A, B);
Z, prob = solve_double_rotation_sdp(cost);
Rx, Ry = extract_double_rotation(Z);
t = translations_from_rotations(Rx, Ry, Q_full);
display(t[1:3] - X[1:3, 4])
display(t[4:6] - X[1:3, 4])
display(Rx * X[1:3, 1:3]')
display(Ry * Y[1:3, 1:3]')

