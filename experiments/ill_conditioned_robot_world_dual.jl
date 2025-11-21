using NPZ

include("../src/calibration/robot_world_costs.jl");
include("../src/rotation_sdp_solver.jl");

A = npzread("data/robot_world_ill_conditioned/A.npy");
B = npzread("data/robot_world_ill_conditioned/B.npy");
X = npzread("data/robot_world_ill_conditioned/X.npy");
Y = npzread("data/robot_world_ill_conditioned/Y.npy");

#cost, Q_full = transformation_robot_world_cost(A, B);
Q_full = sparse_robot_world_transformation_cost(A, B, ones(size(A, 1)), ones(size(A, 1)))
MOI = Convex.MOI
solver = MOI.OptimizerWithAttributes(COSMO.Optimizer, "eps_abs" => 5e-11, "eps_rel" =>5e-11)
#Z, prob = solve_sdp_dual(Q_full, 2, true, true, solver);
Z, prob = solve_sdp_dual_using_concat(Q_full, 2, true, true, solver)
solution = extract_solution_from_dual(Z);
display(solution[1:3] - X[1:3, 4])
display(solution[4:6] - Y[1:3, 4])
display(reshape(solution[7:15], 3, 3) * X[1:3, 1:3]')
display(reshape(solution[16:24], 3, 3) * Y[1:3, 1:3]')

