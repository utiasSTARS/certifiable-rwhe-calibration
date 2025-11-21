using JuMP  # Optimization modelling language
import HiGHS  # Our LP solver

# Even's payoff matrix for Two-Finger Morra
A = [2 -3;  
    -3  4]
model = Model(HiGHS.Optimizer)

# Define variables and constraints
@variable(model, α)
@variable(model, x[1:2] >= 0)
@constraint(model, sum(x) == 1)
@constraint(model, A'*x .>= α)

@objective(model, Max, α)
optimize!(model)
solution_summary(model)
println(value.(x))

