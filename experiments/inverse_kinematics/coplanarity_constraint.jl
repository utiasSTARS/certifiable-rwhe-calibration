using LinearAlgebra
using Symbolics

@variables P[1:3, 1:4]

P_homog = [P; 1 1 1 1]

det(P_homog)
