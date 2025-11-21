using LinearAlgebra
using Rotations

include("lie_theory_types.jl")

"""
Functions for calculus with SE(3) based on Barfoot 2017, chapter 7 (pg. 252 particularly).
"""

# Compute the adjoint of an SE(3) matrix.
function adjoint(T::Pose3D)
    R =   T[1:3, 1:3]
    return [R            skew(T[1:3, 4])*R;
            zeros(3, 3)                  R]
end

# Compute the adjoint of an se(3) element (this is the curly up-wedge in Barfoot 2017, pg. 227).
function adjoint(ξ::se3)
    ρ = ξ[1:3]
    ϕ = ξ[4:6]
    return [skew(ϕ)     skew(ρ);
            zeros(3, 3) skew(ϕ)]
end

# Compute the left Jacobian of an se(3) element (see Barfoot 2017 pgs. 234, 253).
function jacobian(ξ::se3)
    ξa = adjoint(ξ)
    ϕ = sqrt(sum(ξ[4:6].^2))
    α1 = (4 - ϕ*sin(ϕ) - 4*cos(ϕ))/(2*ϕ^2)
    α2 = (4*ϕ - 5*sin(ϕ) + ϕ*cos(ϕ))/(2*ϕ^3)
    α3 = (2 - ϕ*sin(ϕ) - 2*cos(ϕ))/(2*ϕ^4)
    α4 = (2*ϕ - 3*sin(ϕ) + ϕ*cos(ϕ))/(2*ϕ^5)
    return I(6) + α1*ξa + α2*ξa^2 + α3*ξa^3 + α4*ξa^4
end

function jacobian(ϕ::so3)
    p = norm(ϕ)
    a = ϕ/p
    return sin(p)*I(3)/p + (1 - sin(p)/p)a*a' - p*skew(a)/2
end

function inv_jacobian(ϕ::so3)
    p = norm(ϕ)
    a = ϕ/p
    return (p/2)*cot(p/2)*I(3) + (1 - (p/2)*cot(p/2))a*a' - p*skew(a)/2
end

# Skew-symmetric cross-product matrix central to Lie-algebraic SO(3) 
function skew(ϕ)
    return [0    -ϕ[3] ϕ[2];
             ϕ[3] 0   -ϕ[1];
            -ϕ[2] ϕ[1] 0]
end

function so3_vee(Φ::Matrix)
    return so3([-Φ[2, 3]; Φ[1, 3]; -Φ[1, 2]])
end

# Vee operator brings a matrix representation of se(3) to its vector form
function vee(Ξ::Matrix)
    ρ = Ξ[1:3, 4]
    ϕ = so3_vee(Ξ[1:3, 1:3])
    return se3([ρ; ϕ])
end

function ln(R::RotMatrix{3})
    # R = RotMatrix{3}(R)
    a = rotation_axis(R)
    ϕ = rotation_angle(R) # I believe it uses atan2 and therefore -π < ϕ < π, otherwise see Barfoot pg. 221
    return skew(ϕ*a)
end

# Logarithm for SE(3) matrices
function ln(T::Pose3D)
    Φ = ln(RotMatrix{3}(T[1:3, 1:3]))
    ρ = so3(inv_jacobian(so3_vee(Φ))*T[1:3, 4])
    return [Φ      ρ;
            zeros(1, 4)]
end

function exp(ϕ::so3{<:Real})
    p = norm(ϕ)
    a = ϕ/p
    return cos(p)*I(3) + (1 - cos(p))*a*a' + sin(p)*skew(a)
end

# Exponent for se(3) elements
function exp(ξ::se3)::Pose3D
    return Pose3D([exp(so3(ξ[4:6])) jacobian(so3(ξ[4:6]))*ξ[1:3]; 0 0 0 1])
end

# Apply the full (not approximate) exponent for this perturbation
function apply_perturbation(X::Pose3D, dz::se3)
    return exp(dz)*X
end