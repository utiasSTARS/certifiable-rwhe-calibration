using Rotations
using LinearAlgebra
# using Distributions
using SpecialFunctions

function langevin_sample(κ)
    axis = random_direction()
    angle = von_mises(0, 2*κ)
    return AngleAxis(angle, axis[1], axis[2], axis[3])
end

function random_rotation_sample_normal_magnitude(mean_angle, std_dev)
    axis = random_direction()
    angle = mean_angle + randn()*std_dev
    return AngleAxis(angle, axis[1], axis[2], axis[3])
end

# Source: http://corysimon.github.io/articles/uniformdistn-on-sphere/
function random_direction()
    theta = 2 * π * rand()
    phi = acos(1 - 2 * rand())
    x = sin(phi) * cos(theta)
    y = sin(phi) * sin(theta)
    z = cos(phi)
    return [x; y; z]
end

# Source: https://dlwhittenbury.github.io/ds-1-sampling-and-visualising-the-von-mises-distribution.html
function von_mises(μ, κ)
    a = 1.0 + sqrt(1.0 + 4.0 * κ^2)
    b = (a - sqrt(2.0 * a)) / (2.0 * κ)
    r = (1.0 + b^2) / (2.0 * b)
    
    while true
        U1 = rand()
        U2 = rand()
        U3 = rand()
        z = cos(π * U1)
        f = (1.0 + r *z) / (r + z)
        c = κ * (r - f)

        if ((c * (2.0 - c) - U2) > 0.0) || (log(c/U2) + 1.0 - c) > 0.0 
            return mod(sign(U3 - 0.5) * acos(f) + μ, 2*π)     
        end
    end
end

# Uses Chi distribution in 3 dimensions (also called Maxwell-Boltzmann distribution).
# Also assumes that local Lie-algebra perturbation is in Euler angles.
function convert_langevins_concentration_to_lie_variance(κ)
    return π/(2*κ*(3*π - 8))
end

# For convenience (since the conversion is the same)
function convert_lie_variance_to_langevins_concentration(var)
    return convert_langevins_concentration_to_lie_variance(var)
end

R1 = langevin_sample(0.2)
R2 = langevin_sample(50)

# function langevin_normalization_constant(κ, dim)
#     if dim == 3
#         return exp(κ)*(besselj(0, 2*κ) - besselj(1, 2*κ))
#     elseif dim == 2
#         return besselj(0, 2*κ)
#     else
#         return Base._throw_argerror("Only accepts dim ∈ {2, 3}.")
#     end
# end