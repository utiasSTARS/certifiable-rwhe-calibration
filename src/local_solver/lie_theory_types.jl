using StaticArrays
using Rotations
using LinearAlgebra

# Define an abstract type for our Lie algebras' vector spaces (N=3 or 6 for SO3 and SE3)
# TODO: consider SVector and SMatrix instead for better performance (removes mutability)
const LieVector{N, T} = MVector{N, T}
const se3{T} = LieVector{6, T}
const so3{T} = LieVector{3, T}

# Define transformation matrices
# abstract type Transformation{N, T} <: StaticMatrix{N,N,T} end
const Pose{N, T} = MMatrix{N, N, T}
const Pose3D{T} = Pose{4, T}

# Basic operations
function inv_SE3(H::Pose3D)
    return Pose3D([H[1:3, 1:3]'    -H[1:3, 1:3]'*H[1:3, 4];
                 zeros(1, 3)                          1])
end


# H = [rand(RotMatrix{3}) rand(3); 0 0 0 1]
# H = Pose3D(H)
# inv(H)*H