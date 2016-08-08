# we generate code in this module, so precompile where possible
#__precompile__(true)

module Rotations

using Compat
using StaticArrays

import Base: convert, eltype, size, length, getindex, inv, *, Tuple, eye

include("core_types.jl")
include("quaternion_types.jl")
include("angleaxis_types.jl")
include("euler_types.jl")
include("mean.jl")
#include("rot_covariances.jl")
include("derivatives.jl")


export
    Rotation, RotMatrix,
    Quat, SPQuat,
    AngleAxis, RodriguesVec,
    RotX, RotY, RotZ,
    RotXY, RotYX, RotZX, RotXZ, RotYZ, RotZY,
    RotXYX, RotYXY, RotZXZ, RotXZX, RotYZY, RotZYZ,
    RotXYZ, RotYXZ, RotZXY, RotXZY, RotYZX, RotZYX,

    # check validity of the rotation (is it close to unitary?)
    isrotation,

    # angle and axis introspection
    rotation_angle,
    rotation_axis,

    # derivatives (names clash with ForwarDiff?)
    jacobian, hessian

end # module
