# we generate code in this module, so precompile where possible
 VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module Rotations

using FixedSizeArrays
using Quaternions  # TODO: Switch to the UnitQuaternions package? https://github.com/kam3k/UnitQuaternions.jl

import Base.convert, Base.mean, Base.eltype, Base.isnan, Base.length, Base.call, Base.getindex, Base.inv, Base.*
import Quaternions.Quaternion

src_files = ["rotation_types", "rotation_conversions", "euler_conversions", "type_methods", "rotate", "quaternion_ops", "math_funcs", "rot_covariances"]
for file in src_files
    #println("$(file).jl")
    include("$(file).jl")
end


export

        # Euler angle ordering
        EulerOrder,
        EulerZXY,

        # representations
        RotMatrix,
        SpQuat,
        AngleAxis,
        EulerAngles,         # order is in the type specification
        ProperEulerAngles,   # order is in the type specification

        # quaternions as well
        Quaternion,

        # retrieve the euler order
        euler_order

        # check validity of the rotation
        valid_rotation,

        # extra quaternion operations
        rot_angle  # N.B. angle(Quaternions) is different and defined in Quaternions

end # end rotations
