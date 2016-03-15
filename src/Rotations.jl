# we generate code in this module, so precompile where possible
 VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module Rotations

using FixedSizeArrays
using Quaternions  # TODO: Switch to the UnitQuaternions package? https://github.com/kam3k/UnitQuaternions.jl

import Base.convert, Base.mean, Base.eltype, Base.isnan, Base.length, Base.call, Base.getindex
import Quaternions.Quaternion


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

        # to convert between representations
        convert_rotation,
        rotate_point,

        # check validity of the rotation
        valid_rotation,

        # extra quaternion operations
        rot_angle,  # N.B. angle(Quaternions) is different and defined in Quaternions

        # Misc
        projective                              # convert a rotation matrix into a transformation matrix


src_files = ["rotation_types", "fsa_helpers", "type_methods", "rotation_conversions", "euler_conversions", "rotate", "quaternion_ops", "rot_covariances"]
for file in src_files 
    include("$(file).jl")
end

end # end rotations
