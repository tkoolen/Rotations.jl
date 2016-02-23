module Rotations

using FixedSizeArrays
using Quaternions

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
	   EulerAngles,  # order is in the type specification

	   # to convert between representations
	   convert_rotation,

	   # check validity of the rotation
	   valid_rotation,

       # extra quaternion operations
       mean, rot_angle,

       # Misc
       projective                              # convert a rotation matrix into a transformation matrix


src_files = ["rotation_types", "euler_types", "fsa_helpers", "type_methods", "rotation_conversions", "euler_conversions", "quaternion_ops"]

for file in src_files 
	println(file)
	include("$(file).jl")
end

end # end rotations
