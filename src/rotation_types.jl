

# a default element type
macro DefaultElType(); Float64; end  # need this to be a macro for type stability in type_methods.jl promote_type_sp functions

# some trait style stuff
numel{T}(::Type{T}) = length(fieldnames(T))    # function to get the number of elements in each parameterization
n_params{T}(::Type{T}) = length(T.parameters)  # the number of template parameters

###################################################
# Rotation matrix (typeofalias of Mat{3,3,Float64})
# N.B. all other parameterizations T should implement
# a convert_rotation(T, RotMatrix) method
###################################################

@doc """
A type alias for rotation matrices (alias of Mat{3,3,T})
""" ->
typealias RotMatrix{T <: Real}   Mat{3,3,T}
#immutable RotMatrix{T <: AbstractFloat} <: FixedArray{T,2,Tuple{3,3}} # this would be better but much more effort without an implicit conversion
#    _::Tuple{Tuple{T, T, T}, 
#             Tuple{T, T, T},
#             Tuple{T, T, T}}
#end

 
# fix weirdness, try: mat = Mat(eye(3)); @code_warntype convert(typeof(mat), mat)  
convert{T}(::Type{RotMatrix{T}}, X::RotMatrix{T}) = X

# RotMatrix(eye(3)) hangs, so fix it
# can't use the union types because it wont replace whats in fixed size arrays
Base.convert{U}(::Type{RotMatrix}, mat::Matrix{U}) = RotMatrix{U}(mat)
Base.call{U}(::Type{RotMatrix}, mat::Matrix{U}) = convert(RotMatrix{U}, mat)

# allow converting element types
convert_rotation{T, U}(::Type{RotMatrix{T}}, X::RotMatrix{U}) = convert(RotMatrix{T}, X)

# default parameters
default_params{T <: RotMatrix}(::Type{T}) = (@DefaultElType(), ) 

numel(::Type{RotMatrix}) = 9     # special case
n_params(::Type{RotMatrix})  = 1 # The behaviour of length(RotMatrix.parameters) seems unstable...



###################################################
# Quaternions
# (Quaternions are defined in Quaternions.jl)
###################################################

# add an indexing scheme
# getindex(X::Quaternion, idx::Integer) = getfield(X, idx)

# an extra constructor
call{T}(::Type{Quaternion{T}}, a::Real, b::Real, c::Real, d::Real) = Quaternion(T(a), T(b), T(c), T(d))

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{RotMatrix{T}}, q::Quaternion{T}) = quat_to_rot(q)
convert_rotation{T <: Real}(::Type{Quaternion{T}}, R::RotMatrix{T}) = rot_to_quat(R)

# allow converting element types
convert_rotation{T, U}(::Type{Quaternion{T}}, X::Quaternion{U}) = convert(Quaternion{T}, X)

# more general
convert_rotation{T,U}(::Type{T}, q::Quaternion{U}) = convert_rotation(T, quat_to_rot(promote_eltype(q, eltype(T))))  # go via the rotation representation

# allow construction from a 3 element fixed size array as well, idk about this but it can be done with a Vector so...
convert{T}(::Type{Quaternion}, X::Vec{3,T}) = Quaternion(0, X[1], X[2], X[3])
convert{T}(::Type{Quaternion{T}}, X::Vec{3,T}) = Quaternion(T(0), X[1], X[2], X[3])
convert{T, U}(::Type{Quaternion{T}}, X::Vec{3,U}) = Quaternion(T(0), T(X[1]), T(X[2]), T(X[3]))

# default parameters
default_params{T <: Quaternion}(::Type{T}) = (@DefaultElType(), ) 

# Quaternions have an extra bool field that type_methods.jl doesn't want to know about
numel(::Type{Quaternion}) = 4 



###################################################
# Stereographic projection of Quaternions
###################################################

@doc """

A 3 element parametrization of a unit quaternion Q formed by the intersection of a line from [-1,0,0,0] to Q, with a plane containing the origin and with normal direction [1,0,0,0]. This
is a compact representation of rotations where the derivitives of the rotation matrix's elements wrt the SpQuat parameters are rational functions (good for optimization)

See:

    [1] "A Recipe on the Parameterization of Rotation Matrices for Non-Linear Optimization using Quaternions",
    Terzakis, Culverhouse, Bugmann, Sharma, Sutton, 
    MIDAS technical report, 2011
    http://www.tech.plymouth.ac.uk/sme/springerusv/2011/publications_files/Terzakis%20et%20al%202012,%20A%20Recipe%20on%20the%20Parameterization%20of%20Rotation%20Matrices...MIDAS.SME.2012.TR.004.pdf

    Note 1: the singularity (origin) has been moved from [0,0,0,-1] in Ref[1] to [-1,0,0,0], so the 0 rotation quaternion [1,0,0,0] maps to [0,0,0] as opposed of to [1,0,0]. 
    Note 2: the quaternion rotation ambiguity q = -q means there will be a rotation with ||SpQuat|| <= 1 and an equivilent rotation with ||SpQuat|| > 1.  This package
            uses the solution with ||SpQuat|| <= 1 

""" ->
immutable SpQuat{T <: Real}
    x::T
    y::T
    z::T
end


# convert to a known one
convert_rotation{T <: Real}(::Type{Quaternion{T}}, spq::SpQuat{T}) = spquat_to_quat(spq)
convert_rotation{T <: Real}(::Type{SpQuat{T}}, q::Quaternion{T}) = quat_to_spquat(q)

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{SpQuat{T}}, R::RotMatrix{T}) = quat_to_spquat(convert_rotation(Quaternion{T}, R))
convert_rotation{T,U}(::Type{T}, spq::SpQuat{U}) = convert_rotation(T, spquat_to_quat(promote_eltype(spq, T))) # this has template overload problems without the U????

# allow converting element types
convert_rotation{T, U}(::Type{SpQuat{T}}, X::SpQuat{U}) = convert(SpQuat{T}, X)

# default parameters
default_params{T <: SpQuat}(::Type{T}) = (@DefaultElType(), )

# define an inverse for the SPQuat since its trivial.  Inverse in the sense of the corresponding inverse rotation 
inv(X::SpQuat) = SpQuat(-X.x, -X.y, -X.z) 

# define multiplication for SpQuat corresponding to combining rotations
*(lhs::SpQuat, rhs::SpQuat) = SpQuat(Quaternion(lhs) * Quaternion(rhs))

# pull the rotation angle as well
rot_angle(X::SpQuat) = rot_angle(Quaternion(X))



###################################################
# Arbitrary axis rotation representations
# (limited support for now)
###################################################


immutable AngleAxis{T <: Real}  
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T
end

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{Quaternion{T}}, aa::AngleAxis{T}) = arbaxis_to_quat(aa)
convert_rotation{T <: Real}(::Type{AngleAxis{T}}, q::Quaternion{T}) = quat_to_arbaxis(q)

# go via the quaternion representation
convert_rotation{T <: Real}(::Type{AngleAxis{T}}, R::RotMatrix{T}) = quat_to_arbaxis(convert_rotation(Quaternion{T}, R))
convert_rotation{T,U}(::Type{T}, aa::AngleAxis{U}) = convert_rotation(T, arbaxis_to_quat(promote_eltype(aa, eltype(T))))  # this has template overload problems without the U????

# allow converting element types
convert_rotation{T, U}(::Type{AngleAxis{T}}, X::AngleAxis{U}) = convert(AngleAxis{T}, X)

# default parameters
default_params{T <: AngleAxis}(::Type{T}) = (@DefaultElType(), )

# accessors
rot_angle(aa::AngleAxis) = aa.theta - floor((aa.theta+pi) / (2*pi)) * 2*pi  # named to match the quaternion equivilent
axis(aa::AngleAxis) = Vec(aa.axis_x, aa.axis_y, aa.axis_z)

# construct with funny inputs
call(::Type{AngleAxis}, theta::Real, x::Real, y::Real, z::Real) = add_params(AngleAxis, promote_type(typeof(theta), typeof(x), typeof(y), typeof(z)))(theta,x,y,z)




###################################################
# Euler Angles
###################################################

# this file defines the types:
# EulerAngles       - for tait bryan ordering
# ProperEulerAngles - for proper Euler ordering
include("euler_types.jl") 

default_params{T <: EulerAngles}(::Type{T}) = (@DefaultEulerOrder(), @DefaultElType()) 
default_params{T <: ProperEulerAngles}(::Type{T}) = (@DefaultProperEulerOrder(), @DefaultElType()) # there's no reason why this order was chosen as the default

# Rotations to Euler angles
convert_rotation{T <: EulerAngles}(::Type{T}, R::RotMatrix) = rot_to_euler(add_params(T, R), R)
convert_rotation{T <: ProperEulerAngles}(::Type{T}, R::RotMatrix) = rot_to_euler(add_params(T, R), R)

# Euler angles to rotations
convert_rotation{T <: RotMatrix}(::Type{T}, ea::EulerAngles) = euler_to_rot(add_params(RotMatrix, ea), ea)
convert_rotation{T <: RotMatrix}(::Type{T}, ea::ProperEulerAngles) = euler_to_rot(add_params(RotMatrix, ea), ea)

# Euler angles to other
convert_rotation{T,ORD,U}(::Type{T}, ea::EulerAngles{ORD, U}) = convert_rotation(add_params(T, ea), euler_to_rot(RotMatrix, promote_eltype(ea, eltype(T))))
convert_rotation{T,ORD,U}(::Type{T}, ea::ProperEulerAngles{ORD, U}) = convert_rotation(add_params(T, ea), euler_to_rot(RotMatrix, promote_eltype(ea, eltype(T))))

# allow converting element types
convert_rotation{ORD, T, U}(::Type{EulerAngles{ORD, T}}, X::EulerAngles{ORD, U}) = convert(EulerAngles{ORD, T}, X)
convert_rotation{ORD, T, U}(::Type{ProperEulerAngles{ORD, T}}, X::ProperEulerAngles{ORD, U}) = convert(ProperEulerAngles{ORD, T}, X)
















