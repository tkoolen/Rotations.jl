

# a default element type
DefaultElType = Float64

###################################################
# Rotation matrix (typeofalias of Mat{3,3,Float64})
# N.B. all other parameterizations T should implement
# a convert_rotation(T, RotMatrix) method
###################################################

@doc """
A type alias for rotation matrices (alias of Mat{3,3,T})
""" ->
typealias RotMatrix{T <: AbstractFloat}   Mat{3,3,T}
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
default_params{T <: RotMatrix}(::Type{T}) = (DefaultElType, ) 


###################################################
# Quaternions
# (Quaternions are defined in Quaternions.jl)
###################################################

# element type is handy
eltype(::Type{Quaternion}) = Any
eltype{T}(::Type{Quaternion{T}}) = T


# add an indexing scheme
getindex(X::Quaternion, idx::Integer) = getfield(X, idx)

# an extra constructor (its annoying this is missing)
call{T}(::Type{Quaternion{T}}, a::Real, b::Real, c::Real, d::Real) = Quaternion(T(a), T(b), T(c), T(d))

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{RotMatrix{T}}, q::Quaternion{T}) = quattorot(q)
convert_rotation{T <: Real}(::Type{Quaternion{T}}, R::RotMatrix{T}) = rottoquat(R)

# allow converting element types
convert_rotation{T, U}(::Type{Quaternion{T}}, X::Quaternion{U}) = convert(Quaternion{T}, X)

# more general
convert_rotation{T,U}(::Type{T}, q::Quaternion{U}) = convert_rotation(T, quattorot(promote_eltype(q, eltype(T))))  # go via the rotation representation

# default parameters
default_params{T <: Quaternion}(::Type{T}) = (DefaultElType, ) 




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

    Note the singularity (origin) has been moved from [0,0,0,-1] in Ref[1] to [-1,0,0,0], so the 0 rotation quaternion [1,0,0,0] maps to [0,0,0] as opposed of to [1,0,0]. 
    This is done because Spquats have 2 norm <= 1, so [1,0,0] lies on an inequality constraint

""" ->
immutable SpQuat{T <: AbstractFloat} <: FixedVectorNoTuple{3, T}
    x::T
    y::T
    z::T
end

# convert to a known one
convert_rotation{T <: Real}(::Type{Quaternion{T}}, spq::SpQuat{T}) = spquattoquat(spq)
convert_rotation{T <: Real}(::Type{SpQuat{T}}, q::Quaternion{T}) = quattospquat(q)

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{SpQuat{T}}, R::RotMatrix{T}) = quattospquat(convert_rotation(Quaternion{T}, R))
convert_rotation{T,U}(::Type{T}, spq::SpQuat{U}) = convert_rotation(T, spquattoquat(promote_eltype(spq, T))) # this has template overload problems without the U????

# allow converting element types
convert_rotation{T, U}(::Type{SpQuat{T}}, X::SpQuat{U}) = convert(SpQuat{T}, X)

# default parameters
default_params{T <: SpQuat}(::Type{T}) = (DefaultElType, )



###################################################
# Arbitrary axis rotation representations
# (limited support for now)
###################################################


immutable AngleAxis{T <: AbstractFloat}  <: FixedVectorNoTuple{4, T}
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T
end

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{Quaternion{T}}, aa::AngleAxis{T}) = arbaxistoquat(a)
convert_rotation{T <: Real}(::Type{AngleAxis{T}}, q::Quaternion{T}) = quattoarbaxis(q)

# go via the quaternion representation
convert_rotation{T <: Real}(::Type{AngleAxis{T}}, R::RotMatrix{T}) = quattoarbaxis(convert_rotation(Quaternion{T}, R))
convert_rotation{T,U}(::Type{T}, aa::AngleAxis{U}) = convert_rotation(T, arbaxistoquat(promote_eltype(aa, eltype(T))))  # this has template overload problems without the U????

# allow converting element types
convert_rotation{T, U}(::Type{AngleAxis{T}}, X::AngleAxis{U}) = convert(AngleAxis{T}, X)

# default parameters
default_params{T <: AngleAxis}(::Type{T}) = (DefaultElType, )

# accessors
rot_angle(aa::AngleAxis) = aa.theta  # named to match the quaternion equivilent
axis(aa::AngleAxis) = Vec(aa.axis_x, aa.axis_y, aa.axis_z)




###################################################
# Euler Angles
###################################################

# this file defines the types:
# EulerAngles       - for tait bryan ordering
# ProperEulerAngles - for proper Euler ordering
include("euler_types.jl") 

DefaultEulerOrder = EulerZXY       # there's no reason why this order was chosen as the default
DefaultProperEulerOrder = EulerXZX # there's no reason why this order was chosen as the default

default_params{T <: EulerAngles}(::Type{T}) = (DefaultEulerOrder, DefaultElType) 
default_params{T <: ProperEulerAngles}(::Type{T}) = (DefaultProperEulerOrder, DefaultElType) # there's no reason why this order was chosen as the default

# add default values to them
call(::Type{EulerAngles}, x::Integer, y::Integer, z::Integer) = add_params(EulerAngles)(x,y,z)
call(::Type{EulerAngles}, x::Real, y::Real, z::Real) = add_params(EulerAngles, promote_type(typeof(x), typeof(y), typeof(z)))(x,y,z)

call(::Type{ProperEulerAngles}, x::Integer, y::Integer, z::Integer) = add_params(ProperEulerAngles)(x,y,z)
call(::Type{ProperEulerAngles}, x::Real, y::Real, z::Real) = add_params(ProperEulerAngles, promote_type(typeof(x), typeof(y), typeof(z)))(x,y,z)

# add some extra construction methods
call{T <: TaitByranOrder}(::Type{EulerAngles{T}}, x::Integer, y::Integer, z::Integer) = add_params(EulerAngles{T})(x,y,z)
call{T <: TaitByranOrder}(::Type{EulerAngles{T}}, x::Real, y::Real, z::Real) = add_params(EulerAngles{T}, promote_type(typeof(x), typeof(y), typeof(z)))(x,y,z)

call{T <: ProperEulerOrder}(::Type{ProperEulerAngles{T}}, x::Integer, y::Integer, z::Integer) = add_params(ProperEulerAngles{T})(x,y,z)
call{T <: ProperEulerOrder}(::Type{ProperEulerAngles{T}}, x::Real, y::Real, z::Real) = add_params(ProperEulerAngles{T}, promote_type(typeof(x), typeof(y), typeof(z)))(x,y,z)


# Rotations to Euler angles
convert_rotation{T <: EulerAngles}(::Type{T}, R::RotMatrix) = rottoeuler(add_params(T, R), R)
convert_rotation{T <: ProperEulerAngles}(::Type{T}, R::RotMatrix) = rottoeuler(add_params(T, R), R)

# Euler angles to rotations
convert_rotation{T <: RotMatrix}(::Type{T}, ea::EulerAngles) = eulertorot(add_params(RotMatrix, ea), ea)
convert_rotation{T <: RotMatrix}(::Type{T}, ea::ProperEulerAngles) = eulertorot(add_params(RotMatrix, ea), ea)

# Euler angles to other
convert_rotation{T,ORD,U}(::Type{T}, ea::EulerAngles{ORD, U}) = convert_rotation(add_params(T, ea), eulertorot(RotMatrix, promote_eltype(ea, eltype(T))))
convert_rotation{T,ORD,U}(::Type{T}, ea::ProperEulerAngles{ORD, U}) = convert_rotation(add_params(T, ea), eulertorot(RotMatrix, promote_eltype(ea, eltype(T))))

# allow converting element types
convert_rotation{ORD, T, U}(::Type{EulerAngles{ORD, T}}, X::EulerAngles{ORD, U}) = convert(EulerAngles{ORD, T}, X)
convert_rotation{ORD, T, U}(::Type{ProperEulerAngles{ORD, T}}, X::ProperEulerAngles{ORD, U}) = convert(ProperEulerAngles{ORD, T}, X)
















