
###################################################
# Rotation matrix (typeofalias of Mat{3,3,Float64})
# N.B. all other parameterizations T should implement
# a convert_rotation(T, RotMatrix) method
###################################################

@doc """
A type alias for rotation matrices (alias of Mat{3,3,T})
""" ->
typealias RotMatrix{T <: Real}   Mat{3,3,T}

# go back to unparameterized type
strip_eltype{T <: RotMatrix}(::Type{T}) = RotMatrix
 
# fix weirdness, try: mat = Mat(eye(3)); @code_warntype convert(typeof(mat), mat)  
convert{T <: RotMatrix}(::Type{T}, X::T) = X

# RotMatrix(eye(3)) hangs, so fix it
# can't use the union types because it wont replace whats in fixed size arrays
Base.convert{U}(::Type{RotMatrix}, mat::Matrix{U}) = RotMatrix{U}(mat)
Base.call{U}(::Type{RotMatrix}, mat::Matrix{U}) = convert(RotMatrix{U}, mat)



###################################################
# Quaternions
# (Quaternions are defined in Quaternions.jl)
###################################################

# element type is handy
eltype(::Type{Quaternion}) = Any
eltype{T}(::Type{Quaternion{T}}) = T
strip_eltype{T <: Quaternion}(::Type{T}) = Quaternion

# add an indexing scheme
getindex(X::Quaternion, idx::Integer) = getfield(X, idx)

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{RotMatrix{T}}, q::Quaternion{T}) = quattorot(q)
convert_rotation{T <: Real}(::Type{Quaternion{T}}, R::RotMatrix{T}) = rottoquat(R)
convert_rotation{T}(::Type{T}, q::Quaternion) = convert_rotation(T, quattorot(promote_rottype(q, eltype(T))))  # go via the rotation representation




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
immutable SpQuat{T <: Real} <: FixedVectorNoTuple{3, T}
	x::T
	y::T
	z::T
end

strip_eltype{T <: SpQuat}(::Type{T}) = SpQuat

# convert to a known one
convert_rotation{T <: Real}(::Type{Quaternion{T}}, spq::SpQuat{T}) = spquattoquat(spq)
convert_rotation{T <: Real}(::Type{SpQuat{T}}, q::Quaternion{T}) = quattospquat(q)

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{SpQuat{T}}, R::RotMatrix{T}) = quattospquat(convert_rotations(Quaternion{T}, R))
convert_rotation{T}(::Type{T}, spq::SpQuat) = convert_rotation(T, spquattoquat(promote_rottype(spq, eltype(T))))






###################################################
# Arbitrary axis rotation representations
# (limited support for now)
###################################################


immutable AngleAxis{T <: Real}  <: FixedVectorNoTuple{4, T}
	theta::T
	axis_x::T
	axis_y::T
	axis_z::T
end

strip_eltype{T <: AngleAxis}(::Type{T}) = AngleAxis

# define its interaction with other angle representations
convert_rotation{T <: Real}(::Type{Quaternion{T}}, aa::AngleAxis{T}) = arbaxistoquat(a)
convert_rotation{T <: Real}(::Type{AngleAxis{T}}, q::Quaternion{T}) = quattoarbaxis(q)

# go via the quaternion representation
convert_rotation{T <: Real}(::Type{AngleAxis{T}}, R::RotMatrix{T}) = quattosarbaxis(convert_rotations(Quaternion{T}, R))
convert_rotation{T}(::Type{T}, aa::AngleAxis) = convert_rotation(T, arbaxistoquat(promote_rottype(aa, eltype(T))))  













