###################################################
# Euler Representations
###################################################


@doc """
EulerOrder specifies the extrinsic (axes of rotation are fixed) order of rotations.
""" -> 
abstract EulerOrder # for order of rotations


# add other orders here

@doc """
EulerZXY corresponds to:

	Extrinsic: R = R(z) * R(x) * R(y)      

	or equivilently,

    Intrinsic: R = R(y'') * R(x') * R(z)   (where axes of rotation inherit previous rotations)
""" ->
type EulerZXY <: EulerOrder end



@doc """
EulerAngles{EulerOrder, T} where order specifies the extrinsic (axes of rotation are fixed) order of rotations.

Example: order EulerZXY corresponds to:

	Extrinsic: R = R(z) * R(x) * R(y)      

	or equivilently,

    Intrinsic: R = R(y'') * R(x') * R(z)   (where axes of rotation inherit previous rotations)

Check Rotation.EulerOrders for implemented orders

Notes:

	The fields for of this structure (x,y,z) are setup for conceptual simplicity for Trait Bryan angles which
	use rotationa each axes.  Proper Euler angles use repeat rotations about a single axis (e.g. EulerZXZ)

""" ->
immutable EulerAngles{Order <: EulerOrder, T <: Real} <: FixedVectorNoTuple{3, T}  # be nice to get quaternions in here
	x::T
	y::T
	z::T
end

# a default ordering scheme
DefaultEulerOrder = EulerZXY

# get the order from the type
eltype{Order <: EulerOrder, T}(::Type{EulerAngles{Order,T}}) = T
euler_order{Order <: EulerOrder, T}(::Type{EulerAngles{Order,T}}) = Order



# define its interaction with other angle representations
convert_rotation{RTYPE <: RotMatrix}(::Type{RTYPE}, ea::EulerAngles) = convert(RTYPE, eulertorot(ea))
convert_rotation{EA <: EulerAngles}(::Type{EA}, R::RotMatrix) = convert(add_params(EA, R), rottoeuler(EA, R))


convert_rotation{ETYPE <: EulerAngles}(::Type{ETYPE}, R::RotMatrix) = convert(ETYPE, rottoeuler(ETYPE, R))
convert_rotation{T}(::Type{T}, ea::EulerAngles) = convert_rotation(T, eulertorot(ea))  # go via the rotation representation

# allow convert_rotation to also convert element types (its convenient)
convert_rotation{ORD <: EulerOrder, T <: Real}(::Type{EulerAngles{ORD,T}}, ea::EulerAngles{ORD}) = convert(EulerAngles{ORD,T}, ea) 
convert_rotation{ORD <: EulerOrder}(::Type{EulerAngles{ORD}}, ea::EulerAngles{ORD}) = ea  # no conversion needed                    
