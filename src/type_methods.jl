# defines common methods for the various rotation types

# It'd be nice if this was an abstract type
TypeList       = [RotMatrix, Quaternion, SpQuat, EulerAngles, AngleAxis]
RotationTypes  = Union{TypeList...}
MatTypes       = Union{RotMatrix}
VectorTypes    = Union{Quaternion, SpQuat, EulerAngles, AngleAxis}
Mat3x3Types    = Union{RotMatrix}
Vec4Types      = Union{Quaternion, AngleAxis}
Vec3Types      = Union{SpQuat, EulerAngles}
VecNoFSAType   = Union{Quaternion}


#################################################################
# force the output type to have all of its template parameters
# N.B. these need to be generated for type stability because 
# T and U need to be known to the compiler to detemine the output type
#################################################################

DefaultReal = Float64

@generated add_param{T <: RotationTypes}(::Type{T}) = add_param(T, DefaultReal)  # supplying a default arg to the below function broke the types
@generated add_param{T <: RotationTypes, U <: Real}(::Type{T}, ::Type{U}) = isa(T.parameters[end], DataType) ? :(T) : :(T{$(U)})  # use Float64 by default

add_param{T <: RotationTypes, U <: RotationTypes}(::Type{T}, ::U) = add_param(T, U)
@generated add_param{T <: RotationTypes, U <: RotationTypes}(::Type{T}, ::Type{U}) = isa(T.parameters[end], DataType) ? :(T) : :(T{$(eltype(U))})

add_param{T <: EulerAngles, U <: RotationTypes}(::Type{T}, ::U) = add_param(T, U)
@generated add_param{T <: EulerAngles, U <: RotationTypes}(::Type{T}, ::Type{U}) = add_param(T, U.parameters[end])

@generated add_param{T <: EulerAngles}(::Type{T}) = add_param(T, DefaultReal)
@generated function add_param{T <: EulerAngles, U <: Real}(::Type{T}, ::Type{U})
	if isa(T.parameters[1], TypeVar) 
		:(T{$(DefaultEulerOrder), $(U)} )
	elseif isa(T.parameters[2], TypeVar) 
		:(EulerAngles{$(T.parameters[1]), $(U)})
	else
		:(T)
	end
end

# function for element type promotions
promote_type_sp{T <: Real}(::Type{Any}, ::Type{T}) = T
promote_type_sp{T <: Real}(::Type{T}, ::Type{Any}) = T
promote_type_sp{T <: Real, U <: Real}(::Type{T}, ::Type{U}) = promote_type(T, U)

promote_eltype{T <: RotationTypes, U <: RotationTypes}(::Type{T}, ::Type{U}) = promote_type_sp(eltype(U), eltype(V))

promote_rottype{T <: RotationTypes}(::Type{T}, ::Type{Any}) = T
@generated function promote_rottype{T <: RotationTypes, U <: Real}(::Type{T}, ::Type{U})
	elT = promote_type_sp(eltype(T), U)
	:($(strip_eltype(T){elT}))
end

promote_rottype{T <: RotationTypes}(X::T, ::Type{Any}) = X
promote_rottype{T <: RotationTypes, U <: Real}(X::T, ::Type{U}) = convert(promote_rottype(T, U), X)
promote_rottype{T <: RotationTypes, U <: RotationTypes}(X::T, ::Type{U}) = promote_rottype(X, eltype(U))



################################################################
# Extra constructors
################################################################

# allow constuction from any rotation parameterization
call{T <: RotationTypes, U <: RotationTypes}(::Type{T}, X::U) = convert_rotation(add_param(T, U), X)


################################################################
# Allow "convert_rotation" to perform element type conversions
################################################################





################################################################
# conversion to / from fixed size arrays and Vectors/ Matrices
################################################################

# convert to immutable vector
convert{T <: VecNoFSAType}(::Type{Vec}, X::T) = convert(Vec{ T <: Vec4Types ? 4 : 3, eltype{T}}, X)
convert{T <: Real, U <: Vec3Types}(::Type{Vec{3,T}}, X::U) = Vec{3,T}(T(X[1]),T(X[2]),T(X[3]))
convert{T <: Real, U <: Quaternion}(::Type{Vec{4,T}}, X::U) = Vec{4,T}(T(X[1]),T(X[2]),T(X[3]),T(X[4]))

# convert to mutable vector
convert{T <: VectorTypes}(::Type{Vector}, X::T) = convert(Vector{eltype{T}}, X)
convert{T <: Real, U <: Vec3Types}(::Type{Vector{T}},X::U) = vcat(T(X[1]),T(X[2]),T(X[3]))
convert{T <: Real, U <: Vec4Types}(::Type{Vector{T}},X::U) = vcat(T(X[1]),T(X[2]),T(X[3]),T(X[4]))

# convert from immutable arrays
function convert{T <: Vec4Types, U <: Real}(::Type{T}, v::Vec{4, U})
	oT = add_param(T, U)
	elT = eltype(oT)
	els == oT(elT(v[1]),elT(v[2]),elT(v[3]),elT(v[4]))
end
function convert{T <: Vec3Types, U <: Real}(::Type{T}, v::Vec{3, U})
	oT = add_param(T, U)
	elT = eltype(oT)
	els == oT(elT(v[1]),elT(v[2]),elT(v[3]))
end

# convert from mutable vector
function convert{T <: VectorTypes, U <: Real}(::Type{T}, v::Vector{U})
	els = T <: Vec4Types ? 4 : 3
	(length(v)) == els || error("Input vector should have length $(els) for conversion to a $(T)")
	oT = add_param(T, U)
	elT = eltype(oT)
	els == T <: Vec4Types ? oT(elT(v[1]),elT(v[2]),elT(v[3]),elT(v[4])) : oT(elT(v[1]),elT(v[2]),elT(v[3]))
end























################################################################
# And some extra methods
################################################################

# add isanan to fixed size arrays (convenient)
isnan{T <: Mat3x3Types}(R::T) = @fsa_isnan(R, 3, 3)
isnan{T <: Vec4Types}(V::T) = @fsa_isnan_vec(V, 4)
isnan{T <: Vec3Types}(V::T) = @fsa_isnan_vec(V, 3)


@doc """
function to convert an immutable transformation matrix from RN to PN
"""  ->
function projective(Tmat::RotMatrix{Float64})
	return @fsa_projective(Tmat, 3, 3, Float64)
end

