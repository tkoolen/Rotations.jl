"""
    abstract Rotation{T} <: StaticMatrix{T}

An abstract type representing 3D rotations. More abstractly, they represent
3×3 unitary (orthogonal) matrices.
"""
abstract Rotation{T} <: StaticMatrix{T}

Base.size{R <: Rotation}(::Union{R,Type{R}}) = (3,3)
Base.ctranspose(r::Rotation) = inv(r)
Base.transpose{T<:Real}(r::Rotation{T}) = inv(r)

# Rotation matrices should be orthoginal/unitary. Only the operations we define,
# like multiplication, will stay as Rotations, otherwise users will get an
# SMatrix{3,3} (e.g. rot1 + rot2 -> SMatrix)
Base.@pure StaticArrays.similar_type{R <: Rotation}(::Union{R,Type{R}}) = SMatrix{3, 3, eltype(R)}
Base.@pure StaticArrays.similar_type{R <: Rotation, T}(::Union{R,Type{R}}, ::Type{T}) = SMatrix{3, 3, T}

################################################################################
################################################################################
"""
    immutable RotMatrix{T} <: Rotation{T}

A statically-sized, 3×3 unitary (orthogonal) matrix.

Note: the orthonormality of the input matrix is not checked by the constructor.
"""
immutable RotMatrix{T} <: Rotation{T} # which is <: AbstractMatrix{T}
    mat::SMatrix{3, 3, T, 9} # The final parameter to SMatrix is the "length" of the matrix, 3 × 3 = 9
end

# These 2 functions (plus size) are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RM}){RM<:RotMatrix}(t::NTuple{9}) = RM(SMatrix{3,3}(t))
Base.@propagate_inbounds Base.getindex(r::RotMatrix, i::Integer) = r.mat[i]

# A rotation is more-or-less defined as being an orthogonal (or unitary) matrix
Base.inv(r::RotMatrix) = RotMatrix(r.mat')

# A useful constructor for identity rotation (eye is already provided by StaticArrays, but needs an eltype)
@inline Base.eye(::Type{RotMatrix}) = eye(RotMatrix{Float64})

# By default, composition of rotations will go through RotMatrix, unless overridden
@inline *(r1::Rotation, r2::Rotation) = RotMatrix(r1) * RotMatrix(r2)
@inline *(r1::RotMatrix, r2::Rotation) =r1 * RotMatrix(r2)
@inline *(r1::Rotation, r2::RotMatrix) = RotMatrix(r1) * r2
@inline *(r1::RotMatrix, r2::RotMatrix) = RotMatrix(r1.mat * r2.mat) # TODO check that this doesn't involve extra copying.

# Removes module name from output, to match other types
function Base.summary(r::RotMatrix)
    "3×3 RotMatrix{$(eltype(r))})"
end

################################################################################
################################################################################
"""
    Quat{T} <: Rotation{T}
    Quat(w, x, y, z)

The `Quat` type is a 3×3 matrix representation of a normalized quaternion.
They allow you to transparently use (fast) quaternion algebra to store, compose
and invert 3D rotations, while at the same time letting you apply rotations
through matrix-vector multiplication.

Note: the constructor will always renormalize the input so that the quaternion
has length 1 (w² + x² + y² + z² = 1), and the rotation matrix is orthogonal.
"""
immutable Quat{T} <: Rotation{T}
    w::T
    x::T
    y::T
    z::T

    # For the moment we ensure that it is normalized upon construction.
    function Quat(w, x, y, z)
        norm = copysign(sqrt(w*w + x*x + y*y + z*z), w)
        # Should this be an error or warning, if it isn't approximately normalized? E.g.:
        #if norm !≈ 1
        #    error("Expected a normalized quaternion") # or warn() ?
        #end
        new(w/norm, x/norm, y/norm, z/norm)
    end
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline (::Type{Quat}){W,X,Y,Z}(w::W, x::X, y::Y, z::Z) = Quat{promote_type(promote_type(promote_type(W, X), Y), Z)}(w, x, y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
function (::Type{RQ}){RQ<:Quat}(t::NTuple{9})
    q = Quat(sqrt(abs(1  + t[1] + t[5] + t[9])) / 2,
                copysign(sqrt(abs(1 + t[1] - t[5] - t[9]))/2, t[6] - t[8]),
                copysign(sqrt(abs(1 - t[1] + t[5] - t[9]))/2, t[7] - t[3]),
                copysign(sqrt(abs(1 - t[1] - t[5] + t[9]))/2, t[2] - t[4]))
end

function Base.getindex(q::Quat, i::Integer)
    if i == 1
        ww = (q.w * q.w)
        xx = (q.x * q.x)
        yy = (q.y * q.y)
        zz = (q.z * q.z)

        ww + xx - yy - zz
    elseif i == 2
        xy = (q.x * q.y)
        zw = (q.w * q.z)

        2 * (xy + zw)
    elseif i == 3
        xz = (q.x * q.z)
        yw = (q.y * q.w)

        2 * (xz - yw)
    elseif i == 4
        xy = (q.x * q.y)
        zw = (q.w * q.z)

        2 * (xy - zw)
    elseif i == 5
        ww = (q.w * q.w)
        xx = (q.x * q.x)
        yy = (q.y * q.y)
        zz = (q.z * q.z)

        ww - xx + yy - zz
    elseif i == 6
        yz = (q.y * q.z)
        xw = (q.w * q.x)

        2 * (yz + xw)
    elseif i == 7
        xz = (q.x * q.z)
        yw = (q.y * q.w)

        2 * (xz + yw)
    elseif i == 8
        yz = (q.y * q.z)
        xw = (q.w * q.x)

        2 * (yz - xw)
    elseif i == 9
        ww = (q.w * q.w)
        xx = (q.x * q.x)
        yy = (q.y * q.y)
        zz = (q.z * q.z)

        ww - xx - yy + zz
    else
        throw(BoundsError(r,i))
    end
end

# This speeds up some StaticArray methods (such as conversion to RotMatrix)
function Base.convert(::Type{Tuple}, q::Quat)
    ww = (q.w * q.w)
    xx = (q.x * q.x)
    yy = (q.y * q.y)
    zz = (q.z * q.z)
    xy = (q.x * q.y)
    zw = (q.w * q.z)
    xz = (q.x * q.z)
    yw = (q.y * q.w)
    yz = (q.y * q.z)
    xw = (q.w * q.x)

    # initialize rotation part
    return (ww + xx - yy - zz,
            2 * (xy + zw),
            2 * (xz - yw),
            2 * (xy - zw),
            ww - xx + yy - zz,
            2 * (yz + xw),
            2 * (xz + yw),
            2 * (yz - xw),
            ww - xx - yy + zz)
end

function Base.:*(q::Quat, x::AbstractVector)
    if length(X) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(X))")
    end

    qo = (-q.x * X[1] - q.y * X[2] - q.z * X[3],
           q.w * X[1] + q.y * X[3] - q.z * X[2],
           q.w * X[2] - q.x * X[3] + q.z * X[1],
           q.w * X[3] + q.x * X[2] - q.y * X[1])

    Xo = SVector(-qo[1] * q.x + qo[2] * q.w - qo[3] * q.z + qo[4] * q.y,
                 -qo[1] * q.y + qo[2] * q.z + qo[3] * q.w - qo[4] * q.x,
                 -qo[1] * q.z - qo[2] * q.y + qo[3] * q.x + qo[4] * q.w)
end

# TODO consider Ac_mul_B, etc...
function Base.:*(q1::Quat, q2::Quat)
    Quat(q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z,
            q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y,
            q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x,
            q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w)
end

function Base.inv(q::Quat)
    Quat(q.w, -q.x, -q.y, -q.z)
end

@inline Base.eye(::Type{Quat}) = Quat(1.0, 0.0, 0.0, 0.0)
@inline Base.eye{T}(::Type{Quat{T}}) = Quat{T}(one(T), zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(q::Quat)
    "3×3 Quat{$(eltype(q))}($(q.w), $(q.x), $(q.y), $(q.z))"
end

#=
function Base.show(io::IO, mime::(MIME"text/plain"), q::Quat)
    if get(io, :compact, false)
        print(io, "Quat(")
        print(io, q.w) # Shorter printing of numbers than "$(q.w)"
        print(io, ", ")
        print(io, q.x)
        print(io, ", ")
        print(io, q.y)
        print(io, ", ")
        print(io, q.z)
        print(io, ")")
    else
        # Force the matrix-header to include quaternion parameters.
        println(io, "3×3 Quat{$(eltype(q))}($(q.w), $(q.x), $(q.y), $(q.z)):")
        Base.print_matrix(IOContext(io, :compact => true), q)
        println(io)
    end
end
=#

################################################################################
################################################################################
"""
    immutable SPQuat{T} <: Rotation{T}
    SPQuat(x, y, z)

A SPQuat is a 3D rotation represented by the "stereographic projection" of a normalized quaternion (shortened to "SPQuat"), which is
a 3-element parametrization of a unit quaternion Q formed by the intersection of a line from [-1,0,0,0] to Q, with a plane containing the origin and with normal direction [1,0,0,0]. This
is a compact representation of rotations where the derivitives of the rotation matrix's elements w.r.t. the SPQuat parameters are rational functions (making them useful for optimization).

See:

    [1] "A Recipe on the Parameterization of Rotation Matrices for Non-Linear Optimization using Quaternions",
    Terzakis, Culverhouse, Bugmann, Sharma, Sutton,
    MIDAS technical report, 2011
    http://www.tech.plymouth.ac.uk/sme/springerusv/2011/publications_files/Terzakis%20et%20al%202012,%20A%20Recipe%20on%20the%20Parameterization%20of%20Rotation%20Matrices...MIDAS.SME.2012.TR.004.pdf

    Note 1: the singularity (origin) has been moved from [0,0,0,-1] in Ref[1] to [-1,0,0,0], so the 0 rotation quaternion [1,0,0,0] maps to [0,0,0] as opposed of to [1,0,0].
    Note 2: the quaternion rotation ambiguity q = -q means there will be a rotation with ||SPQuat|| <= 1 and an equivilent rotation with ||SPQuat|| > 1.  This package
            uses the solution with ||SPQuat|| <= 1
    Note 3: it is safe to assume that the corresponding matrix is orthogonal/unitary for any input x, y, z.

"""
immutable SPQuat{T} <: Rotation{T}
    x::T
    y::T
    z::T

    # TODO should we enforce norm <= 1?
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{SPQuat}){X,Y,Z}(x::X, y::Y, z::Z) = SPQuat{promote_type(promote_type(X, Y), Z)}(x, y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{SPQuat})(t::NTuple{9}) = SPQuat(Quat(t))
@inline Base.getindex(spq::SPQuat, i::Integer) = Quat(spq)[i]

@inline function Base.convert{Q <: Quat}(::Type{Q}, spq::SPQuat)
    # Both the sign and norm of the Quat is automatically dealt with in its inner constructor
    return Q(1 - (spq.x*spq.x + spq.y*spq.y + spq.z*spq.z), 2*spq.x, 2*spq.y, 2*spq.z)
end

@inline function Base.convert{SPQ <: SPQuat}(::Type{SPQ}, q::Quat)
    alpha2 = (1 - q.w) / (1 + q.w) # <= 1 since q.w >= 0
    spq = SPQ(q.x * (alpha2 + 1)*0.5,  q.y * (alpha2 + 1)*0.5, q.z * (alpha2 + 1)*0.5)
end

@inline Base.convert(::Type{Tuple}, spq::SPQuat) = Tuple(Quat(spq))

@inline Base.:*(spq::SPQuat, x::AbstractVector) = Quat(spq) * x

@inline Base.:*(spq::SPQuat, r::Rotation) = Quat(spq) * r
@inline Base.:*(r::Rotation, spq::SPQuat) = r * Quat(spq)
@inline Base.:*(spq1::SPQuat, spq2::SPQuat) = SPQuat(Quat(spq1) * Quat(spq2))

@inline Base.inv(spq::SPQuat) = SPQuat(-spq.x, -spq.y, -spq.z)

@inline Base.eye(::Type{SPQuat}) = SPQuat(0.0, 0.0, 0.0)
@inline Base.eye{T}(::Type{SPQuat{T}}) = SPQuat{T}(zero(T), zero(T), zero(T))

# rotation properties
@inline rotation_angle(spq::SPQuat) = rotation_angle(Quaternion(spq))
@inline rotation_axis(spq::SPQuat) = rotation_axis(Quaternion(spq))

# Need SPQuat parameters in the output representation
function Base.summary(spq::SPQuat)
    "3×3 SPQuat{$(eltype(spq))}($(spq.x), $(spq.y), $(spq.z))"
end

#=
function Base.show(io::IO, mime::(MIME"text/plain"), spq::SPQuat)
    if get(io, :compact, false)
        print(io, "SPQuat(")
        print(io, spq.x)
        print(io, ", ")
        print(io, spq.y)
        print(io, ", ")
        print(io, spq.z)
        print(io, ")")
    else
        # Force the matrix-header to include quaternion parameters.
        println(io, "3×3 SPQuat{$(eltype(spq))}($(spq.x), $(spq.y), $(spq.z)):")
        Base.print_matrix(IOContext(io, :compact => true), spq)
        println(io)
    end
end
=#








#=
#import Base: convert, getindex, *, eltype, eye, norm, inv
#import Quaternions: axis, angle

@deprecate(rot_angle, rotation_angle)

# defines common methods for the various rotation types
RotTypeList = Vector{Type}(0)

# build a dictionary to build transformation paths from one representation to the next
# I hate this but trying to make a generic overload for convert causes no end of problems...
conversion_path = Dict()
defined_conversions = Vector{NTuple{2,Type}}(0)  # use this to keep track of all conversions that have been defined (fill in the blanks later)


numel{T}(::Type{T}) = length(fieldnames(T))  # for getting the number of elements in the representation
                                             # e.g. length(Quaternion) == 1                # its a number
                                             # e.g. length(fieldnames(Quaternion)) == 5    # it has a normalized field
                                             # e.g. numel(Quaternion) == 4                 # what we want



"""
A type alias for rotation matrices (alias of Mat{3,3,T})
"""
typealias RotMatrix{T}   Mat{3,3,T}

# add to the type list
push!(RotTypeList, RotMatrix)

# Fixes for FixedSizeArray stuff (based on version 0.1.0)
@compat @inline (::Type{RotMatrix}){T}(mat::RotMatrix{T}) = mat    # weird result otherwise
@inline convert{T}(::Type{RotMatrix}, mat::Matrix{T}) = mat          # this hangs othewise
@compat @inline (::Type{RotMatrix}){T}(mat::Matrix{T}) = convert(RotMatrix{T}, mat)    # this hangs othewise

@compat @inline (::Type{RotMatrix}){T}(xt::NTuple{9,T}) = RotMatrix{T}(xt[1:3], xt[4:6], xt[7:9])

@compat @inline (::Type{RotMatrix}){T}(c1::NTuple{3,T}, c2::NTuple{3,T}, c3::NTuple{3,T}) = RotMatrix{T}(c1, c2, c3)
@compat @inline (::Type{RotMatrix}){T1,T2,T3}(c1::NTuple{3,T1}, c2::NTuple{3,T2}, c3::NTuple{3,T3}) = RotMatrix{promote_type(T1,T2,T3)}(c1, c2, c3)

@compat (::Type{RotMatrix}){T}(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T) = RotMatrix{T}((x1,x2,x3), (x4,x5,x6), (x7,x8,x9)) # Fixed at v0.0.1, but still needed so the below works
@compat @inline (::Type{RotMatrix{T}}){T}(x1, x2, x3, x4, x5, x6, x7, x8, x9) =
                    RotMatrix((T(x1),T(x2),T(x3)), (T(x4),T(x5),T(x6)), (T(x7),T(x8),T(x9)))

conversion_path[RotMatrix] = nothing  # the blessed type, everything goes through RotMatrixs by default

# number of numeric elements
numel(::Type{RotMatrix}) = 9

strip_eltype{T <: RotMatrix}(::Type{T}) = RotMatrix

# define null rotations for convenience
@inline eye(::Type{RotMatrix}) = eye(RotMatrix{Float64})
=#
#=

###################################################
# Quaternions
# (Quaternions are defined in Quaternions.jl)
###################################################

# add to the type list
push!(RotTypeList, Quaternion)

# an extra constructor so we don't need to specify whether its normalized
@compat @inline (::Type{Quaternion{T}}){T}(a, b, c, d) = Quaternion(T(a), T(b), T(c), T(d))

# define its interaction with RotMatrixs
@inline convert(::Type{RotMatrix}, q::Quaternion) = quat_to_rot(q)
@inline convert{T}(::Type{RotMatrix{T}}, q::Quaternion) = convert(RotMatrix{T}, quat_to_rot(q))

@inline convert(::Type{Quaternion}, R::RotMatrix) = rot_to_quat(R)
@inline convert{T}(::Type{Quaternion{T}}, R::RotMatrix) = convert(Quaternion{T}, rot_to_quat(R))

# go from a Quaternion to any other representation via the Rotation Matrix
# @inline convert{T}(::Type{T}, q::Quaternion) = convert(T, quat_to_rot(q))  # go via the rotation representation
conversion_path[Quaternion] = RotMatrix
append!(defined_conversions, [(Quaternion, RotMatrix), (RotMatrix, Quaternion)])

# define its interaction with other angle representations
@inline eltype{T}(::Type{Quaternion{T}}) = T
@inline eltype{T}(::Quaternion{T}) = T

# Quaternions have an extra bool field that type_methods.jl doesn't want to know about
numel(::Type{Quaternion}) = 4
@inline getindex(X::Quaternion, i::Integer) = getfield(X, i)

# angle and axis functions
@inline function rotation_angle(q::Quaternion) # I think normalizing rounding errors will make things worse

    # this version to get the angle of the unit quaternion q_hat for arbitrary scaled q, q = s * q_hat
    theta =  2 * atan2(sqrt(q.x*q.x + q.y*q.y + q.z*q.z), q.w)

    #= If we want it to throw a domain error for non-unit quaternions
    if (abs(q.w) > 1)
        thresh = 1e-9  # choosen by voodoo
        if (q.w > 1)
            theta = (q.w - 1 < thresh) ?  2 * acos(one(q.w)) : 2 * acos(q.w)  # 2nd case will throw
        else
            theta = (q.w + 1 < thresh) ?  2 * acos(-one(q.w)) : 2 * acos(q.w) # 2nd case will throw
        end
    else
        theta = 2 * acos(q.w)
    end

end

@inline rotation_axis(q::Quaternion) = rotation_axis(AngleAxis(q))

strip_eltype{T <: Quaternion}(::Type{T}) = Quaternion

# define null rotations for convenience
@inline eye(::Type{Quaternion}) = Quaternion(1.0, 0.0, 0.0, 0.0)
@inline eye{T}(::Type{Quaternion{T}}) = Quaternion{T}(one(T), zero(T), zero(T), zero(T))



###################################################
# Stereographic projection of Quaternions
###################################################

"""

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

"""
immutable SpQuat{T}
    x::T
    y::T
    z::T
end

# add to the type list
push!(RotTypeList, SpQuat)

# enable mixed element construction...
@compat @inline (::Type{SpQuat})(x1, x2, x3) = SpQuat(promote(x1,x2,x3)...)

@inline getindex(X::SpQuat, i::Integer) = getfield(X, i)

# element conversions
@inline convert{T}(::Type{SpQuat{T}}, spq::SpQuat{T}) = spq
@inline convert{T}(::Type{SpQuat{T}}, spq::SpQuat) = SpQuat{T}(T(spq.x), T(spq.y), T(spq.z))

# convert to and from quaternions
@inline convert(::Type{Quaternion}, spq::SpQuat) = spquat_to_quat(spq)
@inline convert{T}(::Type{Quaternion{T}}, spq::SpQuat) = convert(Quaternion{T}, spquat_to_quat(spq))

@inline convert(::Type{SpQuat}, q::Quaternion) = quat_to_spquat(q)
@inline convert{T}(::Type{SpQuat{T}}, q::Quaternion) = convert(SpQuat{T}, quat_to_spquat(q))

# go from a SpQuat to any other representation via Quaternions
# @inline convert{T}(::Type{T}, spq::SpQuat) = convert(T, spquat_to_quat(spq))
conversion_path[SpQuat] = Quaternion
append!(defined_conversions, [(Quaternion, SpQuat), (SpQuat, Quaternion)])

# define an inverse for the SPQuat since its trivial.  Inverse in the sense of the corresponding inverse rotation
@inline inv(X::SpQuat) = SpQuat(-X.x, -X.y, -X.z)

# define multiplication for SpQuat corresponding to combining rotations
@inline *(lhs::SpQuat, rhs::SpQuat) = SpQuat(Quaternion(lhs) * Quaternion(rhs))

# rotation properties
@inline rotation_angle(spq::SpQuat) = rotation_angle(Quaternion(spq))
@inline rotation_axis(spq::SpQuat) = rotation_axis(Quaternion(spq))

# element type is handy
@inline eltype{T}(::Type{SpQuat{T}}) = T
@inline eltype{T}(::SpQuat{T}) = T
strip_eltype{T <: SpQuat}(::Type{T}) = SpQuat

# define null rotations for convenience
@inline eye(::Type{SpQuat}) = SpQuat(0.0, 0.0, 0.0)
@inline eye{T}(::Type{SpQuat{T}}) = SpQuat{T}(zero(T), zero(T), zero(T))



###################################################
# Arbitrary axis rotation representations
# (limited support for now)
###################################################

"""
Rotation by theta about an arbitrary axis [axis_x, axis_y, axis_z]
"""
immutable AngleAxis{T}
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T
end

# add to the type list
push!(RotTypeList, AngleAxis)

@inline getindex(X::AngleAxis, i::Integer) = getfield(X, i)

# enable mixed element construction...
@compat @inline (::Type{AngleAxis})(x1, x2, x3, x4) = AngleAxis(promote(x1,x2,x3,x4)...)

@inline convert{T}(::Type{AngleAxis{T}}, aa::AngleAxis{T}) = aa
@inline convert{T}(::Type{AngleAxis{T}}, aa::AngleAxis) = AngleAxis{T}(T(aa.theta), T(aa.axis_x), T(aa.axis_y), T(aa.axis_z))

# enable construction via an angle and a vector
@compat @inline (::Type{AngleAxis})(theta, axis::FixedVector{3}) = AngleAxis(theta, axis...)
@compat @inline function (::Type{AngleAxis})(theta, axis::AbstractVector)
    (length(axis) == 3) || throw(DimensionMismatch("The axis vector should have 3 elements (length(axis) == $(length(axis))"))
    AngleAxis(theta, axis[1], axis[2], axis[3])
end

# define its interaction with other angle representations
@inline convert(::Type{Quaternion}, aa::AngleAxis) = angleaxis_to_quat(aa)
@inline convert{T}(::Type{Quaternion{T}}, aa::AngleAxis) = convert(Quaternion{T}, angleaxis_to_quat(aa))

@inline convert(::Type{AngleAxis}, q::Quaternion) = quat_to_angleaxis(q)
@inline convert{T}(::Type{AngleAxis{T}}, q::Quaternion) = convert(AngleAxis{T}, quat_to_angleaxis(q))

@inline inv(aa::AngleAxis) = AngleAxis(-aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)

# go from an AngleAxis to any other representation via Quaternions
# @inline convert{T}(::Type{T}, aa::AngleAxis) = convert(T, angleaxis_to_quat(spq))
conversion_path[AngleAxis] = Quaternion
append!(defined_conversions, [(Quaternion, AngleAxis), (AngleAxis, Quaternion)])

# accessors
@inline rotation_angle(aa::AngleAxis) = aa.theta #  - floor((aa.theta+pi) / (2*pi)) * 2*pi
@inline rotation_axis(aa::AngleAxis) = Vec(aa.axis_x, aa.axis_y, aa.axis_z)

# element type is handy
@inline eltype{T}(::Type{AngleAxis{T}}) = T
@inline eltype{T}(::AngleAxis{T}) = T

strip_eltype{T <: AngleAxis}(::Type{T}) = AngleAxis


# define null rotations for convenience
@inline eye(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0.0, 0.0)
@inline eye{T}(::Type{AngleAxis{T}}) = AngleAxis{T}(zero(T), one(T), zero(T), zero(T))


###################################################
# Rodrigues Vector theta *axis
# Don't trust autodiff on these!
###################################################

"""
Rodrigues vector parameterization (stored as theta * axis)
"""
immutable RodriguesVec{T}
    sx::T
    sy::T
    sz::T
end

# add to the type list
push!(RotTypeList, RodriguesVec)

@inline getindex(X::RodriguesVec, i::Integer) = getfield(X, i)

# enable mixed element construction...
@compat @inline (::Type{RodriguesVec})(x1, x2, x3) = RodriguesVec(promote(x1,x2,x3)...)

@inline convert{T}(::Type{RodriguesVec{T}}, rv::RodriguesVec{T}) = rv
@inline convert{T}(::Type{RodriguesVec{T}}, rv::RodriguesVec) = RodriguesVec{T}(T(rv.sx), T(rv.sy), T(rv.sz))

# define its interaction with other angle representations
@inline convert(::Type{Quaternion}, rv::RodriguesVec) = rodrigues_to_quat(rv)
@inline convert{T}(::Type{Quaternion{T}}, rv::RodriguesVec) = convert(Quaternion{T}, rodrigues_to_quat(rv))
@inline convert(::Type{AngleAxis}, rv::RodriguesVec) = rodrigues_to_angleaxis(rv)
@inline convert{T}(::Type{AngleAxis{T}}, rv::RodriguesVec) = convert(AngleAxis{T}, rodrigues_to_angleaxis(rv))

@inline convert(::Type{RodriguesVec}, q::Quaternion) = quat_to_rodrigues(q)
@inline convert{T}(::Type{RodriguesVec{T}}, q::Quaternion) = convert(RodriguesVec{T}, quat_to_rodrigues(q))
@inline convert(::Type{RodriguesVec}, aa::AngleAxis) = angleaxis_to_rodrigues(aa)
@inline convert{T}(::Type{RodriguesVec{T}}, aa::AngleAxis) = convert(RodriguesVec{T}, angleaxis_to_rodrigues(aa))

# get here from a rotation matrix
@inline inv(rv::RodriguesVec) = RodriguesVec(-rv.sx, -rv.sy, -rv.sz)

# go from an RodriguesVec to any other representation via Quaternions
conversion_path[RodriguesVec] = Quaternion
append!(defined_conversions, [(Quaternion, RodriguesVec), (RodriguesVec, Quaternion)])
append!(defined_conversions, [(AngleAxis, RodriguesVec), (RodriguesVec, AngleAxis)])

# rotation properties
@inline norm(rv::RodriguesVec) = sqrt(rv.sx * rv.sx + rv.sy * rv.sy + rv.sz * rv.sz) # norm is meaningful for these things
@inline rotation_angle(rv::RodriguesVec) = norm(rv)
function rotation_axis(rv::RodriguesVec)     # what should this return for theta = 0?
    theta = norm(rv)
    return (theta > 0 ? Vec(rv.sx / theta, rv.sy / theta, rv.sz / theta) : Vec(one(theta), zero(theta), zero(theta)))
end

# element type is handy
@inline eltype{T}(::Type{RodriguesVec{T}}) = T
@inline eltype{T}(::RodriguesVec{T}) = T

strip_eltype{T <: RodriguesVec}(::Type{T}) = RodriguesVec

# define null rotations for convenience
@inline eye(::Type{RodriguesVec}) = RodriguesVec(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RodriguesVec{T}}) = RodriguesVec{T}(zero(T), zero(T), zero(T))


###################################################
# Euler Angles
###################################################

# this file defines the types:
# EulerAngles       - for tait bryan ordering
# ProperEulerAngles - for proper Euler ordering
include("euler_types.jl")

# add to the type list
append!(RotTypeList, [EulerAngles, ProperEulerAngles])

#
# Euler Angles
#
@inline convert(::Type{EulerAngles}, x, y, z) = EulerAngles{DefaultEulerOrder()}(x,y,z)

@inline getindex(X::EulerAngles, i::Integer) = getfield(X, i)

@inline convert(::Type{EulerAngles}, ea::EulerAngles) = ea
@inline convert{ORD <: TaitByranOrder,T}(::Type{EulerAngles{ORD}}, ea::EulerAngles{ORD,T}) = ea
@inline convert{ORD <: TaitByranOrder,T}(::Type{EulerAngles{ORD,T}}, ea::EulerAngles{ORD,T}) = ea
@inline convert{ORD <: TaitByranOrder,T, U}(::Type{EulerAngles{ORD,T}}, ea::EulerAngles{ORD,U}) = EulerAngles{ORD,T}(T(ea.theta_x), T(ea.theta_y), T(ea.theta_z))


# define its interaction with RotMatrixs
@inline convert(::Type{RotMatrix}, ea::EulerAngles) = euler_to_rot(ea)
@inline convert{T}(::Type{RotMatrix{T}}, ea::EulerAngles) = convert(RotMatrix{T}, euler_to_rot(ea))

@inline convert(::Type{EulerAngles}, R::RotMatrix) = rot_to_euler(EulerAngles{DefaultEulerOrder()}, R)
@inline convert{ORD}(::Type{EulerAngles{ORD}}, R::RotMatrix) = rot_to_euler(EulerAngles{ORD}, R)
@inline convert{ORD, T}(::Type{EulerAngles{ORD,T}}, R::RotMatrix) = rot_to_euler(EulerAngles{ORD,T}, R)

# go from a EulerAngles to any other representation via a rotation matrix
# @inline convert{T}(::Type{T}, ea::EulerAngles) = convert(T, euler_to_rot(ea))
conversion_path[EulerAngles] = RotMatrix
append!(defined_conversions, [(RotMatrix, EulerAngles), (EulerAngles, RotMatrix)])

# element type is handy
@inline eltype{ORD,T}(::Type{EulerAngles{ORD,T}}) = T
@inline eltype{ORD,T}(::EulerAngles{ORD,T}) = T

# order as well
@inline euler_order(::Type{EulerAngles}) = DefaultEulerOrder()
@inline euler_order{ORD}(::Type{EulerAngles{ORD}}) = ORD
@inline euler_order{ORD,T}(::Type{EulerAngles{ORD,T}}) = ORD
@inline euler_order{ORD,T}(::EulerAngles{ORD,T}) = ORD

# need special constructors to fill in the template parametera
@compat @inline (::Type{EulerAngles})(x, y, z) = EulerAngles{DefaultEulerOrder()}(promote(x, y, z)...)
@compat @inline (::Type{EulerAngles{ORD}}){ORD, T}(x::T, y::T, z::T) = EulerAngles{ORD,T}(x, y, z)
@compat @inline (::Type{EulerAngles{ORD}}){ORD}(x, y, z) = EulerAngles{ORD}(promote(x, y, z)...)

# convert from one order to another
@inline convert{ORD1, ORD2, eT}(::Type{EulerAngles{ORD1}}, ea::EulerAngles{ORD2, eT}) = convert(EulerAngles{ORD1}, euler_to_rot(ea))
@inline convert{ORD1, ORD2, eT1, eT2}(::Type{EulerAngles{ORD1, eT1}}, ea::EulerAngles{ORD2, eT2}) = convert(EulerAngles{ORD1, eT1}, euler_to_rot(ea))

strip_eltype{T <: EulerAngles}(::Type{T}) = EulerAngles{euler_order(T)}


#
# Proper Euler Angles
#
@inline convert(::Type{ProperEulerAngles}, x, y, z) = EulerAngles{DefaultProperEulerOrder()}(x,y,z)

@inline getindex(X::ProperEulerAngles, i::Integer) = getfield(X, i)

@inline convert(::Type{ProperEulerAngles}, ea::ProperEulerAngles) = ea
@inline convert{ORD <: ProperEulerOrder, T}(::Type{ProperEulerAngles{ORD}}, ea::ProperEulerAngles{ORD,T}) = ea
@inline convert{ORD <: ProperEulerOrder, T}(::Type{ProperEulerAngles{ORD,T}}, ea::ProperEulerAngles{ORD,T}) = ea
@inline convert{ORD <: ProperEulerOrder, T, U}(::Type{ProperEulerAngles{ORD,T}}, ea::ProperEulerAngles{ORD,U}) = ProperEulerAngles{ORD,T}(T(ea.theta_1), T(ea.theta_2), T(ea.theta_3))

# define its interaction with RotMatrixs
@inline convert(::Type{RotMatrix}, ea::ProperEulerAngles) = euler_to_rot(ea)
@inline convert{T}(::Type{RotMatrix{T}}, ea::ProperEulerAngles) = convert(RotMatrix{T}, euler_to_rot(ea))

@inline convert(::Type{ProperEulerAngles}, R::RotMatrix) = rot_to_euler(ProperEulerAngles{DefaultProperEulerOrder()}, R)
@inline convert{ORD}(::Type{ProperEulerAngles{ORD}}, R::RotMatrix) = rot_to_euler(ProperEulerAngles{ORD}, R)
@inline convert{ORD,T}(::Type{ProperEulerAngles{ORD,T}}, R::RotMatrix) = rot_to_euler(ProperEulerAngles{ORD,T}, R)

# go from a ProperEulerAngles to any other representation via a rotation matrix
# @inline convert{T}(::Type{T}, ea::ProperEulerAngles) = convert(T, euler_to_rot(ea))
conversion_path[ProperEulerAngles] = RotMatrix
append!(defined_conversions, [(RotMatrix, ProperEulerAngles), (ProperEulerAngles, RotMatrix)])

# element type is handy
@inline eltype{ORD,T}(::Type{ProperEulerAngles{ORD,T}}) = T
@inline eltype{ORD,T}(::ProperEulerAngles{ORD,T}) = T

# order as well
@inline euler_order(::Type{ProperEulerAngles}) = DefaultProperEulerOrder()
@inline euler_order{ORD}(::Type{ProperEulerAngles{ORD}}) = ORD
@inline euler_order{ORD,T}(::Type{ProperEulerAngles{ORD,T}}) = ORD
@inline euler_order{ORD,T}(::ProperEulerAngles{ORD,T}) = ORD

# need special constructors to fill in the template parametera
@compat @inline (::Type{ProperEulerAngles})(x, y, z) = ProperEulerAngles{DefaultProperEulerOrder()}(promote(x, y, z)...)
@compat @inline (::Type{ProperEulerAngles{ORD}}){ORD, T}(x::T, y::T, z::T) = ProperEulerAngles{ORD,T}(x, y, z)
@compat @inline (::Type{ProperEulerAngles{ORD}}){ORD}(x, y, z) = ProperEulerAngles{ORD}(promote(x, y, z)...)

# convert from one order to another
@inline convert{ORD1, ORD2, eT}(::Type{ProperEulerAngles{ORD1}}, ea::ProperEulerAngles{ORD2, eT}) = convert(ProperEulerAngles{ORD1}, euler_to_rot(ea))
@inline convert{ORD1, ORD2, eT1, eT2}(::Type{ProperEulerAngles{ORD1, eT1}}, ea::ProperEulerAngles{ORD2, eT2}) = convert(ProperEulerAngles{ORD1, eT1}, euler_to_rot(ea))

strip_eltype{T <: ProperEulerAngles}(::Type{T}) = ProperEulerAngles{euler_order(T)}

# define null rotations for convenience
@inline eye(::Type{EulerAngles}) = EulerAngles(0.0, 0.0, 0.0)
@inline eye{ORD}(::Type{EulerAngles{ORD}}) = EulerAngles{ORD}(0.0, 0.0, 0.0)
@inline eye{ORD, T}(::Type{EulerAngles{ORD, T}}) = EulerAngles{ORD, T}(zero(T), zero(T), zero(T))

# define null rotations for convenience
@inline eye(::Type{ProperEulerAngles}) = ProperEulerAngles(0.0, 0.0, 0.0)
@inline eye{ORD}(::Type{ProperEulerAngles{ORD}}) = ProperEulerAngles{ORD}(0.0, 0.0, 0.0)
@inline eye{ORD, T}(::Type{ProperEulerAngles{ORD, T}}) = ProperEulerAngles{ORD, T}(zero(T), zero(T), zero(T))
=#
=#
