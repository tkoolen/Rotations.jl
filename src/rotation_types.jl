"""
    abstract Rotation{T} <: StaticMatrix{T}

An abstract type representing 3D rotations. More abstractly, they represent
3×3 unitary (orthogonal) matrices.
"""
abstract Rotation{T} <: StaticMatrix{T}

Base.size{R <: Rotation}(::Union{R,Type{R}}) = (3,3)
Base.ctranspose(r::Rotation) = inv(r)
Base.transpose{T<:Real}(r::Rotation{T}) = inv(r)

rotation_angle(r::Rotation) = rotation_angle(AngleAxis(r))
rotation_axis(r::Rotation) = rotation_axis(AngleAxis(r))

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
    "3×3 RotMatrix{$(eltype(r))}"
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
function (::Type{Q}){Q<:Quat}(t::NTuple{9})
    q = Q(sqrt(abs(1  + t[1] + t[5] + t[9])) / 2,
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

function Base.:*(q::Quat, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    qo = (-q.x * v[1] - q.y * v[2] - q.z * v[3],
           q.w * v[1] + q.y * v[3] - q.z * v[2],
           q.w * v[2] - q.x * v[3] + q.z * v[1],
           q.w * v[3] + q.x * v[2] - q.y * v[1])

    T = promote_type(eltype(q), eltype(v))

    return similar_type(v, T)(-qo[1] * q.x + qo[2] * q.w - qo[3] * q.z + qo[4] * q.y,
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

An `SPQuat` is a 3D rotation matrix represented by the "stereographic projection" of a normalized quaternion (shortened to "SPQuat"), which is
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
@inline Base.:*(spq1::SPQuat, spq2::SPQuat) = Quat(spq1) * Quat(spq2)

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


################################################################################
################################################################################
"""
    immutable AngleAxis{T} <: Rotation{T}
    AngleAxis(Θ, x, y, z)

A 3×3 rotation matrix parameterized by a 3D rotation by angle θ about an
arbitrary axis `[x, y, z]`.

Note that the axis is not unique for θ = 0, and that this parameterization does
not continuously map the neighbourhood of the null rotation (and therefore
might not be suitable for autodifferentation and optimization purposes).
"""
immutable AngleAxis{T} <: Rotation{T}
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T

    # Ensure axis is normalized
    function AngleAxis(θ, x, y, z)
        norm = sqrt(x*x + y*y + z*z)
        # Not sure what to do with theta?? Should it become theta * norm ?
        new(θ, x/norm, y/norm, z/norm)
    end
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline (::Type{AngleAxis}){Θ,X,Y,Z}(θ::Θ, x::X, y::Y, z::Z) = AngleAxis{promote_type(promote_type(promote_type(Θ, X), Y), Z)}(θ, x, y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{AngleAxis})(t::NTuple{9}) = AngleAxis(Quat(t))
@inline Base.getindex(aa::AngleAxis, i::Integer) = Quat(aa)[i]

@inline function Base.convert{Q <: Quat}(::Type{Q}, aa::AngleAxis)
    qtheta = cos(aa.theta / 2)
    s = sin(aa.theta / 2) / sqrt(aa.axis_x * aa.axis_x + aa.axis_y * aa.axis_y + aa.axis_z * aa.axis_z)
    return Q(qtheta, s * aa.axis_x, s * aa.axis_y, s * aa.axis_z)
end

@inline function Base.convert{AA <: AngleAxis}(::Type{AA}, q::Quat)
    # TODO: consider how to deal with derivative near theta = 0
    s = sqrt(q.x*q.x + q.y*q.y + q.z*q.z)
    theta =  2 * atan2(s, q.w)
    return s > 0 ? AA(theta, q.x / s, q.y / s, q.z / s) : AA(theta, one(theta), zero(theta), zero(theta))
end

@inline Base.convert(::Type{Tuple}, aa::AngleAxis) = Tuple(Quat(aa))

# Using Rodrigues formula on an AngleAxis parameterization (assume unit axis length) to do the rotation
# (implementation from: https://ceres-solver.googlesource.com/ceres-solver/+/1.10.0/include/ceres/rotation.h)
function Base.:*(aa::AngleAxis, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    w = rotation_axis(aa)
    ct, st = cos(aa.theta), sin(aa.theta)
    w_cross_pt = cross(w, v)
    m = dot(v, w) * (one(w_cross_pt[1]) - ct)
    T = promote_type(eltype(aa), eltype(v))
    return similar_type(v,T)(v[1] * ct + w_cross_pt[1] * st + w[1] * m,
                             v[2] * ct + w_cross_pt[2] * st + w[2] * m,
                             v[3] * ct + w_cross_pt[3] * st + w[3] * m)
end

@inline Base.:*(aa::AngleAxis, r::Rotation) = Quat(aa) * r
@inline Base.:*(r::Rotation, aa::AngleAxis) = r * Quat(aa)
@inline Base.:*(aa1::AngleAxis, aa2::AngleAxis) = Quat(aa1) * Quat(aa2)

@inline inv(aa::AngleAxis) = AngleAxis(-aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)

# define null rotations for convenience
@inline eye(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0.0, 0.0)
@inline eye{T}(::Type{AngleAxis{T}}) = AngleAxis{T}(zero(T), one(T), zero(T), zero(T))

# accessors
@inline rotation_angle(aa::AngleAxis) = aa.theta #  - floor((aa.theta+pi) / (2*pi)) * 2*pi
@inline rotation_axis(aa::AngleAxis) = SVector(aa.axis_x, aa.axis_y, aa.axis_z)

# Need SPQuat parameters in the output representation
function Base.summary(aa::AngleAxis)
    "3×3 AngleAxis{$(eltype(aa))}($(aa.theta) rad, [$(aa.axis_x), $(aa.axis_y), $(aa.axis_z)])"
end


################################################################################
################################################################################
"""
    immutable RodriguesVec{T} <: Rotation{T}
    RodriguesVec(sx, sy, sz)

Rodrigues vector parameterization of a 3×3 rotation matrix. The direction of the
vector [sx, sy, sz] defines the axis of rotation, and the rotation angle is
given by its norm.
"""
immutable RodriguesVec{T} <: Rotation{T}
    sx::T
    sy::T
    sz::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline (::Type{RodriguesVec}){X,Y,Z}(x::X, y::Y, z::Z) = RodriguesVec{promote_type(promote_type(X, Y), Z)}(x, y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RodriguesVec})(t::NTuple{9}) = RodriguesVec(Quat(t))
@inline Base.getindex(aa::RodriguesVec, i::Integer) = Quat(aa)[i]

# define its interaction with other angle representations

function Base.convert{AA <: AngleAxis}(::Type{AA}, rv::RodriguesVec)
    # TODO: consider how to deal with derivative near theta = 0. There should be a first-order expansion here.
    theta = rotation_angle(rv)
    return theta > 0 ? AA(theta, rv.sx / theta, rv.sy / theta, rv.sz / theta) : AA(zero(theta), one(theta), zero(theta), zero(theta))
end

function Base.convert{RV <: RodriguesVec}(::Type{RV}, aa::AngleAxis)
    return RV(aa.theta * aa.axis_x, aa.theta * aa.axis_y, aa.theta * aa.axis_z)
end

function Base.convert{Q <: Quat}(::Type{Q}, rv::RodriguesVec)
    theta = rotation_angle(rv)
    qtheta = cos(theta / 2)
    #s = abs(1/2 * sinc((theta / 2) / pi))
    s = (1/2 * sinc((theta / 2) / pi)) # TODO check this (I removed an abs)
    return Q(qtheta, s * rv.sx, s * rv.sy, s * rv.sz)
end

function Base.convert{RV <: RodriguesVec}(::Type{RV}, q::Quat)
    s2 = q.x*q.x + q.y*q.y + q.z*q.z
    if (s2 > 0)
        cos_t2 = sqrt(s2)
        theta = 2 * atan2(cos_t2, q.w)
        sc = theta / cos_t2
    else
        sc = 2                 # N.B. the 2 "should" match the derivitive as cos_t2 -> 0
    end
    return RV(sc * q.x, sc * q.y, sc * q.z )
end

@inline Base.convert(::Type{Tuple}, rv::RodriguesVec) = Tuple(Quat(rv))

function Base.:*{T1,T2}(rv::RodriguesVec{T1}, v::StaticVector{T2})
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    theta = rotation_angle(rv)
    if (theta > eps(T1)) # use eps here because we have the 1st order series expansion defined
        return AngleAxis(rv) * v
    else
        return similar_type(typeof(v), promote_type(T1,T2))(
                    v[1] + rv[2] * v[3] - rv[3] * v[2],
                    v[2] + rv[3] * v[1] - rv[1] * v[3],
                    v[3] + rv[1] * v[2] - rv[2] * v[1])
    end
end

@inline Base.:*(rv::RodriguesVec, r::Rotation) = Quat(rv) * r
@inline Base.:*(r::Rotation, rv::RodriguesVec) = r * Quat(rv)
@inline Base.:*(rv1::RodriguesVec, rv2::RodriguesVec) = Quat(rv1) * Quat(rv2)

@inline inv(rv::RodriguesVec) = RodriguesVec(-rv.sx, -rv.sy, -rv.sz)

# rotation properties
@inline rotation_angle(rv::RodriguesVec) = sqrt(rv.sx * rv.sx + rv.sy * rv.sy + rv.sz * rv.sz)
function rotation_axis(rv::RodriguesVec)     # what should this return for theta = 0?
    theta = rotation_angle(rv)
    return (theta > 0 ? SVector(rv.sx / theta, rv.sy / theta, rv.sz / theta) : SVector(one(theta), zero(theta), zero(theta)))
end

# define null rotations for convenience
@inline eye(::Type{RodriguesVec}) = RodriguesVec(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RodriguesVec{T}}) = RodriguesVec{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(rv::RodriguesVec)
    "3×3 RodriguesVec{$(eltype(rv))}($(rv.sx), $(rv.sy), $(rv.sz))"
end

################################################################################
################################################################################

#########################
# Single axis rotations #
#########################
"""
    immutable RotX{T} <: Rotation{T}
    RotX(theta)

A 3×3 rotation matrix which represents a rotation by `theta` about the X axis.
"""
immutable RotX{T} <: Rotation{T}
    theta::T
end

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RotX})(t::NTuple{9}) = error("Cannot construct a cardinal axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotX{T}, i::Integer)
    if i == 1
        one(T)
    elseif i < 5
        zero(T)
    elseif i == 5
        cos(r.theta)
    elseif i == 6
        sin(r.theta)
    elseif i == 7
        zero(T)
    elseif i == 8
        -sin(r.theta)
    elseif i == 9
        cos(r.theta)
    else
        throw(BoundsError(r,i))
    end
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotX{T})
    (one(T),  zero(T),      zero(T),   # transposed representation
     zero(T), cos(r.theta), sin(r.theta),
     zero(T), -sin(r.theta), cos(r.theta))
end

@inline function Base.:*(r::RotX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    ct, st = cos(r.theta), sin(r.theta)
    T = promote_type(eltype(r), eltype(v))
    return similar_type(v,T)(v[1],
                             v[2] * ct - v[3] * st,
                             v[3] * ct + v[2] * st)
end

@inline Base.:*(r1::RotX, r2::RotX) = RotX(r1.theta + r2.theta)

@inline inv(r::RotX) = RotX(-r.theta)

# define null rotations for convenience
@inline eye(::Type{RotX}) = RotX(0.0)
@inline eye{T}(::Type{RotX{T}}) = RotX{T}(zero(T))

# Need a good output representation
function Base.summary(r::RotX)
    "3×3 RotX{$(eltype(r))}($(r.theta))"
end

"""
    immutable RotY{T} <: Rotation{T}
    RotY(theta)

A 3×3 rotation matrix which represents a rotation by `theta` about the Y axis.
"""
immutable RotY{T} <: Rotation{T}
    theta::T
end

@inline (::Type{RotY})(t::NTuple{9}) = error("Cannot construct a cardinal axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotY{T}, i::Integer)
    if i == 1
        cos(r.theta)
    elseif i == 2
        zero(T)
    elseif i == 3
        -sin(r.theta)
    elseif i == 4
        zero(T)
    elseif i == 5
        one(T)
    elseif i == 6
        zero(T)
    elseif i == 7
        sin(r.theta)
    elseif i == 8
        zero(T)
    elseif i == 9
        cos(r.theta)
    else
        throw(BoundsError(r,i))
    end
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotY{T})
    (cos(r.theta), zero(T), -sin(r.theta),   # transposed representation
     zero(T),      one(T),   zero(T),
     sin(r.theta), zero(T),  cos(r.theta))
end

@inline function Base.:*(r::RotY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    ct, st = cos(r.theta), sin(r.theta)
    T = promote_type(eltype(r), eltype(v))
    return similar_type(v,T)(v[1] * ct + v[3] * st,
                             v[2],
                             v[3] * ct - v[1] * st)
end

@inline Base.:*(r1::RotY, r2::RotY) = RotY(r1.theta + r2.theta)

@inline inv(r::RotY) = RotY(-r.theta)

# define null rotations for convenience
@inline eye(::Type{RotY}) = RotY(0.0)
@inline eye{T}(::Type{RotY{T}}) = RotY{T}(zero(T))

# Need a good output representation
function Base.summary(r::RotY)
    "3×3 RotY{$(eltype(r))}($(r.theta))"
end

"""
    immutable RotZ{T} <: Rotation{T}
    RotZ(theta)


A 3×3 rotation matrix which represents a rotation by `theta` about the Z axis.
"""
immutable RotZ{T} <: Rotation{T}
    theta::T
end

@inline (::Type{RotZ})(t::NTuple{9}) = error("Cannot construct a cardinal axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotZ{T}, i::Integer)
    if i == 1
        cos(r.theta)
    elseif i == 2
        sin(r.theta)
    elseif i == 3
        zero(T)
    elseif i == 4
        -sin(r.theta)
    elseif i == 5
        cos(r.theta)
    elseif i < 9
        zero(T)
    elseif i == 9
        one(T)
    else
        throw(BoundsError(r,i))
    end
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotZ{T})
    ( cos(r.theta), sin(r.theta), zero(T),   # transposed representation
     -sin(r.theta), cos(r.theta), zero(T),
      zero(T),      zero(T),      one(T))
end

@inline function Base.:*(r::RotZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    ct, st = cos(r.theta), sin(r.theta)
    T = promote_type(eltype(r), eltype(v))
    return similar_type(v,T)(v[1] * ct - v[2] * st,
                             v[2] * ct + v[1] * st,
                             v[3])
end

@inline Base.:*(r1::RotZ, r2::RotZ) = RotZ(r1.theta + r2.theta)

# define null rotations for convenience
@inline eye(::Type{RotZ}) = RotZ(0.0)
@inline eye{T}(::Type{RotZ{T}}) = RotZ{T}(zero(T))

@inline inv(r::RotZ) = RotZ(-r.theta)

# Need a good output representation
function Base.summary(r::RotZ)
    "3×3 RotZ{$(eltype(r))}($(r.theta))"
end

################################################################################
################################################################################

######################
# Two axis rotations #
######################

"""
    immutable RotXY{T} <: Rotation{T}
    RotXY(theta_x, theta_y)

A 3×3 rotation matrix which represents a rotation by `theta_y` about the Y axis,
followed by a rotation by `theta_x` about the X axis.
"""
immutable RotXY{T} <: Rotation{T}
    theta1::T
    theta2::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 2 inputs (not 9).
@inline (::Type{RotXY}){X,Y}(x::X, y::Y) = RotXY{promote_type(X, Y)}(x, y)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RotXY})(t::NTuple{9}) = error("Cannot construct a two-axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotXY{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotXY{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    # transposed representation
    (cosθ₂,  sinθ₁*sinθ₂,  cosθ₁*-sinθ₂,
     zero(T),      cosθ₁,          sinθ₁,
     sinθ₂,  -sinθ₁*cosθ₂,   cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotXY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    return similar_type(v,T)(cosθ₂*v[1] + sinθ₂*v[3],
                             sinθ₁*sinθ₂*v[1] + cosθ₁*v[2] + -sinθ₁*cosθ₂*v[3],
                             cosθ₁*-sinθ₂*v[1] + sinθ₁*v[2] + cosθ₁*cosθ₂*v[3])
end

@inline Base.:*(r1::RotX, r2::RotY) = RotXY(r1.theta, r2.theta)
@inline Base.:*(r1::RotXY, r2::RotY) = RotXY(r1.theta1, r1.theta2 + r2.theta)
@inline Base.:*(r1::RotX, r2::RotXY) = RotXY(r1.theta + r2.theta1, r2.theta2)

@inline inv(r::RotXY) = RotYX(-r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotXY}) = RotXY(0.0, 0.0)
@inline eye{T}(::Type{RotXY{T}}) = RotXY{T}(zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotXY)
    "3×3 RotXY{$(eltype(r))}($(r.theta1), $(r.theta2))"
end


"""
    immutable RotYX{T} <: Rotation{T}
    RotYX(theta_y, theta_x)

A 3×3 rotation matrix which represents a rotation by `theta_x` about the X axis,
followed by a rotation by `theta_y` about the Y axis.
"""
immutable RotYX{T} <: Rotation{T}
    theta1::T
    theta2::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 2 inputs (not 9).
@inline (::Type{RotYX}){Y,X}(y::Y, x::X) = RotYX{promote_type(Y, X)}(y, x)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RotYX})(t::NTuple{9}) = error("Cannot construct a two-axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotYX{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotYX{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    # transposed representation
    (cosθ₁,        zero(T),       -sinθ₁,
     sinθ₁*sinθ₂,  cosθ₂,   cosθ₁*sinθ₂,
     sinθ₁*cosθ₂,  -sinθ₂,  cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotYX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    return similar_type(v,T)(cosθ₁*v[1] + sinθ₁*sinθ₂*v[2] + sinθ₁*cosθ₂*v[3],
                             cosθ₂*v[2] + -sinθ₂*v[3],
                             -sinθ₁*v[1] + cosθ₁*sinθ₂*v[2] + cosθ₁*cosθ₂*v[3])
end

@inline Base.:*(r1::RotY, r2::RotX) = RotYX(r1.theta, r2.theta)
@inline Base.:*(r1::RotYX, r2::RotX) = RotYX(r1.theta1, r1.theta2 + r2.theta)
@inline Base.:*(r1::RotY, r2::RotYX) = RotYX(r1.theta + r2.theta1, r2.theta2)

@inline inv(r::RotYX) = RotXY(-r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotYX}) = RotYX(0.0, 0.0)
@inline eye{T}(::Type{RotYX{T}}) = RotYX{T}(zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotYX)
    "3×3 RotYX{$(eltype(r))}($(r.theta1), $(r.theta2))"
end


"""
    immutable RotXZ{T} <: Rotation{T}
    RotXZ(theta_x, theta_z)

A 3×3 rotation matrix which represents a rotation by `theta_z` about the Z axis,
followed by a rotation by `theta_x` about the X axis.
"""
immutable RotXZ{T} <: Rotation{T}
    theta1::T
    theta2::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 2 inputs (not 9).
@inline (::Type{RotXZ}){X,Z}(x::X, z::Z) = RotXZ{promote_type(X, Z)}(x, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RotXZ})(t::NTuple{9}) = error("Cannot construct a two-axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotXZ{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotXZ{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    # transposed representation
    (cosθ₂,   cosθ₁*sinθ₂,  sinθ₁*sinθ₂,
     -sinθ₂,  cosθ₁*cosθ₂,  sinθ₁*cosθ₂,
     zero(T),       -sinθ₁,       cosθ₁)
end

@inline function Base.:*(r::RotXZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    return similar_type(v,T)(cosθ₂*v[1] + -sinθ₂*v[2],
                             cosθ₁*sinθ₂*v[1] + cosθ₁*cosθ₂*v[2] + -sinθ₁*v[3],
                             sinθ₁*sinθ₂*v[1] + sinθ₁*cosθ₂*v[2] + cosθ₁*v[3])
end

@inline Base.:*(r1::RotX, r2::RotZ) = RotXZ(r1.theta, r2.theta)
@inline Base.:*(r1::RotXZ, r2::RotZ) = RotXZ(r1.theta1, r1.theta2 + r2.theta)
@inline Base.:*(r1::RotX, r2::RotXZ) = RotXZ(r1.theta + r2.theta1, r2.theta2)

@inline inv(r::RotXZ) = RotZX(-r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotXZ}) = RotXZ(0.0, 0.0)
@inline eye{T}(::Type{RotXZ{T}}) = RotXZ{T}(zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotXZ)
    "3×3 RotXZ{$(eltype(r))}($(r.theta1), $(r.theta2))"
end


"""
    immutable RotZX{T} <: Rotation{T}
    RotZX(theta_z, theta_x)

A 3×3 rotation matrix which represents a rotation by `theta_x` about the X axis,
followed by a rotation by `theta_z` about the Z axis.
"""
immutable RotZX{T} <: Rotation{T}
    theta1::T
    theta2::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 2 inputs (not 9).
@inline (::Type{RotZX}){Z,X}(z::Z, x::X) = RotZX{promote_type(Z, X)}(z, x)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RotZX})(t::NTuple{9}) = error("Cannot construct a two-axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotZX{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotZX{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    # transposed representation
    (cosθ₂,   cosθ₁*sinθ₂,  sinθ₁*sinθ₂,
     -sinθ₂,  cosθ₁*cosθ₂,  sinθ₁*cosθ₂,
     zero(T),       -sinθ₁,       cosθ₁)
end

@inline function Base.:*(r::RotZX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    return similar_type(v,T)(cosθ₂*v[1] + -sinθ₂*v[2],
                             cosθ₁*sinθ₂*v[1] + cosθ₁*cosθ₂*v[2] + -sinθ₁*v[3],
                             sinθ₁*sinθ₂*v[1] + sinθ₁*cosθ₂*v[2] + cosθ₁*v[3])
end

@inline Base.:*(r1::RotZ, r2::RotX) = RotZX(r1.theta, r2.theta)
@inline Base.:*(r1::RotZX, r2::RotX) = RotZX(r1.theta1, r1.theta2 + r2.theta)
@inline Base.:*(r1::RotZ, r2::RotZX) = RotZX(r1.theta + r2.theta1, r2.theta2)

@inline inv(r::RotZX) = RotXZ(-r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotZX}) = RotZX(0.0, 0.0)
@inline eye{T}(::Type{RotZX{T}}) = RotZX{T}(zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotZX)
    "3×3 RotZX{$(eltype(r))}($(r.theta1), $(r.theta2))"
end


"""
    immutable RotZY{T} <: Rotation{T}
    RotZY(theta_z, theta_y)

A 3×3 rotation matrix which represents a rotation by `theta_y` about the Y axis,
followed by a rotation by `theta_z` about the Z axis.
"""
immutable RotZY{T} <: Rotation{T}
    theta1::T
    theta2::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 2 inputs (not 9).
@inline (::Type{RotZY}){Z,Y}(z::Z, y::Y) = RotZY{promote_type(Z, Y)}(z, y)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RotZY})(t::NTuple{9}) = error("Cannot construct a two-axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotZY{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotZY{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    # transposed representation
    ( cosθ₁*cosθ₂,  sinθ₁*cosθ₂, -sinθ₂,
     -sinθ₁,        cosθ₁,        zero(T),
      cosθ₁*sinθ₂,  sinθ₁*sinθ₂,  cosθ₂)
end

@inline function Base.:*(r::RotZY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    return similar_type(v,T)(cosθ₁*cosθ₂*v[1] + -sinθ₁*v[2] + cosθ₁*sinθ₂*v[3],
                             sinθ₁*cosθ₂*v[1] + cosθ₁*v[2] + sinθ₁*sinθ₂*v[3],
                             -sinθ₂*v[1] + cosθ₂*v[3])
end

@inline Base.:*(r1::RotZ, r2::RotY) = RotZY(r1.theta, r2.theta)
@inline Base.:*(r1::RotZY, r2::RotY) = RotZY(r1.theta1, r1.theta2 + r2.theta)
@inline Base.:*(r1::RotZ, r2::RotZY) = RotZY(r1.theta + r2.theta1, r2.theta2)

@inline inv(r::RotZY) = RotYZ(-r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotZY}) = RotZY(0.0, 0.0)
@inline eye{T}(::Type{RotZY{T}}) = RotZY{T}(zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotZY)
    "3×3 RotZY{$(eltype(r))}($(r.theta1), $(r.theta2))"
end


"""
    immutable RotYZ{T} <: Rotation{T}
    RotYZ(theta_y, theta_z)

A 3×3 rotation matrix which represents a rotation by `theta_z` about the Z axis,
followed by a rotation by `theta_y` about the Y axis.
"""
immutable RotYZ{T} <: Rotation{T}
    theta1::T
    theta2::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 2 inputs (not 9).
@inline (::Type{RotYZ}){Y,Z}(y::Y, z::Z) = RotYZ{promote_type(Y, Z)}(y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RotYZ})(t::NTuple{9}) = error("Cannot construct a two-axis rotation from a matrix")
@inline function Base.getindex{T}(r::RotYZ{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotYZ{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    # transposed representation
    (cosθ₁*cosθ₂,   sinθ₂,    -sinθ₁*cosθ₂,
     cosθ₁*-sinθ₂,  cosθ₂,     sinθ₁*sinθ₂,
     sinθ₁,         zero(T),      cosθ₁)
end

@inline function Base.:*(r::RotYZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)

    return similar_type(v,T)(cosθ₁*cosθ₂*v[1] + cosθ₁*-sinθ₂*v[2] + sinθ₁*v[3],
                             sinθ₂*v[1] + cosθ₂*v[2],
                             -sinθ₁*cosθ₂*v[1] + sinθ₁*sinθ₂*v[2] + cosθ₁*v[3])
end

@inline Base.:*(r1::RotY, r2::RotZ) = RotYZ(r1.theta, r2.theta)
@inline Base.:*(r1::RotYZ, r2::RotZ) = RotYZ(r1.theta1, r1.theta2 + r2.theta)
@inline Base.:*(r1::RotY, r2::RotYZ) = RotYZ(r1.theta + r2.theta1, r2.theta2)

# define null rotations for convenience
@inline eye(::Type{RotYZ}) = RotYZ(0.0, 0.0)
@inline eye{T}(::Type{RotYZ{T}}) = RotYZ{T}(zero(T), zero(T))

@inline inv(r::RotYZ) = RotZY(-r.theta2, -r.theta1)

# Need a good output representation
function Base.summary(r::RotYZ)
    "3×3 RotYZ{$(eltype(r))}($(r.theta1), $(r.theta2))"
end

################################################################################
################################################################################

##########################
# Proper Euler Rotations #
##########################

"""
    immutable RotXYX{T} <: Rotation{T}
    RotXYX(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" XYX Euler angle convention,
consisting of first a rotation about the X axis by `theta3`, followed by a
rotation about the Y axis by `theta2`, and finally a rotation about the X axis
by `theta1`.
"""
immutable RotXYX{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotXYX}){X,Y,Z}(x::X, y::Y, z::Z) = RotXYX{promote_type(promote_type(X, Y), Z)}(x, y, z)


# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotXYX}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[2, 1], (-R[3, 1] + eps(t[1])) - eps(t[1]))  # TODO: handle denormal numbers better, as atan2(0,0) != atan2(0,-0)
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2((R[1, 2] * R[1, 2] + R[1, 3] * R[1, 3])^(1/2), R[1, 1]),
        atan2(- R[2, 3]*ct1 - R[3, 3]*st1, R[2, 2]*ct1 + R[3, 2]*st1))
end
@inline function Base.getindex{T}(r::RotXYX{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotXYX{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₂,        sinθ₁*sinθ₂,                        cosθ₁*-sinθ₂,
     sinθ₂*sinθ₃,  cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃,   sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃,
     sinθ₂*cosθ₃,  cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃,  sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotXYX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        cosθ₂*v[1] + sinθ₂*sinθ₃*v[2] + sinθ₂*cosθ₃*v[3],
        sinθ₁*sinθ₂*v[1] + (cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃)*v[2] + (cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃)*v[3],
        cosθ₁*-sinθ₂*v[1] + (sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃)*v[2] + (sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃)*v[3])
end

@inline Base.:*(r1::RotX, r2::RotYX) = RotXYX(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotXY, r2::RotX) = RotXYX(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotXYX, r2::RotX) = RotXYX(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotX, r2::RotXYX) = RotXYX(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotXYX) = RotXYX(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotXYX}) = RotXYX(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotXYX{T}}) = RotXYX{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotXYX)
    "3×3 RotXYX{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotXZX{T} <: Rotation{T}
    RotXZX(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" XZX Euler angle convention,
consisting of first a rotation about the X axis by `theta3`, followed by a
rotation about the Z axis by `theta2`, and finally a rotation about the X axis
by `theta1`.
"""
immutable RotXZX{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotXZX}){X,Y,Z}(x::X, y::Y, z::Z) = RotXZX{promote_type(promote_type(X, Y), Z)}(x, y, z)


# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotXZX}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[3, 1], R[2, 1])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2((R[1, 2] * R[1, 2] + R[1, 3] * R[1, 3])^(1/2), R[1, 1]),
        atan2(R[3, 2]*ct1 - R[2, 2]*st1, R[3, 3]*ct1 - R[2, 3]*st1))
end
@inline function Base.getindex{T}(r::RotXZX{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotXZX{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₂,         cosθ₁*sinθ₂,                        sinθ₁*sinθ₂,
     -sinθ₂*cosθ₃,  cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃,   sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃,
     sinθ₂*sinθ₃,   cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃,  sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotXZX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        cosθ₂*v[1] + -sinθ₂*cosθ₃*v[2] + sinθ₂*sinθ₃*v[3],
        cosθ₁*sinθ₂*v[1] + (cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃)*v[2] + (cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃)*v[3],
        sinθ₁*sinθ₂*v[1] + (sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[2] + (sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[3])
end

@inline Base.:*(r1::RotX, r2::RotZX) = RotXZX(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotXZ, r2::RotX) = RotXZX(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotXZX, r2::RotX) = RotXZX(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotX, r2::RotXZX) = RotXZX(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotXZX) = RotXZX(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotXZX}) = RotXZX(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotXZX{T}}) = RotXZX{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotXZX)
    "3×3 RotXZX{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotYXY{T} <: Rotation{T}
    RotYXY(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" YXY Euler angle convention,
consisting of first a rotation about the Y axis by `theta3`, followed by a
rotation about the X axis by `theta2`, and finally a rotation about the Y axis
by `theta1`.
"""
immutable RotYXY{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotYXY}){X,Y,Z}(x::X, y::Y, z::Z) = RotYXY{promote_type(promote_type(X, Y), Z)}(x, y, z)


# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotYXY}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[1, 2], R[3, 2])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2((R[2, 1] * R[2, 1] + R[2, 3] * R[2, 3])^(1/2), R[2, 2]),
        atan2(R[1, 3]*ct1 - R[3, 3]*st1, R[1, 1]*ct1 - R[3, 1]*st1))
end
@inline function Base.getindex{T}(r::RotYXY{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotYXY{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₃ + sinθ₁*cosθ₂*-sinθ₃,  sinθ₂*sinθ₃,   -sinθ₁*cosθ₃ + cosθ₁*cosθ₂*-sinθ₃,
     sinθ₁*sinθ₂,                       cosθ₂,         cosθ₁*sinθ₂,
     cosθ₁*sinθ₃ + sinθ₁*cosθ₂*cosθ₃,   -sinθ₂*cosθ₃,  -sinθ₁*sinθ₃ + cosθ₁*cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotYXY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        (cosθ₁*cosθ₃ + sinθ₁*cosθ₂*-sinθ₃)*v[1] + sinθ₁*sinθ₂*v[2] + (cosθ₁*sinθ₃ + sinθ₁*cosθ₂*cosθ₃)*v[3],
        sinθ₂*sinθ₃*v[1] + cosθ₂*v[2] + -sinθ₂*cosθ₃*v[3],
        (-sinθ₁*cosθ₃ + cosθ₁*cosθ₂*-sinθ₃)*v[1] + cosθ₁*sinθ₂*v[2] + (-sinθ₁*sinθ₃ + cosθ₁*cosθ₂*cosθ₃)*v[3])
end

@inline Base.:*(r1::RotY, r2::RotXY) = RotYXY(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotYX, r2::RotY) = RotYXY(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotYXY, r2::RotY) = RotYXY(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotY, r2::RotYXY) = RotYXY(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotYXY) = RotYXY(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotYXY}) = RotYXY(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotYXY{T}}) = RotYXY{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotYXY)
    "3×3 RotYXY{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotYZY{T} <: Rotation{T}
    RotYZY(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" YXY Euler angle convention,
consisting of first a rotation about the Y axis by `theta3`, followed by a
rotation about the Z axis by `theta2`, and finally a rotation about the Y axis
by `theta1`.
"""
immutable RotYZY{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotYZY}){X,Y,Z}(x::X, y::Y, z::Z) = RotYZY{promote_type(promote_type(X, Y), Z)}(x, y, z)


# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotYZY}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[3, 2], -R[1, 2])  # TODO: handle denormal numbers better, as atan2(0,0) != atan2(0,-0)
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2((R[2, 1] * R[2, 1] + R[2, 3] * R[2, 3])^(1/2), R[2, 2]),
        atan2(- R[3, 1]*ct1 - R[1, 1]*st1, R[3, 3]*ct1 + R[1, 3]*st1))
end
@inline function Base.getindex{T}(r::RotYZY{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotYZY{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₂*cosθ₃ + sinθ₁*-sinθ₃,  sinθ₂*cosθ₃,  -sinθ₁*cosθ₂*cosθ₃ + cosθ₁*-sinθ₃,
     cosθ₁*-sinθ₂,                      cosθ₂,        sinθ₁*sinθ₂,
     cosθ₁*cosθ₂*sinθ₃ + sinθ₁*cosθ₃,   sinθ₂*sinθ₃,  -sinθ₁*cosθ₂*sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotYZY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        (cosθ₁*cosθ₂*cosθ₃ + sinθ₁*-sinθ₃)*v[1] + cosθ₁*-sinθ₂*v[2] + (cosθ₁*cosθ₂*sinθ₃ + sinθ₁*cosθ₃)*v[3],
        sinθ₂*cosθ₃*v[1] + cosθ₂*v[2] + sinθ₂*sinθ₃*v[3],
        (-sinθ₁*cosθ₂*cosθ₃ + cosθ₁*-sinθ₃)*v[1] + sinθ₁*sinθ₂*v[2] + (-sinθ₁*cosθ₂*sinθ₃ + cosθ₁*cosθ₃)*v[3])
end

@inline Base.:*(r1::RotY, r2::RotZY) = RotYZY(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotYZ, r2::RotY) = RotYZY(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotYZY, r2::RotY) = RotYZY(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotY, r2::RotYZY) = RotYZY(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotYZY) = RotYZY(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotYZY}) = RotYZY(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotYZY{T}}) = RotYZY{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotYZY)
    "3×3 RotYZY{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotZXZ{T} <: Rotation{T}
    RotZXZ(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" ZXZ Euler angle convention,
consisting of first a rotation about the Z axis by `theta3`, followed by a
rotation about the X axis by `theta2`, and finally a rotation about the Z axis
by `theta1`.
"""
immutable RotZXZ{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotZXZ}){X,Y,Z}(x::X, y::Y, z::Z) = RotZXZ{promote_type(promote_type(X, Y), Z)}(x, y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotZXZ}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[1, 3], (-R[2, 3] + eps()) - eps())  # TODO: handle denormal numbers better, as atan2(0,0) != atan2(0,-0)
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2((R[3, 1] * R[3, 1] + R[3, 2] * R[3, 2])^(1/2), R[3, 3]),
        atan2(- R[1, 2]*ct1 - R[2, 2]*st1, R[1, 1]*ct1 + R[2, 1]*st1))
end

@inline function Base.getindex{T}(r::RotZXZ{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotZXZ{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃,   sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃,   sinθ₂*sinθ₃,
     cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃,  sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃,  sinθ₂*cosθ₃,
     sinθ₁*sinθ₂,                        cosθ₁*-sinθ₂,                      cosθ₂)
end

@inline function Base.:*(r::RotZXZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
         (cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃)*v[1] + (cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃)*v[2] + sinθ₁*sinθ₂*v[3],
         (sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃)*v[1] + (sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃)*v[2] + cosθ₁*-sinθ₂*v[3],
         sinθ₂*sinθ₃*v[1] + sinθ₂*cosθ₃*v[2] + cosθ₂*v[3])
end

@inline Base.:*(r1::RotZ, r2::RotXZ) = RotZXZ(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotZX, r2::RotZ) = RotZXZ(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotZXZ, r2::RotZ) = RotZXZ(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotZ, r2::RotZXZ) = RotZXZ(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotZXZ) = RotZXZ(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotZXZ}) = RotZXZ(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotZXZ{T}}) = RotZXZ{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotZXZ)
    "3×3 RotZXZ{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotZYZ{T} <: Rotation{T}
    RotZYZ(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" ZXZ Euler angle convention,
consisting of first a rotation about the Z axis by `theta3`, followed by a
rotation about the Y axis by `theta2`, and finally a rotation about the Z axis
by `theta1`.
"""
immutable RotZYZ{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotZYZ}){X,Y,Z}(x::X, y::Y, z::Z) = RotZYZ{promote_type(promote_type(X, Y), Z)}(x, y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotZYZ}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[2, 3], R[1, 3])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2((R[3, 1] * R[3, 1] + R[3, 2] * R[3, 2])^(1/2), R[3, 3]),
        atan2(R[2, 1]*ct1 - R[1, 1]*st1, R[2, 2]*ct1 - R[1, 2]*st1))
end

@inline function Base.getindex{T}(r::RotZYZ{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotZYZ{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃,   sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃,   -sinθ₂*cosθ₃,
     cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃,  sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃,  sinθ₂*sinθ₃,
     cosθ₁*sinθ₂,                        sinθ₁*sinθ₂,                       cosθ₂)
end

@inline function Base.:*(r::RotZYZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        (cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃)*v[1] + (cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃)*v[2] + cosθ₁*sinθ₂*v[3],
        (sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[1] + (sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[2] + sinθ₁*sinθ₂*v[3],
        -sinθ₂*cosθ₃*v[1] + sinθ₂*sinθ₃*v[2] + cosθ₂*v[3])
end

@inline Base.:*(r1::RotZ, r2::RotYZ) = RotZYZ(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotZY, r2::RotZ) = RotZYZ(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotZYZ, r2::RotZ) = RotZYZ(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotZ, r2::RotZYZ) = RotZYZ(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotZYZ) = RotZYZ(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotZYZ}) = RotZYZ(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotZYZ{T}}) = RotZYZ{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotZYZ)
    "3×3 RotZYZ{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end



###############################
# Tait-Bryant Euler Rotations #
###############################

"""
    immutable RotXYZ{T} <: Rotation{T}
    RotXYZ(theta1, theta2, theta3)
    RotXYZ(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" XYZ Euler angle
convention, consisting of first a rotation about the Z axis by `theta3`,
followed by a rotation about the Y axis by `theta2`, and finally a rotation
about the X axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in XYZ order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
immutable RotXYZ{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotXYZ}){X,Y,Z}(x::X, y::Y, z::Z) = RotXYZ{promote_type(promote_type(X, Y), Z)}(x, y, z)
@inline (::Type{Rot}){Rot<:RotXYZ}(; roll=0, pitch=0, yaw=0) = Rot(roll, pitch, yaw)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotXYZ}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(-R[2, 3], R[3, 3])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2(R[1, 3], (R[1, 1] * R[1, 1] + R[1, 2] * R[1, 2])^(1/2)),
        atan2(R[2, 1]*ct1 + R[3, 1]*st1, R[2, 2]*ct1 + R[3, 2]*st1))
end

@inline function Base.getindex{T}(r::RotXYZ{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotXYZ{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₂*cosθ₃,   sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃,   cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃,
     cosθ₂*-sinθ₃,  sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃,  cosθ₁*-sinθ₂*-sinθ₃ + sinθ₁*cosθ₃,
     sinθ₂,         -sinθ₁*cosθ₂,                      cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotXYZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        cosθ₂*cosθ₃*v[1] + cosθ₂*-sinθ₃*v[2] + sinθ₂*v[3],
        (sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[1] + (sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[2] + -sinθ₁*cosθ₂*v[3],
        (cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃)*v[1] + (cosθ₁*sinθ₂*sinθ₃ + sinθ₁*cosθ₃)*v[2] + cosθ₁*cosθ₂*v[3])
end

@inline Base.:*(r1::RotX, r2::RotYZ) = RotXYZ(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotXY, r2::RotZ) = RotXYZ(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotXYZ, r2::RotZ) = RotXYZ(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotX, r2::RotXYZ) = RotXYZ(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotXYZ) = RotZYX(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotXYZ}) = RotXYZ(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotXYZ{T}}) = RotXYZ{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotXYZ)
    "3×3 RotXYZ{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotZYX{T} <: Rotation{T}
    RotZYX(theta1, theta2, theta3)
    RotZYX(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" ZYX Euler angle
convention, consisting of first a rotation about the X axis by `theta3`,
followed by a rotation about the Y axis by `theta2`, and finally a rotation
about the Z axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in ZYX order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
immutable RotZYX{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotZYX}){X,Y,Z}(x::X, y::Y, z::Z) = RotZYX{promote_type(promote_type(X, Y), Z)}(x, y, z)
@inline (::Type{Rot}){Rot<:RotZYX}(; roll=0, pitch=0, yaw=0) = Rot(yaw, pitch, roll)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotZYX}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[2, 1], R[1, 1])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2(-R[3, 1], (R[3, 2] * R[3, 2] + R[3, 3] * R[3, 3])^(1/2)),
        atan2(R[1, 3]*st1 - R[2, 3]*ct1, R[2, 2]*ct1 - R[1, 2]*st1))
end

@inline function Base.getindex{T}(r::RotZYX{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotZYX{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    ( cosθ₁*cosθ₂,                       sinθ₁*cosθ₂,                      -sinθ₂,
     -sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃,   cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃,   cosθ₂*sinθ₃,
      sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃,   cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃,  cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotZYX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        cosθ₁*cosθ₂*v[1] + (-sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃)*v[2] + (sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃)*v[3],
        sinθ₁*cosθ₂*v[1] + (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃)*v[2] + (cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃)*v[3],
        -sinθ₂*v[1] + cosθ₂*sinθ₃*v[2] + cosθ₂*cosθ₃*v[3])
end

@inline Base.:*(r1::RotZ, r2::RotYX) = RotZYX(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotZY, r2::RotX) = RotZYX(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotZYX, r2::RotX) = RotZYX(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotZ, r2::RotZYX) = RotZYX(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotZYX) = RotXYZ(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotZYX}) = RotZYX(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotZYX{T}}) = RotZYX{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotZYX)
    "3×3 RotZYX{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotXZY{T} <: Rotation{T}
    RotXZY(theta1, theta2, theta3)
    RotXZY(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" XZY Euler angle
convention, consisting of first a rotation about the Y axis by `theta3`,
followed by a rotation about the Z axis by `theta2`, and finally a rotation
about the X axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in XZY order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
immutable RotXZY{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotXZY}){X,Y,Z}(x::X, y::Y, z::Z) = RotXZY{promote_type(promote_type(X, Y), Z)}(x, y, z)
@inline (::Type{Rot}){Rot<:RotXZY}(; roll=0, pitch=0, yaw=0) = Rot(roll, yaw, pitch)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotXZY}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[3, 2], R[2, 2])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2(-R[1, 2], (R[1, 1] * R[1, 1] + R[1, 3] * R[1, 3])^(1/2)),
        atan2(R[2, 1]*st1 - R[3, 1]*ct1, R[3, 3]*ct1 - R[2, 3]*st1))
end

@inline function Base.getindex{T}(r::RotXZY{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotXZY{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    ( cosθ₂*cosθ₃,  cosθ₁*sinθ₂*cosθ₃ + sinθ₁*sinθ₃,   sinθ₁*sinθ₂*cosθ₃ + cosθ₁*-sinθ₃,
     -sinθ₂,        cosθ₁*cosθ₂,                       sinθ₁*cosθ₂,
      cosθ₂*sinθ₃,  cosθ₁*sinθ₂*sinθ₃ + -sinθ₁*cosθ₃,  sinθ₁*sinθ₂*sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotXZY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        cosθ₂*cosθ₃*v[1] + -sinθ₂*v[2] + cosθ₂*sinθ₃*v[3],
        (cosθ₁*sinθ₂*cosθ₃ + sinθ₁*sinθ₃)*v[1] + cosθ₁*cosθ₂*v[2] + (cosθ₁*sinθ₂*sinθ₃ + -sinθ₁*cosθ₃)*v[3],
        (sinθ₁*sinθ₂*cosθ₃ + cosθ₁*-sinθ₃)*v[1] + sinθ₁*cosθ₂*v[2] + (sinθ₁*sinθ₂*sinθ₃ + cosθ₁*cosθ₃)*v[3])
end

@inline Base.:*(r1::RotX, r2::RotZY) = RotXZY(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotXZ, r2::RotY) = RotXZY(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotXZY, r2::RotY) = RotXZY(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotX, r2::RotXZY) = RotXZY(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotXZY) = RotYZX(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotXZY}) = RotXZY(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotXZY{T}}) = RotXZY{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotXZY)
    "3×3 RotXZY{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotYZX{T} <: Rotation{T}
    RotYZX(theta1, theta2, theta3)
    RotYZX(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" YZX Euler angle
convention, consisting of first a rotation about the X axis by `theta3`,
followed by a rotation about the Z axis by `theta2`, and finally a rotation
about the Y axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in YZX order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
immutable RotYZX{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotYZX}){X,Y,Z}(x::X, y::Y, z::Z) = RotYZX{promote_type(promote_type(X, Y), Z)}(x, y, z)
@inline (::Type{Rot}){Rot<:RotYZX}(; roll=0, pitch=0, yaw=0) = Rot(pitch, yaw, roll)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotYZX}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(-R[3, 1], R[1, 1])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2(R[2, 1], (R[2, 2] * R[2, 2] + R[2, 3] * R[2, 3])^(1/2)),
        atan2(R[3, 2]*ct1 + R[1, 2]*st1, R[3, 3]*ct1 + R[1, 3]*st1))
end

@inline function Base.getindex{T}(r::RotYZX{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotYZX{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₂,                        sinθ₂,         -sinθ₁*cosθ₂,
     cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃,   cosθ₂*cosθ₃,   sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃,
     cosθ₁*sinθ₂*sinθ₃ + sinθ₁*cosθ₃,    cosθ₂*-sinθ₃,  sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotYZX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        cosθ₁*cosθ₂*v[1] + (cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃)*v[2] + (cosθ₁*-sinθ₂*-sinθ₃ + sinθ₁*cosθ₃)*v[3],
        sinθ₂*v[1] + cosθ₂*cosθ₃*v[2] + cosθ₂*-sinθ₃*v[3],
        -sinθ₁*cosθ₂*v[1] + (sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[2] + (sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[3])
end

@inline Base.:*(r1::RotY, r2::RotZX) = RotYZX(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotYZ, r2::RotX) = RotYZX(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotYZX, r2::RotX) = RotYZX(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotY, r2::RotYZX) = RotYZX(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotYZX) = RotXZY(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotYZX}) = RotYZX(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotYZX{T}}) = RotYZX{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotYZX)
    "3×3 RotYZX{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotYXZ{T} <: Rotation{T}
    RotYXZ(theta1, theta2, theta3)
    RotYXZ(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" YXZ Euler angle
convention, consisting of first a rotation about the Z axis by `theta3`,
followed by a rotation about the X axis by `theta2`, and finally a rotation
about the Y axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in YXZ order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
immutable RotYXZ{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotYXZ}){X,Y,Z}(x::X, y::Y, z::Z) = RotYXZ{promote_type(promote_type(X, Y), Z)}(x, y, z)
@inline (::Type{Rot}){Rot<:RotYXZ}(; roll=0, pitch=0, yaw=0) = Rot(pitch, roll, yaw)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotYXZ}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(R[1, 3], R[3, 3])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2(-R[2, 3], (R[2, 1] * R[2, 1] + R[2, 2] * R[2, 2])^(1/2)),
        atan2(R[3, 2]*st1 - R[1, 2]*ct1, R[1, 1]*ct1 - R[3, 1]*st1))
end

@inline function Base.getindex{T}(r::RotYXZ{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotYXZ{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃,   cosθ₂*sinθ₃,  -sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃,
     cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃,  cosθ₂*cosθ₃,  sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃,
     sinθ₁*cosθ₂,                       -sinθ₂,       cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotYXZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃)*v[1] + (cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃)*v[2] + sinθ₁*cosθ₂*v[3],
        cosθ₂*sinθ₃*v[1] + cosθ₂*cosθ₃*v[2] + -sinθ₂*v[3],
        (-sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃)*v[1] + (sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃)*v[2] + cosθ₁*cosθ₂*v[3])
end

@inline Base.:*(r1::RotY, r2::RotXZ) = RotYXZ(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotYX, r2::RotZ) = RotYXZ(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotYXZ, r2::RotZ) = RotYXZ(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotY, r2::RotYXZ) = RotYXZ(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotYXZ) = RotZXY(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotYXZ}) = RotYXZ(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotYXZ{T}}) = RotYXZ{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotYXZ)
    "3×3 RotYXZ{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end


"""
    immutable RotZXY{T} <: Rotation{T}
    RotZXY(theta1, theta2, theta3)
    RotZXY(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" ZXY Euler angle
convention, consisting of first a rotation about the Y axis by `theta3`,
followed by a rotation about the X axis by `theta2`, and finally a rotation
about the Z axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in ZXY order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
immutable RotZXY{T} <: Rotation{T}
    theta1::T
    theta2::T
    theta3::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 3 inputs (not 9).
@inline (::Type{RotZXY}){X,Y,Z}(x::X, y::Y, z::Z) = RotZXY{promote_type(promote_type(X, Y), Z)}(x, y, z)
@inline (::Type{Rot}){Rot<:RotZXY}(; roll=0, pitch=0, yaw=0) = Rot(yaw, roll, pitch)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline function (::Type{Rot}){Rot <: RotZXY}(t::NTuple{9})
    R = SMatrix{3,3}(t)

    t1 = atan2(-R[1, 2], R[2, 2])
    ct1, st1 = cos(t1), sin(t1)

    Rot(t1,
        atan2(R[3, 2], (R[3, 1] * R[3, 1] + R[3, 3] * R[3, 3])^(1/2)),
        atan2(R[1, 3]*ct1 + R[2, 3]*st1, R[1, 1]*ct1 + R[2, 1]*st1))
end

@inline function Base.getindex{T}(r::RotZXY{T}, i::Integer)
    Tuple(r)[i] # Slow...
end

@inline function Base.convert{T}(::Type{Tuple}, r::RotZXY{T})
    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    # transposed representation
    ( cosθ₁*cosθ₃ + sinθ₁*sinθ₂*-sinθ₃,  sinθ₁*cosθ₃ + cosθ₁*-sinθ₂*-sinθ₃,  cosθ₂*-sinθ₃,
     -sinθ₁*cosθ₂,                       cosθ₁*cosθ₂,                        sinθ₂,
      cosθ₁*sinθ₃ + sinθ₁*sinθ₂*cosθ₃,   sinθ₁*sinθ₃ + cosθ₁*-sinθ₂*cosθ₃,   cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotZXY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    T = promote_type(eltype(r), eltype(v))

    sinθ₁ = sin(r.theta1)
    cosθ₁ = cos(r.theta1)
    sinθ₂ = sin(r.theta2)
    cosθ₂ = cos(r.theta2)
    sinθ₃ = sin(r.theta3)
    cosθ₃ = cos(r.theta3)

    return similar_type(v,T)(
        (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*-sinθ₃)*v[1] + -sinθ₁*cosθ₂*v[2] + (cosθ₁*sinθ₃ + sinθ₁*sinθ₂*cosθ₃)*v[3],
        (sinθ₁*cosθ₃ + cosθ₁*-sinθ₂*-sinθ₃)*v[1] + cosθ₁*cosθ₂*v[2] + (sinθ₁*sinθ₃ + cosθ₁*-sinθ₂*cosθ₃)*v[3],
         cosθ₂*-sinθ₃*v[1] + sinθ₂*v[2] + cosθ₂*cosθ₃*v[3])
end

@inline Base.:*(r1::RotZ, r2::RotXY) = RotZXY(r1.theta, r2.theta1, r2.theta2)
@inline Base.:*(r1::RotZX, r2::RotY) = RotZXY(r1.theta1, r1.theta2, r2.theta)
@inline Base.:*(r1::RotZXY, r2::RotY) = RotZXY(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
@inline Base.:*(r1::RotZ, r2::RotZXY) = RotZXY(r1.theta + r2.theta1, r2.theta2, r2.theta3)

@inline inv(r::RotZXY) = RotYXZ(-r.theta3, -r.theta2, -r.theta1)

# define null rotations for convenience
@inline eye(::Type{RotZXY}) = RotZXY(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RotZXY{T}}) = RotZXY{T}(zero(T), zero(T), zero(T))

# Need a good output representation
function Base.summary(r::RotZXY)
    "3×3 RotZXY{$(eltype(r))}($(r.theta1), $(r.theta2), $(r.theta3))"
end

################################################################################
################################################################################



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
