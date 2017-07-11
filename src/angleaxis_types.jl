################################################################################
################################################################################
"""
    struct AngleAxis{T} <: Rotation{3,T}
    AngleAxis(Θ, x, y, z)

A 3×3 rotation matrix parameterized by a 3D rotation by angle θ about an
arbitrary axis `[x, y, z]`.

Note that the axis is not unique for θ = 0, and that this parameterization does
not continuously map the neighbourhood of the null rotation (and therefore
might not be suitable for autodifferentation and optimization purposes).

Note: by default, the constructor will renormalize the input so that the axis
has length 1 (x² + y² + z² = 1).

Renormalization can be skipped by passing `Val{false}()` as an additional constructor
argument, in which case the user provides the guarantee that the input arguments
represent a normalized rotation axis. Operations on an `AngleAxis` with a rotation
axis that does not have unit norm, created by skipping renormalization in this fashion,
are not guaranteed to do anything sensible.
"""
struct AngleAxis{T} <: Rotation{3,T}
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T

    @inline function AngleAxis{T}(θ, x, y, z, ::Val{Normalize} = Val{true}()) where {T, Normalize}
        if Normalize
            # Not sure what to do with theta?? Should it become theta * norm ?
            norm = sqrt(x*x + y*y + z*z)
            new(θ, x/norm, y/norm, z/norm)
        else
            new(θ, x, y, z)
        end
    end
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline function AngleAxis(θ::Θ, x::X, y::Y, z::Z, normalize::Val = Val{true}()) where {Θ,X,Y,Z}
    AngleAxis{promote_type(promote_type(promote_type(Θ, X), Y), Z)}(θ, x, y, z, normalize)
end

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{AA}){AA <: AngleAxis}(t::NTuple{9}) = AA(Quat(t))
@inline Base.getindex(aa::AngleAxis, i::Int) = Quat(aa)[i]

@inline function Base.convert{R <: RotMatrix}(::Type{R}, aa::AngleAxis)
    # Rodrigues' rotation formula.
    T = eltype(aa)

    s = sin(aa.theta)
    c = cos(aa.theta)
    c1 = one(T) - c

    c1x2 = c1 * aa.axis_x^2
    c1y2 = c1 * aa.axis_y^2
    c1z2 = c1 * aa.axis_z^2

    c1xy = c1 * aa.axis_x * aa.axis_y
    c1xz = c1 * aa.axis_x * aa.axis_z
    c1yz = c1 * aa.axis_y * aa.axis_z

    sx = s * aa.axis_x
    sy = s * aa.axis_y
    sz = s * aa.axis_z

    # Note that the RotMatrix constructor argument order makes this look transposed:
    R(one(T) - c1y2 - c1z2, c1xy + sz, c1xz - sy,
      c1xy - sz, one(T) - c1x2 - c1z2, c1yz + sx,
      c1xz + sy, c1yz - sx, one(T) - c1x2 - c1y2)
end

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
@inline Base.:*(aa::AngleAxis, r::RotMatrix) = Quat(aa) * r
@inline Base.:*(aa::AngleAxis, r::SPQuat) = Quat(aa) * r
@inline Base.:*(r::Rotation, aa::AngleAxis) = r * Quat(aa)
@inline Base.:*(r::RotMatrix, aa::AngleAxis) = r * Quat(aa)
@inline Base.:*(r::SPQuat, aa::AngleAxis) = r * Quat(aa)
@inline Base.:*(aa1::AngleAxis, aa2::AngleAxis) = Quat(aa1) * Quat(aa2)

@inline inv(aa::AngleAxis) = AngleAxis(-aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)
@inline Base.:^(aa::AngleAxis, t::Real) = AngleAxis(aa.theta*t, aa.axis_x, aa.axis_y, aa.axis_z)
@inline Base.:^(aa::AngleAxis, t::Integer) = AngleAxis(aa.theta*t, aa.axis_x, aa.axis_y, aa.axis_z) # to avoid ambiguity


# define null rotations for convenience
@inline eye(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0.0, 0.0)
@inline eye{T}(::Type{AngleAxis{T}}) = AngleAxis{T}(zero(T), one(T), zero(T), zero(T))

# accessors
@inline rotation_angle(aa::AngleAxis) = aa.theta #  - floor((aa.theta+pi) / (2*pi)) * 2*pi
@inline rotation_axis(aa::AngleAxis) = SVector(aa.axis_x, aa.axis_y, aa.axis_z)


################################################################################
################################################################################
"""
    struct RodriguesVec{T} <: Rotation{3,T}
    RodriguesVec(sx, sy, sz)

Rodrigues vector parameterization of a 3×3 rotation matrix. The direction of the
vector [sx, sy, sz] defines the axis of rotation, and the rotation angle is
given by its norm.
"""
struct RodriguesVec{T} <: Rotation{3,T}
    sx::T
    sy::T
    sz::T
end

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline (::Type{RodriguesVec}){X,Y,Z}(x::X, y::Y, z::Z) = RodriguesVec{promote_type(promote_type(X, Y), Z)}(x, y, z)

# These 2 functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RV}){RV <: RodriguesVec}(t::NTuple{9}) = RV(Quat(t))
@inline Base.getindex(aa::RodriguesVec, i::Int) = Quat(aa)[i]

# define its interaction with other angle representations
@inline Base.convert{R <: RotMatrix}(::Type{R}, rv::RodriguesVec) = convert(R, AngleAxis(rv))

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
@inline Base.:*(rv::RodriguesVec, r::RotMatrix) = Quat(rv) * r
@inline Base.:*(rv::RodriguesVec, r::SPQuat) = Quat(rv) * r
@inline Base.:*(rv::RodriguesVec, r::AngleAxis) = Quat(rv) * r
@inline Base.:*(r::Rotation, rv::RodriguesVec) = r * Quat(rv)
@inline Base.:*(r::RotMatrix, rv::RodriguesVec) = r * Quat(rv)
@inline Base.:*(r::SPQuat, rv::RodriguesVec) = r * Quat(rv)
@inline Base.:*(r::AngleAxis, rv::RodriguesVec) = r * Quat(rv)
@inline Base.:*(rv1::RodriguesVec, rv2::RodriguesVec) = Quat(rv1) * Quat(rv2)

@inline inv(rv::RodriguesVec) = RodriguesVec(-rv.sx, -rv.sy, -rv.sz)
@inline Base.:^(rv::RodriguesVec, t::Real) = RodriguesVec(rv.sx*t, rv.sy*t, rv.sz*t)
@inline Base.:^(rv::RodriguesVec, t::Integer) = RodriguesVec(rv.sx*t, rv.sy*t, rv.sz*t) # to avoid ambiguity




# rotation properties
@inline rotation_angle(rv::RodriguesVec) = sqrt(rv.sx * rv.sx + rv.sy * rv.sy + rv.sz * rv.sz)
function rotation_axis(rv::RodriguesVec)     # what should this return for theta = 0?
    theta = rotation_angle(rv)
    return (theta > 0 ? SVector(rv.sx / theta, rv.sy / theta, rv.sz / theta) : SVector(one(theta), zero(theta), zero(theta)))
end

# define null rotations for convenience
@inline eye(::Type{RodriguesVec}) = RodriguesVec(0.0, 0.0, 0.0)
@inline eye{T}(::Type{RodriguesVec{T}}) = RodriguesVec{T}(zero(T), zero(T), zero(T))
