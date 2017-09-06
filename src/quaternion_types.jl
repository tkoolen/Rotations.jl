"""
    struct Quat{T} <: Rotation{3,T}
    Quat(w, x, y, z)

The `Quat` type is a 3×3 matrix representation of a normalized quaternion.
They allow you to transparently use (fast) quaternion algebra to store, compose
and invert 3D rotations, while at the same time letting you apply rotations
through matrix-vector multiplication.

Note: by default, the constructor will renormalize the input so that the quaternion
has length 1 (w² + x² + y² + z² = 1), and the rotation matrix is orthogonal.

Renormalization can be skipped by passing `false` as an additional constructor
argument, in which case the user provides the guarantee that the input arguments
represent a unit quaternion. Operations on an unnormalized `Quat`, created by
skipping renormalization in this fashion, are not guaranteed to do anything sensible.
"""
struct Quat{T} <: Rotation{3,T}
    w::T
    x::T
    y::T
    z::T

    @inline function Quat{T}(w, x, y, z, normalize::Bool = true) where {T}
        if normalize
            norm = sqrt(w*w + x*x + y*y + z*z)
            new(w/norm, x/norm, y/norm, z/norm)
        else
            new(w, x, y, z)
        end
    end

    Quat{T}(q::Quat) where {T} = new{T}(q.w, q.x, q.y, q.z)
end

@inline function Quat(w::W, x::X, y::Y, z::Z, normalize::Bool = true) where {W, X, Y, Z}
    Quat{promote_type(promote_type(promote_type(W, X), Y), Z)}(w, x, y, z, normalize)
end
@inline Quat(q::Quat{T}) where {T} = Quat{T}(q)

@inline convert(::Type{Q}, q::Quat) where {Q<:Quat} = Q(q)
@inline convert(::Type{Q}, q::Q) where {Q<:Quat} = q

# These 3 functions are enough to satisfy the entire StaticArrays interface:
function (::Type{Q})(t::NTuple{9}) where Q<:Quat
    q = Q(sqrt(abs(1  + t[1] + t[5] + t[9])) / 2,
               copysign(sqrt(abs(1 + t[1] - t[5] - t[9]))/2, t[6] - t[8]),
               copysign(sqrt(abs(1 - t[1] + t[5] - t[9]))/2, t[7] - t[3]),
               copysign(sqrt(abs(1 - t[1] - t[5] + t[9]))/2, t[2] - t[4]))
end

function Base.getindex(q::Quat, i::Int)
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

function Tuple(q::Quat)
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
@inline Base.eye(::Type{Quat{T}}) where {T} = Quat{T}(one(T), zero(T), zero(T), zero(T))

"""
    rotation_between(from, to)

Compute the quaternion that rotates vector `from` so that it aligns with vector
`to`, along the geodesic (shortest path).
"""
rotation_between(from::AbstractVector, to::AbstractVector) = rotation_between(SVector{3}(from), SVector{3}(to))
function rotation_between(from::SVector{3}, to::SVector{3})
    # Robustified version of implementation from https://www.gamedev.net/topic/429507-finding-the-quaternion-betwee-two-vectors/#entry3856228
    normprod = sqrt(dot(from, from) * dot(to, to))
    T = typeof(normprod)
    normprod < eps(T) && throw(ArgumentError("Input vectors must be nonzero."))
    w = normprod + dot(from, to)
    v = abs(w) < 100 * eps(T) ? perpendicular_vector(from) : cross(from, to)
    @inbounds return Quat(w, v[1], v[2], v[3]) # relies on normalization in constructor
end

################################################################################
################################################################################
"""
    struct SPQuat{T} <: Rotation{3,T}
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
struct SPQuat{T} <: Rotation{3,T}
    x::T
    y::T
    z::T

    # TODO should we enforce norm <= 1?
    SPQuat{T}(x, y, z) where {T} = new{T}(x, y, z)
    SPQuat{T}(spq::SPQuat) where {T} = new{T}(spq.x, spq.y, spq.z)
end

@inline SPQuat(x::X, y::Y, z::Z) where {X,Y,Z} = SPQuat{promote_type(promote_type(X, Y), Z)}(x, y, z)
@inline SPQuat(spq::SPQuat{T}) where {T} = SPQuat{T}(spq)

@inline convert(::Type{SPQ}, spq::SPQuat) where {SPQ<:SPQuat} = SPQ(spq)
@inline convert(::Type{SPQ}, spq::SPQ) where {SPQ<:SPQuat} = spq

# These functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{SPQ})(t::NTuple{9}) where {SPQ <: SPQuat} = SPQ(Quat(t))
@inline Base.getindex(spq::SPQuat, i::Int) = Quat(spq)[i]
@inline Tuple(spq::SPQuat) = Tuple(Quat(spq))

@inline function Base.convert(::Type{Q}, spq::SPQuat) where Q <: Quat
    # Both the sign and norm of the Quat is automatically dealt with in its inner constructor
    return Q(1 - (spq.x*spq.x + spq.y*spq.y + spq.z*spq.z), 2*spq.x, 2*spq.y, 2*spq.z)
end

@inline function Base.convert(::Type{SPQ}, q::Quat) where SPQ <: SPQuat
    alpha2 = (1 - q.w) / (1 + q.w) # <= 1 since q.w >= 0
    spq = SPQ(q.x * (alpha2 + 1)*0.5,  q.y * (alpha2 + 1)*0.5, q.z * (alpha2 + 1)*0.5)
end

@inline Base.:*(spq::SPQuat, x::StaticVector) = Quat(spq) * x

@inline Base.:*(spq::SPQuat, r::Rotation) = Quat(spq) * r
@inline Base.:*(spq::SPQuat, r::RotMatrix) = Quat(spq) * r
@inline Base.:*(r::Rotation, spq::SPQuat) = r * Quat(spq)
@inline Base.:*(r::RotMatrix, spq::SPQuat) = r * Quat(spq)
@inline Base.:*(spq1::SPQuat, spq2::SPQuat) = Quat(spq1) * Quat(spq2)

@inline Base.inv(spq::SPQuat) = SPQuat(-spq.x, -spq.y, -spq.z)

@inline Base.eye(::Type{SPQuat}) = SPQuat(0.0, 0.0, 0.0)
@inline Base.eye(::Type{SPQuat{T}}) where {T} = SPQuat{T}(zero(T), zero(T), zero(T))

# rotation properties
@inline rotation_angle(spq::SPQuat) = rotation_angle(Quat(spq))
@inline rotation_axis(spq::SPQuat) = rotation_axis(Quat(spq))
