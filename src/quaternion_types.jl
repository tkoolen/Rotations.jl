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
@inline (::Type{SPQ}){SPQ <: SPQuat}(t::NTuple{9}) = SPQ(Quat(t))
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

@inline Base.:*(spq::SPQuat, x::StaticVector) = Quat(spq) * x

@inline Base.:*(spq::SPQuat, r::Rotation) = Quat(spq) * r
@inline Base.:*(spq::SPQuat, r::RotMatrix) = Quat(spq) * r
@inline Base.:*(r::Rotation, spq::SPQuat) = r * Quat(spq)
@inline Base.:*(r::RotMatrix, spq::SPQuat) = r * Quat(spq)
@inline Base.:*(spq1::SPQuat, spq2::SPQuat) = Quat(spq1) * Quat(spq2)

@inline Base.inv(spq::SPQuat) = SPQuat(-spq.x, -spq.y, -spq.z)

@inline Base.eye(::Type{SPQuat}) = SPQuat(0.0, 0.0, 0.0)
@inline Base.eye{T}(::Type{SPQuat{T}}) = SPQuat{T}(zero(T), zero(T), zero(T))

# rotation properties
@inline rotation_angle(spq::SPQuat) = rotation_angle(Quat(spq))
@inline rotation_axis(spq::SPQuat) = rotation_axis(Quat(spq))
