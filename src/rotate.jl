
#
# definition for rotating points
#

"""
Rotate the point X using the rotation parameters specified in R (i.e. R * X)

rotate(R, X)

"""
@inline rotate{T <: RotationTypes, U <: Real}(R::T, X::FixedVector{3,U})  =  rotate(RotMatrix(R), X)    # generic version, go through a rotation matrix

# in case of mutable input (not the recomended form)
@inline rotate{T <: RotationTypes, U <: Real}(R::T, X::Vector{U}) = Array(rotate(R, Vec{3,U}(X[1], X[2], X[3])))   # generic version, go through a rotation matrix

# cases where we start with a rotation matrix with matched element types
@inline rotate{T <: Real}(R::RotMatrix{T}, X::FixedVector{3,T}) =   R * X

# this will actually throw
# @inline rotate_point{T <: Real, U <: Real}(R::RotMatrix{T}, X::Vec{3,U}) =   R * X  # Revert to this once changes to FixedSizeArrays are in
@inline rotate{T, U <: Real}(R::RotMatrix{T}, X::FixedVector{3,U}) = convert(RotMatrix{promote_type(T,U)}, R) * convert(Vec{3, promote_type(T,U)}, X)



################################################
# Using Quaternions to do the rotation
# (benchmark says not worth)
################################################

#=
rotate{T <: AbstractFloat}(q::Quaternion{T}, X::Vector{T}) = imag(q*Quaternion(X)*conj(q))
function rotate{T <: AbstractFloat}(q::Quaternion{T}, X::Vec{3,T})
    qo = (q*Quaternion(X)*conj(q))
    return Vec{3,T}(qo.v1, qo.v2, qo.v3)
end

function benchmark_mutable(n::Int=1_000_000)
    function rotate_mutable(R, X, n)
        Xb, Xo = zeros(3), zeros(3)
        for i = 1:n
            A_mul_B!(Xb, R, X)
            # @inbounds Xo[1] += Xb[1]; @inbounds Xo[2] += Xb[2]; @inbounds Xo[3] += Xb[3];
        end
        return Xo
    end

    function rotate_immutable(R, X, n)
        Xo = Vec(0.0,0,0)
        for i = 1:n
            Xb = R * X
            # Xo += Xb
        end
        return Xo
    end

    # Initialise
    R_mute, R_immute = eye(3), RotMatrix(eye(3))
    X_mute, X_immute = zeros(3), Vec(0.0,0,0)

    # and test
    rotate_mutable(R_mute, X_mute, 1)
    println("Rotating using mutables")
    @time Xo = rotate_mutable(R_mute, X_mute, n)

    rotate_immutable(R_immute, X_immute, 1)
    println("Rotating using immutables")
    @time Xo = rotate_immutable(R_immute, X_immute, n)

end


using Quaternions
function benchmark_quaternion_rotation(n::Int=1_000_000)

    Qv = [nquatrand()::Quaternion{Float64} for i in 1:n]

    X = Vec{3,Float64}(randn(3))
    Xo = Vector{Vec{3,Float64}}(n)

    function quat_test!(Xo, Qv, X)
        for (i, q) in enumerate(Qv)
            Xo[i] = rotate_point(q, X)
        end
    end

    function rot_test!(Xo, Qv, X)
        for (i, q) in enumerate(Qv)
            Xo[i] = rotate_point(RotMatrix(q), X)
        end
    end

    println("Rotate via RotMatrix")
    rot_test!(Xo, Qv[1:1], X)
    @time rot_test!(Xo, Qv, X)

    println("Rotate via Quaternions")
    quat_test!(Xo, Qv[1:1], X)
    @time quat_test!(Xo, Qv, X)

end
=#

