
#
# definition for rotating points
# 

"""
Rotate the point X using the rotation parameters specified in R (i.e. R * X)

rotate_point(R::T, X)

"""
rotate_point{T <: RotationTypes, U <: AbstractFloat}(R::T, X::Vec{3,U})  =  rotate_point(RotMatrix(R), X)    # generic version, go through a rotation matrix

# in case of mutable input (not the recomended form)
rotate_point{T <: RotationTypes, U <: AbstractFloat}(R::T, X::Vector{U}) = rotate_point(R, Vec{3,U}(X[1], X[2], X[3]))   # generic version, go through a rotation matrix

# cases where we start with a rotation matrix with matched element types
rotate_point{T <: AbstractFloat}(R::RotMatrix{T}, X::Vec{3,T}) =   R * X

# cases with unmatched types
rotate_point{T <: AbstractFloat, U <: AbstractFloat}(R::RotMatrix{T}, X::Vec{3,U}) =  R * X  # should this have return type U?

#=

################################################
# Using Quaternions to do the rotation
# (benchmark says not worth)
################################################

rotate_point{T <: AbstractFloat}(q::Quaternion{T}, X::Vector{T}) = imag(q*Quaternion(X)*conj(q))
function rotate_point{T <: AbstractFloat}(q::Quaternion{T}, X::Vec{3,T}) 
    qo = (q*Quaternion(X)*conj(q))
    return Vec{3,T}(qo.v1, qo.v2, qo.v3)
end



using Quaternions
function benchmark_rotation(n::Int=1_000_000)

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
    @time rot_test!(Xo, Qv, X)

    println("Rotate via Quaternions")
    @time quat_test!(Xo, Qv, X)

end
=#
