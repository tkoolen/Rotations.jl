
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
rotate_point{T <: AbstractFloat, U <: AbstractFloat}(R::RotMatrix{T}, X::Vec{3,U}) =  convert(RotMatrix{U}, R) * X

# ending up with X as Ints seems inevitable
rotate_point{T <: AbstractFloat, U <: Integer}(R::RotMatrix{T}, X::Vec{3,U}) = R * convert(Vec{3, T}, X)





# TODO: add methods for quaternions



