
"""
Function to convert bewteen Rotation parameterizations

Args:

    T         -  the rotation type to convert to

    X         -  the input rotation parameterisation

"""
convert_rotation{T}(::Type{T}, X::T) = X   # do nothing yay!


# if we're in this version we probably need to match up parameters of the two types
function convert_rotation{T,U}(::Type{T}, X::U) 
    oT = add_params(T, U)       # output type
    X = promote_eltype(X, T)    # working type
    wT = promote_eltype(oT, U)  # working type
    Xo = convert(oT, convert_rotation(wT, X))
end

    




#######################################
# implementations
#######################################


@doc """
function to convert the input quaternion to a rotation matrix
""" ->
function quat_to_rot(q::Quaternion)

    # get rotation matrix from quaternion
    xx = q.v1 * q.v1
    yy = q.v2 * q.v2
    zz = q.v3 * q.v3
    xy = q.v1 * q.v2
    zw = q.s  * q.v3
    xz = q.v1 * q.v3
    yw = q.v2 * q.s
    yz = q.v2 * q.v3
    xw = q.s  * q.v1

    # initialize rotation part
    return @fsa([1 - 2 * (yy + zz)    2 * (xy - zw)       2 * (xz + yw);
                 2 * (xy + zw)        1 - 2 * (xx + zz)   2 * (yz - xw); 
                 2 * (xz - yw)        2 * (yz + xw)       1 - 2 * (xx + yy)])    

end

"""
 a mutable version filling the 3x3 buffer R using the 4 element vector q

 quat_to_rot!(R, q)   

"""
function quat_to_rot!(R, q)

    # get rotation matrix from quaternion
    xx = q[2] * q[2]
    yy = q[3] * q[3]
    zz = q[4] * q[4]
    xy = q[2] * q[3]
    zw = q[1] * q[4]
    xz = q[2] * q[4]
    yw = q[3] * q[1]
    yz = q[3] * q[4]
    xw = q[1] * q[2]

    # fill the rotation matrix
    R[1,1], R[1,2], R[1,3] = 1 - 2 * (yy + zz),        2 * (xy - zw),       2 * (xz + yw)
    R[2,1], R[2,2], R[2,3] =     2 * (xy + zw),    1 - 2 * (xx + zz),       2 * (yz - xw)
    R[3,1], R[3,2], R[3,3] =     2 * (xz - yw),        2 * (yz + xw),   1 - 2 * (xx + yy)
                

end



@doc """
function to convert the input rotation matrix into a quaternion.
""" ->
function rot_to_quat(R::RotMatrix)
    q = Quaternion(sqrt(abs(1  + R[1,1] + R[2,2] + R[3,3])) / 2.0, 
                   copysign(sqrt(abs(1 + R[1,1] - R[2,2] - R[3,3]))/2, R[3,2] - R[2,3]),
                   copysign(sqrt(abs(1 - R[1,1] + R[2,2] - R[3,3]))/2, R[1,3] - R[3,1]),
                   copysign(sqrt(abs(1 - R[1,1] - R[2,2] + R[3,3]))/2, R[2,1] - R[1,2]),
                   true)

    return q

end

"""
 a mutable version filling the 4 element buffer q using the 3x3 matrix R

 rot_to_quat!(q, R)

"""
function rot_to_quat!(q, R)
    q[1] = sqrt(abs(1  + R[1,1] + R[2,2] + R[3,3])) / 2.0
    q[2] = copysign(sqrt(abs(1 + R[1,1] - R[2,2] - R[3,3]))/2, R[3,2] - R[2,3])
    q[3] = copysign(sqrt(abs(1 - R[1,1] + R[2,2] - R[3,3]))/2, R[1,3] - R[3,1])
    q[4] = copysign(sqrt(abs(1 - R[1,1] - R[2,2] + R[3,3]))/2, R[2,1] - R[1,2])
end

@doc """
function to find the stereographic projection of a quaternion (see ?SpQuat).
""" ->
function quat_to_spquat(q::Quaternion)

    # N.B. there's a singularity at q = [-1.0,0,0,0] - avoid it
    # alpha2 = (q.s > -.99999999) ? (1 - q.s) / (1 + q.s) : 0.0
    q *= sign(q.s)  # this way will keep the 2norm of the output spquat <= 1, I think they're better behaved in this range.  Also avoids the singularity
    alpha2 = (1 - q.s) / (1 + q.s)
    spq = SpQuat(q.v1 * (alpha2 + 1)/2,  q.v2 * (alpha2 + 1)/2, q.v3 * (alpha2 + 1)/2)
    return spq

end

"""
 a mutable version filling the 3 element buffer spq using the 4 element vector q

 quat_to_spquat!(spq, q)

"""
function quat_to_spquat!(spq, q)
    qs = sign(q[1])  # this way will keep the 2norm of the output spquat <= 1, I think they're better behaved in this range.  Also avoids the singularity
    alpha2 = (1 - qs * q[1]) / (1 + qs * q[1])
    spq[1], spq[2], spq[3]  = qs * q[2] * (alpha2 + 1)/2,  qs * q[3] * (alpha2 + 1)/2, qs * q[4] * (alpha2 + 1)/2
end


@doc """
function to convert a the stereographic projection of a unit quaternion back into a quaternion (see ?SpQuat).
""" ->
function spquat_to_quat(spq::SpQuat)
    alpha2 = sum(spq.x .* spq.x + spq.y .* spq.y + spq.z .* spq.z)
    q = Quaternion((1.0-alpha2) / (alpha2 + 1), 2*spq.x   / (alpha2 + 1),   2*spq.y  / (alpha2 + 1), 2*spq.z / (alpha2 + 1), true)  
    q *= sign(q.s)
    return q::Quaternion{Float64}
end

"""
 a mutable version filling the 4 element vector q using the 3 element buffer spq

 sqquat_to_quat!(q, spq)

"""
function spquat_to_quat!(q, spq)
    alpha2 = spq[1] * spq[1] + spq[2] * spq[2] + spq[3] * spq[3]
    q[1] = (1.0-alpha2) / (alpha2 + 1)
    qs = sign(q[1])
    q[1] *= qs
    q[2], q[3], q[4] = qs * 2 * spq[1] / (alpha2 + 1),   qs * 2 * spq[2]  / (alpha2 + 1),  qs * 2 * spq[3] / (alpha2 + 1)  
end


# define an spquat -> rot because its common
"""
 a mutable version filling the 3 x 3 buffer R using the 3 element buffer spq

 sqquat_to_rot!(R, spq)

"""
function spquat_to_rot!(R, spq)

    # use the R buffer to store the quaternion
    spquat_to_quat!(R, spq)
    q = (R[1], R[2], R[3], R[4])  # and grab it

    # now fill R using the quaternion we just grabbed
    quat_to_rot!(R, q)
    
end




@doc """
function to convert an arbitrary axis rotation to a quaternion
""" ->
function arbaxis_to_quat(ar::AngleAxis)

    qtheta = cos(ar.theta / 2.0)

    # scale the axis _to_ make the quaternion normalized
    axis = Vec(ar.axis_x, ar.axis_y, ar.axis_z)
    s =  sqrt((1.0 - qtheta*qtheta) / sum(axis .* axis))
    
    #keep q[1] > 0 
    qsign = sign(qtheta)
    s *= qsign
    return Quaternion(qsign * qtheta, s * ar.axis_x, s * ar.axis_y, s * ar.axis_z)

end

@doc """
function to convert an a quaternion to an arbitrary axis rotation
""" ->
function quat_to_arbaxis(q::Quaternion)
    theta = rot_angle(q)
    s = sqrt(q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3)
    hasnorm = s > 1e-12
    return hasnorm ? AngleAxis(theta, q.v1 / s, q.v2 / s, q.v3 / s) : AngleAxis(theta, 1.0, 0.0, 0.0)
end




@doc """
Function to test if R is a valid rotation matrix.  This checks I - R * R'
""" ->
function valid_rotation(R, thresh::Real=1e-9)

    if !isa(R, Mat)
        R = convert_rotation(RotMatrix, R)
    end
    res = R * R'
    et = eltype(res)
    dX = res - Mat{size(res,1), size(res,1), et}(eye(et, size(res,1)))
    return vecnorm(dX) < thresh

end


