# File created: 08-Mar-2016 15:00:09
# using Matlab 2013b with Symbolic Math Toolbox version 5.11
# Working for Matrix -> EulerAngles follows the (hopefully) numerical stable method described in: https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2012/07/euler-angles1.pdf
# Note though that unlike the reference this package uses active (right handed) rotations



############################################################
# Rotation matrix converions for Euler Angles (tait-byran)
############################################################


#
# EulerZYX
#

euler_to_rot{T}(::Type{RotMatrix}, ea::EulerAngles{EulerZYX,T}) = euler_to_rot(ea::EulerAngles{EulerZYX,T})

function euler_to_rot{T}(ea::EulerAngles{EulerZYX,T})


    ct1, st1 = cos(ea.theta_x), sin(ea.theta_x)
    ct2, st2 = cos(ea.theta_y), sin(ea.theta_y)
    ct3, st3 = cos(ea.theta_z), sin(ea.theta_z)

    @fsa([ ct2*ct3  ct3*st1*st2 - ct1*st3    st1*st3 + ct1*ct3*st2  ;
           ct2*st3  ct1*ct3 + st1*st2*st3    ct1*st2*st3 - ct3*st1  ;
           -st2     ct2*st1                  ct1*ct2                ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::EulerAngles{EulerZYX,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_x)), sin(wT(ea.theta_x))
    ct2, st2 = cos(wT(ea.theta_y)), sin(wT(ea.theta_y))
    ct3, st3 = cos(wT(ea.theta_z)), sin(wT(ea.theta_z))

    @fsa([ T(ct2*ct3)  T(ct3*st1*st2 - ct1*st3)    T(st1*st3 + ct1*ct3*st2)  ;
           T(ct2*st3)  T(ct1*ct3 + st1*st2*st3)    T(ct1*st2*st3 - ct3*st1)  ;
           T(-st2)     T(ct2*st1)                  T(ct1*ct2)                ])

end


rot_to_euler{T}(::Type{EulerAngles{EulerZYX}}, R::RotMatrix{T}) = rot_to_euler(EulerAngles{EulerZYX, T}, R)

function rot_to_euler{T}(::Type{EulerAngles{EulerZYX,T}}, R::RotMatrix)

    t1 = atan2(R[2, 1], R[1, 1])
    ct1, st1 = cos(t1), sin(t1)

    EulerAngles{EulerZYX,T}(
        T(atan2(R[1, 3]*st1 - R[2, 3]*ct1, R[2, 2]*ct1 - R[1, 2]*st1)),
        T(atan2(-R[3, 1], (R[3, 2] * R[3, 2] + R[3, 3] * R[3, 3])^(1/2))),
        T(t1)
        )
end


#
# EulerZXY
#

euler_to_rot{T}(::Type{RotMatrix}, ea::EulerAngles{EulerZXY,T}) = euler_to_rot(ea::EulerAngles{EulerZXY,T})

function euler_to_rot{T}(ea::EulerAngles{EulerZXY,T})


    ct1, st1 = cos(ea.theta_x), sin(ea.theta_x)
    ct2, st2 = cos(ea.theta_y), sin(ea.theta_y)
    ct3, st3 = cos(ea.theta_z), sin(ea.theta_z)

    @fsa([ ct2*ct3 - st1*st2*st3  -ct1*st3    ct3*st2 + ct2*st1*st3  ;
           ct2*st3 + ct3*st1*st2  ct1*ct3     st2*st3 - ct2*ct3*st1  ;
           -ct1*st2               st1         ct1*ct2                ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::EulerAngles{EulerZXY,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_x)), sin(wT(ea.theta_x))
    ct2, st2 = cos(wT(ea.theta_y)), sin(wT(ea.theta_y))
    ct3, st3 = cos(wT(ea.theta_z)), sin(wT(ea.theta_z))

    @fsa([ T(ct2*ct3 - st1*st2*st3)  T(-ct1*st3)    T(ct3*st2 + ct2*st1*st3)  ;
           T(ct2*st3 + ct3*st1*st2)  T(ct1*ct3)     T(st2*st3 - ct2*ct3*st1)  ;
           T(-ct1*st2)               T(st1)         T(ct1*ct2)                ])

end


rot_to_euler{T}(::Type{EulerAngles{EulerZXY}}, R::RotMatrix{T}) = rot_to_euler(EulerAngles{EulerZXY, T}, R)

function rot_to_euler{T}(::Type{EulerAngles{EulerZXY,T}}, R::RotMatrix)

    t1 = atan2(-R[1, 2], R[2, 2])
    ct1, st1 = cos(t1), sin(t1)

    EulerAngles{EulerZXY,T}(
        T(atan2(R[3, 2], (R[3, 1] * R[3, 1] + R[3, 3] * R[3, 3])^(1/2))),
        T(atan2(R[1, 3]*ct1 + R[2, 3]*st1, R[1, 1]*ct1 + R[2, 1]*st1)),
        T(t1)
        )
end


#
# EulerYZX
#

euler_to_rot{T}(::Type{RotMatrix}, ea::EulerAngles{EulerYZX,T}) = euler_to_rot(ea::EulerAngles{EulerYZX,T})

function euler_to_rot{T}(ea::EulerAngles{EulerYZX,T})


    ct1, st1 = cos(ea.theta_x), sin(ea.theta_x)
    ct2, st2 = cos(ea.theta_y), sin(ea.theta_y)
    ct3, st3 = cos(ea.theta_z), sin(ea.theta_z)

    @fsa([ ct2*ct3   st1*st2 - ct1*ct2*st3    ct1*st2 + ct2*st1*st3  ;
           st3       ct1*ct3                  -ct3*st1               ;
           -ct3*st2  ct2*st1 + ct1*st2*st3    ct1*ct2 - st1*st2*st3  ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::EulerAngles{EulerYZX,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_x)), sin(wT(ea.theta_x))
    ct2, st2 = cos(wT(ea.theta_y)), sin(wT(ea.theta_y))
    ct3, st3 = cos(wT(ea.theta_z)), sin(wT(ea.theta_z))

    @fsa([ T(ct2*ct3)   T(st1*st2 - ct1*ct2*st3)    T(ct1*st2 + ct2*st1*st3)  ;
           T(st3)       T(ct1*ct3)                  T(-ct3*st1)               ;
           T(-ct3*st2)  T(ct2*st1 + ct1*st2*st3)    T(ct1*ct2 - st1*st2*st3)  ])

end


rot_to_euler{T}(::Type{EulerAngles{EulerYZX}}, R::RotMatrix{T}) = rot_to_euler(EulerAngles{EulerYZX, T}, R)

function rot_to_euler{T}(::Type{EulerAngles{EulerYZX,T}}, R::RotMatrix)

    t1 = atan2(-R[3, 1], R[1, 1])
    ct1, st1 = cos(t1), sin(t1)

    EulerAngles{EulerYZX,T}(
        T(atan2(R[3, 2]*ct1 + R[1, 2]*st1, R[3, 3]*ct1 + R[1, 3]*st1)),
        T(t1),
        T(atan2(R[2, 1], (R[2, 2] * R[2, 2] + R[2, 3] * R[2, 3])^(1/2)))
        )
end


#
# EulerYXZ
#

euler_to_rot{T}(::Type{RotMatrix}, ea::EulerAngles{EulerYXZ,T}) = euler_to_rot(ea::EulerAngles{EulerYXZ,T})

function euler_to_rot{T}(ea::EulerAngles{EulerYXZ,T})


    ct1, st1 = cos(ea.theta_x), sin(ea.theta_x)
    ct2, st2 = cos(ea.theta_y), sin(ea.theta_y)
    ct3, st3 = cos(ea.theta_z), sin(ea.theta_z)

    @fsa([ ct2*ct3 + st1*st2*st3  ct3*st1*st2 - ct2*st3    ct1*st2  ;
           ct1*st3                ct1*ct3                  -st1     ;
           ct2*st1*st3 - ct3*st2  st2*st3 + ct2*ct3*st1    ct1*ct2  ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::EulerAngles{EulerYXZ,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_x)), sin(wT(ea.theta_x))
    ct2, st2 = cos(wT(ea.theta_y)), sin(wT(ea.theta_y))
    ct3, st3 = cos(wT(ea.theta_z)), sin(wT(ea.theta_z))

    @fsa([ T(ct2*ct3 + st1*st2*st3)  T(ct3*st1*st2 - ct2*st3)    T(ct1*st2)  ;
           T(ct1*st3)                T(ct1*ct3)                  T(-st1)     ;
           T(ct2*st1*st3 - ct3*st2)  T(st2*st3 + ct2*ct3*st1)    T(ct1*ct2)  ])

end


rot_to_euler{T}(::Type{EulerAngles{EulerYXZ}}, R::RotMatrix{T}) = rot_to_euler(EulerAngles{EulerYXZ, T}, R)

function rot_to_euler{T}(::Type{EulerAngles{EulerYXZ,T}}, R::RotMatrix)

    t1 = atan2(R[1, 3], R[3, 3])
    ct1, st1 = cos(t1), sin(t1)

    EulerAngles{EulerYXZ,T}(
        T(atan2(-R[2, 3], (R[2, 1] * R[2, 1] + R[2, 2] * R[2, 2])^(1/2))),
        T(t1),
        T(atan2(R[3, 2]*st1 - R[1, 2]*ct1, R[1, 1]*ct1 - R[3, 1]*st1))
        )
end


#
# EulerXYZ
#

euler_to_rot{T}(::Type{RotMatrix}, ea::EulerAngles{EulerXYZ,T}) = euler_to_rot(ea::EulerAngles{EulerXYZ,T})

function euler_to_rot{T}(ea::EulerAngles{EulerXYZ,T})


    ct1, st1 = cos(ea.theta_x), sin(ea.theta_x)
    ct2, st2 = cos(ea.theta_y), sin(ea.theta_y)
    ct3, st3 = cos(ea.theta_z), sin(ea.theta_z)

    @fsa([ ct2*ct3                -ct2*st3                 st2       ;
           ct1*st3 + ct3*st1*st2  ct1*ct3 - st1*st2*st3    -ct2*st1  ;
           st1*st3 - ct1*ct3*st2  ct3*st1 + ct1*st2*st3    ct1*ct2   ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::EulerAngles{EulerXYZ,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_x)), sin(wT(ea.theta_x))
    ct2, st2 = cos(wT(ea.theta_y)), sin(wT(ea.theta_y))
    ct3, st3 = cos(wT(ea.theta_z)), sin(wT(ea.theta_z))

    @fsa([ T(ct2*ct3)                T(-ct2*st3)                 T(st2)       ;
           T(ct1*st3 + ct3*st1*st2)  T(ct1*ct3 - st1*st2*st3)    T(-ct2*st1)  ;
           T(st1*st3 - ct1*ct3*st2)  T(ct3*st1 + ct1*st2*st3)    T(ct1*ct2)   ])

end


rot_to_euler{T}(::Type{EulerAngles{EulerXYZ}}, R::RotMatrix{T}) = rot_to_euler(EulerAngles{EulerXYZ, T}, R)

function rot_to_euler{T}(::Type{EulerAngles{EulerXYZ,T}}, R::RotMatrix)

    t1 = atan2(-R[2, 3], R[3, 3])
    ct1, st1 = cos(t1), sin(t1)

    EulerAngles{EulerXYZ,T}(
        T(t1),
        T(atan2(R[1, 3], (R[1, 1] * R[1, 1] + R[1, 2] * R[1, 2])^(1/2))),
        T(atan2(R[2, 1]*ct1 + R[3, 1]*st1, R[2, 2]*ct1 + R[3, 2]*st1))
        )
end


#
# EulerXZY
#

euler_to_rot{T}(::Type{RotMatrix}, ea::EulerAngles{EulerXZY,T}) = euler_to_rot(ea::EulerAngles{EulerXZY,T})

function euler_to_rot{T}(ea::EulerAngles{EulerXZY,T})


    ct1, st1 = cos(ea.theta_x), sin(ea.theta_x)
    ct2, st2 = cos(ea.theta_y), sin(ea.theta_y)
    ct3, st3 = cos(ea.theta_z), sin(ea.theta_z)

    @fsa([ ct2*ct3                -st3       ct3*st2                ;
           st1*st2 + ct1*ct2*st3  ct1*ct3    ct1*st2*st3 - ct2*st1  ;
           ct2*st1*st3 - ct1*st2  ct3*st1    ct1*ct2 + st1*st2*st3  ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::EulerAngles{EulerXZY,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_x)), sin(wT(ea.theta_x))
    ct2, st2 = cos(wT(ea.theta_y)), sin(wT(ea.theta_y))
    ct3, st3 = cos(wT(ea.theta_z)), sin(wT(ea.theta_z))

    @fsa([ T(ct2*ct3)                T(-st3)       T(ct3*st2)                ;
           T(st1*st2 + ct1*ct2*st3)  T(ct1*ct3)    T(ct1*st2*st3 - ct2*st1)  ;
           T(ct2*st1*st3 - ct1*st2)  T(ct3*st1)    T(ct1*ct2 + st1*st2*st3)  ])

end


rot_to_euler{T}(::Type{EulerAngles{EulerXZY}}, R::RotMatrix{T}) = rot_to_euler(EulerAngles{EulerXZY, T}, R)

function rot_to_euler{T}(::Type{EulerAngles{EulerXZY,T}}, R::RotMatrix)

    t1 = atan2(R[3, 2], R[2, 2])
    ct1, st1 = cos(t1), sin(t1)

    EulerAngles{EulerXZY,T}(
        T(t1),
        T(atan2(R[2, 1]*st1 - R[3, 1]*ct1, R[3, 3]*ct1 - R[2, 3]*st1)),
        T(atan2(-R[1, 2], (R[1, 1] * R[1, 1] + R[1, 3] * R[1, 3])^(1/2)))
        )
end


############################################################
# Rotation matrix converions for Euler Angles (proper-euler)
############################################################


#
# EulerZYZ
#

euler_to_rot{T}(::Type{RotMatrix}, ea::ProperEulerAngles{EulerZYZ,T}) = euler_to_rot(ea::ProperEulerAngles{EulerZYZ,T})

function euler_to_rot{T}(ea::ProperEulerAngles{EulerZYZ,T})


    ct1, st1 = cos(ea.theta_1), sin(ea.theta_1)
    ct2, st2 = cos(ea.theta_2), sin(ea.theta_2)
    ct3, st3 = cos(ea.theta_3), sin(ea.theta_3)

    @fsa([ ct1*ct2*ct3 - st1*st3  -ct3*st1 - ct1*ct2*st3    ct1*st2  ;
           ct1*st3 + ct2*ct3*st1  ct1*ct3 - ct2*st1*st3     st1*st2  ;
           -ct3*st2               st2*st3                   ct2      ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::ProperEulerAngles{EulerZYZ,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_1)), sin(wT(ea.theta_1))
    ct2, st2 = cos(wT(ea.theta_2)), sin(wT(ea.theta_2))
    ct3, st3 = cos(wT(ea.theta_3)), sin(wT(ea.theta_3))

    @fsa([ T(ct1*ct2*ct3 - st1*st3)  T(-ct3*st1 - ct1*ct2*st3)    T(ct1*st2)  ;
           T(ct1*st3 + ct2*ct3*st1)  T(ct1*ct3 - ct2*st1*st3)     T(st1*st2)  ;
           T(-ct3*st2)               T(st2*st3)                   T(ct2)      ])

end


rot_to_euler{T}(::Type{ProperEulerAngles{EulerZYZ}}, R::RotMatrix{T}) = rot_to_euler(ProperEulerAngles{EulerZYZ, T}, R)

function rot_to_euler{T}(::Type{ProperEulerAngles{EulerZYZ,T}}, R::RotMatrix)

    t1 = atan2(R[2, 3], R[1, 3])
    ct1, st1 = cos(t1), sin(t1)

    ProperEulerAngles{EulerZYZ,T}(
        T(t1),
        T(atan2((R[3, 1] * R[3, 1] + R[3, 2] * R[3, 2])^(1/2), R[3, 3])),
        T(atan2(R[2, 1]*ct1 - R[1, 1]*st1, R[2, 2]*ct1 - R[1, 2]*st1))
        )
end


#
# EulerZXZ
#

euler_to_rot{T}(::Type{RotMatrix}, ea::ProperEulerAngles{EulerZXZ,T}) = euler_to_rot(ea::ProperEulerAngles{EulerZXZ,T})

function euler_to_rot{T}(ea::ProperEulerAngles{EulerZXZ,T})


    ct1, st1 = cos(ea.theta_1), sin(ea.theta_1)
    ct2, st2 = cos(ea.theta_2), sin(ea.theta_2)
    ct3, st3 = cos(ea.theta_3), sin(ea.theta_3)

    @fsa([ ct1*ct3 - ct2*st1*st3  -ct1*st3 - ct2*ct3*st1    st1*st2   ;
           ct3*st1 + ct1*ct2*st3  ct1*ct2*ct3 - st1*st3     -ct1*st2  ;
           st2*st3                ct3*st2                   ct2       ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::ProperEulerAngles{EulerZXZ,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_1)), sin(wT(ea.theta_1))
    ct2, st2 = cos(wT(ea.theta_2)), sin(wT(ea.theta_2))
    ct3, st3 = cos(wT(ea.theta_3)), sin(wT(ea.theta_3))

    @fsa([ T(ct1*ct3 - ct2*st1*st3)  T(-ct1*st3 - ct2*ct3*st1)    T(st1*st2)   ;
           T(ct3*st1 + ct1*ct2*st3)  T(ct1*ct2*ct3 - st1*st3)     T(-ct1*st2)  ;
           T(st2*st3)                T(ct3*st2)                   T(ct2)       ])

end


rot_to_euler{T}(::Type{ProperEulerAngles{EulerZXZ}}, R::RotMatrix{T}) = rot_to_euler(ProperEulerAngles{EulerZXZ, T}, R)

function rot_to_euler{T}(::Type{ProperEulerAngles{EulerZXZ,T}}, R::RotMatrix)

    t1 = atan2(R[1, 3], (-R[2, 3] + eps()) - eps())  # TODO: handle denormal numbers better, as atan2(0,0) != atan2(0,-0)
    ct1, st1 = cos(t1), sin(t1)

    ProperEulerAngles{EulerZXZ,T}(
        T(t1),
        T(atan2((R[3, 1] * R[3, 1] + R[3, 2] * R[3, 2])^(1/2), R[3, 3])),
        T(atan2(- R[1, 2]*ct1 - R[2, 2]*st1, R[1, 1]*ct1 + R[2, 1]*st1))
        )
end


#
# EulerYZY
#

euler_to_rot{T}(::Type{RotMatrix}, ea::ProperEulerAngles{EulerYZY,T}) = euler_to_rot(ea::ProperEulerAngles{EulerYZY,T})

function euler_to_rot{T}(ea::ProperEulerAngles{EulerYZY,T})


    ct1, st1 = cos(ea.theta_1), sin(ea.theta_1)
    ct2, st2 = cos(ea.theta_2), sin(ea.theta_2)
    ct3, st3 = cos(ea.theta_3), sin(ea.theta_3)

    @fsa([ ct1*ct2*ct3 - st1*st3   -ct1*st2    ct3*st1 + ct1*ct2*st3  ;
           ct3*st2                 ct2         st2*st3                ;
           -ct1*st3 - ct2*ct3*st1  st1*st2     ct1*ct3 - ct2*st1*st3  ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::ProperEulerAngles{EulerYZY,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_1)), sin(wT(ea.theta_1))
    ct2, st2 = cos(wT(ea.theta_2)), sin(wT(ea.theta_2))
    ct3, st3 = cos(wT(ea.theta_3)), sin(wT(ea.theta_3))

    @fsa([ T(ct1*ct2*ct3 - st1*st3)   T(-ct1*st2)    T(ct3*st1 + ct1*ct2*st3)  ;
           T(ct3*st2)                 T(ct2)         T(st2*st3)                ;
           T(-ct1*st3 - ct2*ct3*st1)  T(st1*st2)     T(ct1*ct3 - ct2*st1*st3)  ])

end


rot_to_euler{T}(::Type{ProperEulerAngles{EulerYZY}}, R::RotMatrix{T}) = rot_to_euler(ProperEulerAngles{EulerYZY, T}, R)

function rot_to_euler{T}(::Type{ProperEulerAngles{EulerYZY,T}}, R::RotMatrix)

    t1 = atan2(R[3, 2], (-R[1, 2] + eps()) - eps())  # TODO: handle denormal numbers better, as atan2(0,0) != atan2(0,-0)
    ct1, st1 = cos(t1), sin(t1)

    ProperEulerAngles{EulerYZY,T}(
        T(t1),
        T(atan2((R[2, 1] * R[2, 1] + R[2, 3] * R[2, 3])^(1/2), R[2, 2])),
        T(atan2(- R[3, 1]*ct1 - R[1, 1]*st1, R[3, 3]*ct1 + R[1, 3]*st1))
        )
end


#
# EulerYXY
#

euler_to_rot{T}(::Type{RotMatrix}, ea::ProperEulerAngles{EulerYXY,T}) = euler_to_rot(ea::ProperEulerAngles{EulerYXY,T})

function euler_to_rot{T}(ea::ProperEulerAngles{EulerYXY,T})


    ct1, st1 = cos(ea.theta_1), sin(ea.theta_1)
    ct2, st2 = cos(ea.theta_2), sin(ea.theta_2)
    ct3, st3 = cos(ea.theta_3), sin(ea.theta_3)

    @fsa([ ct1*ct3 - ct2*st1*st3   st1*st2    ct1*st3 + ct2*ct3*st1  ;
           st2*st3                 ct2        -ct3*st2               ;
           -ct3*st1 - ct1*ct2*st3  ct1*st2    ct1*ct2*ct3 - st1*st3  ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::ProperEulerAngles{EulerYXY,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_1)), sin(wT(ea.theta_1))
    ct2, st2 = cos(wT(ea.theta_2)), sin(wT(ea.theta_2))
    ct3, st3 = cos(wT(ea.theta_3)), sin(wT(ea.theta_3))

    @fsa([ T(ct1*ct3 - ct2*st1*st3)   T(st1*st2)    T(ct1*st3 + ct2*ct3*st1)  ;
           T(st2*st3)                 T(ct2)        T(-ct3*st2)               ;
           T(-ct3*st1 - ct1*ct2*st3)  T(ct1*st2)    T(ct1*ct2*ct3 - st1*st3)  ])

end


rot_to_euler{T}(::Type{ProperEulerAngles{EulerYXY}}, R::RotMatrix{T}) = rot_to_euler(ProperEulerAngles{EulerYXY, T}, R)

function rot_to_euler{T}(::Type{ProperEulerAngles{EulerYXY,T}}, R::RotMatrix)

    t1 = atan2(R[1, 2], R[3, 2])
    ct1, st1 = cos(t1), sin(t1)

    ProperEulerAngles{EulerYXY,T}(
        T(t1),
        T(atan2((R[2, 1] * R[2, 1] + R[2, 3] * R[2, 3])^(1/2), R[2, 2])),
        T(atan2(R[1, 3]*ct1 - R[3, 3]*st1, R[1, 1]*ct1 - R[3, 1]*st1))
        )
end


#
# EulerXYX
#

euler_to_rot{T}(::Type{RotMatrix}, ea::ProperEulerAngles{EulerXYX,T}) = euler_to_rot(ea::ProperEulerAngles{EulerXYX,T})

function euler_to_rot{T}(ea::ProperEulerAngles{EulerXYX,T})


    ct1, st1 = cos(ea.theta_1), sin(ea.theta_1)
    ct2, st2 = cos(ea.theta_2), sin(ea.theta_2)
    ct3, st3 = cos(ea.theta_3), sin(ea.theta_3)

    @fsa([ ct2       st2*st3                  ct3*st2                 ;
           st1*st2   ct1*ct3 - ct2*st1*st3    -ct1*st3 - ct2*ct3*st1  ;
           -ct1*st2  ct3*st1 + ct1*ct2*st3    ct1*ct2*ct3 - st1*st3   ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::ProperEulerAngles{EulerXYX,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_1)), sin(wT(ea.theta_1))
    ct2, st2 = cos(wT(ea.theta_2)), sin(wT(ea.theta_2))
    ct3, st3 = cos(wT(ea.theta_3)), sin(wT(ea.theta_3))

    @fsa([ T(ct2)       T(st2*st3)                  T(ct3*st2)                 ;
           T(st1*st2)   T(ct1*ct3 - ct2*st1*st3)    T(-ct1*st3 - ct2*ct3*st1)  ;
           T(-ct1*st2)  T(ct3*st1 + ct1*ct2*st3)    T(ct1*ct2*ct3 - st1*st3)   ])

end


rot_to_euler{T}(::Type{ProperEulerAngles{EulerXYX}}, R::RotMatrix{T}) = rot_to_euler(ProperEulerAngles{EulerXYX, T}, R)

function rot_to_euler{T}(::Type{ProperEulerAngles{EulerXYX,T}}, R::RotMatrix)

    t1 = atan2(R[2, 1], (-R[3, 1] + eps()) - eps())  # TODO: handle denormal numbers better, as atan2(0,0) != atan2(0,-0)
    ct1, st1 = cos(t1), sin(t1)

    ProperEulerAngles{EulerXYX,T}(
        T(t1),
        T(atan2((R[1, 2] * R[1, 2] + R[1, 3] * R[1, 3])^(1/2), R[1, 1])),
        T(atan2(- R[2, 3]*ct1 - R[3, 3]*st1, R[2, 2]*ct1 + R[3, 2]*st1))
        )
end


#
# EulerXZX
#

euler_to_rot{T}(::Type{RotMatrix}, ea::ProperEulerAngles{EulerXZX,T}) = euler_to_rot(ea::ProperEulerAngles{EulerXZX,T})

function euler_to_rot{T}(ea::ProperEulerAngles{EulerXZX,T})


    ct1, st1 = cos(ea.theta_1), sin(ea.theta_1)
    ct2, st2 = cos(ea.theta_2), sin(ea.theta_2)
    ct3, st3 = cos(ea.theta_3), sin(ea.theta_3)

    @fsa([ ct2      -ct3*st2                 st2*st3                 ;
           ct1*st2  ct1*ct2*ct3 - st1*st3    -ct3*st1 - ct1*ct2*st3  ;
           st1*st2  ct1*st3 + ct2*ct3*st1    ct1*ct3 - ct2*st1*st3   ])

end

function euler_to_rot{T,U}(::Type{RotMatrix{T}}, ea::ProperEulerAngles{EulerXZX,U})

    wT = promote_type(T, U)
    ct1, st1 = cos(wT(ea.theta_1)), sin(wT(ea.theta_1))
    ct2, st2 = cos(wT(ea.theta_2)), sin(wT(ea.theta_2))
    ct3, st3 = cos(wT(ea.theta_3)), sin(wT(ea.theta_3))

    @fsa([ T(ct2)      T(-ct3*st2)                 T(st2*st3)                 ;
           T(ct1*st2)  T(ct1*ct2*ct3 - st1*st3)    T(-ct3*st1 - ct1*ct2*st3)  ;
           T(st1*st2)  T(ct1*st3 + ct2*ct3*st1)    T(ct1*ct3 - ct2*st1*st3)   ])

end


rot_to_euler{T}(::Type{ProperEulerAngles{EulerXZX}}, R::RotMatrix{T}) = rot_to_euler(ProperEulerAngles{EulerXZX, T}, R)

function rot_to_euler{T}(::Type{ProperEulerAngles{EulerXZX,T}}, R::RotMatrix)

    t1 = atan2(R[3, 1], R[2, 1])  
    ct1, st1 = cos(t1), sin(t1)

    ProperEulerAngles{EulerXZX,T}(
        T(t1),
        T(atan2((R[1, 2] * R[1, 2] + R[1, 3] * R[1, 3])^(1/2), R[1, 1])),
        T(atan2(R[3, 2]*ct1 - R[2, 2]*st1, R[3, 3]*ct1 - R[2, 3]*st1))
        )
end








