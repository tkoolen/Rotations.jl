
@doc """
function to convert the input euler angles into a rotation matrix.  Returns a RotMatrix
""" ->
function eulertorot{ORD,T}(att::EulerAngles{ORD,T})
	error("Euler angle order to rotation: $(ORD) is not implemented")
end

# known transforms
function eulertorot{T}(att::EulerAngles{EulerZXY, T})

    # Euler angles to rotation matrix
    si = sin(att[1])
    ci = cos(att[1])
    sj = sin(att[2])
    cj = cos(att[2])
    sk = sin(att[3])
    ck = cos(att[3])

    R = @fsa ([   cj*ck-si*sj*sk   -ci*sk   sj*ck+si*cj*sk;
                  cj*sk+si*sj*ck    ci*ck   sj*sk-si*cj*ck;
                 -ci*sj             si      ci*cj])

    return R

end



@doc """
function to convert the rotation matrix into Euler angles

Args:

	R      - a rotation matrix

	format - a value symbol representing the desired euler angle format Default: Val{:EulerZXY} (see: ?:EulerZXY)

""" ->
function rottoeuler{T <: EulerAngles}(::Type{T}, R)
	error("Rotation matrix to Euler angle order: $(ORD) is not implemented")
end

function rottoeuler{T <: EulerAngles{EulerZXY}}(::Type{T}, R)  
    return convert(T, asin(R[3,2]), atan2(-R[3,1],R[3,3]), atan2(-R[1,2],R[2,2]))
end
