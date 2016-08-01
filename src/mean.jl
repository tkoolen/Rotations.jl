#############################################################
#
#  add extra quaternion operations for Quaternions
#  some should probably be merged back into Quaternions
#
#############################################################


"""
    mean(qvec, [method = 0])

Function to find the mean of a set of rotations.

Args:

    qvec    - Vector of `Rotation`s

    method  - if 0 (default), the returned rotation minimizes sum_i=1:n (||R(Qi) - R(Qbar)||)  with || the Frobenius norm, or
              equivalently minimizes sum_i=1:n (sin(phi_i / 2))^2, where phi_i is the angle of Qe_i, with  Qbar = Q_i x Qe_i

            - if 1, the returned rotation minimizes sum_i=1:n (phi_i ^2) where phi_i is the angle of Qe_i, with  Qbar = Q_i x Qe_i
              this should be about the same as method == 0 if phi_i are all small (note: still not implemented)
"""
function Base.mean{T}(qvec::AbstractVector{Quat{T}}, method::Integer = 0)
    if (method == 0)
        M = zeros(4, 4)
        for i = 1:length(qvec)
            q = qvec[i]
            Qi = @SVector [q.w, q.x, q.y, q.z]  # convert types to ensure we don't get quaternion multiplication
            M .+= Qi * (Qi')
        end
        evec = eigfact(Symmetric(M), 4:4)
        Qbar = Quat(evec.vectors[1], evec.vectors[2], evec.vectors[3], evec.vectors[4]) # This will renormalize the quaternion...
    else
        error("I haven't coded this")
    end

    return Qbar
end

function Base.mean{R<:Rotation}(vec::AbstractVector{R}, method::Integer = 0)
    R(mean(convert(Vector{Quat{eltype(R)}}, vec), method))
end
