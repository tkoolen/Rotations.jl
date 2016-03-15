@doc """

Function to convert a covariance matrix containing one rotation parameterization to another:

Inputs:

    oT        - output rotation type

    iT        - input rotation type

    cX        - covariance matrix (set = nothing to enforce zero mean in input and output)

    xbar      - mean X in the input space

    idx_in   - range with the row indexes of the rotation parameters in the input (they should be contiguous)

    idx_ou   - range with the row indexes of the rotation parameters in the output (defaults to idx_in)


Experimental - This is not a fast implementation

""" ->
function convert_covariance{T <: RotationTypes, U <: RotationTypes}(oT::Type{T}, iT::Type{U}, cX, xbar=zeros(size(cX,1)), idx_in=4:6, idx_out=idx_in)

    # make sure the element type is specified in the output type
    oT = promote_rottype(oT, eltype(iT)) 

    # difference in dimensionality of output and input rotation parameterisation
    ddims = length(idx_out) - length(idx_in)

    # generate the sigma points
    (sigma_points, weights) = getsigmapoints(cX, xbar)

    # convert them to spquat's
    if (0 == ddims)
        sigma_trans = sigma_points  # can transform in place
    else
        # need a different size on output, so create a copy
        sigma_trans = zeros(size(sigma_points,1) + ddims, size(sigma_points,2))
        sidx = 1:idx_in[1] - 1
        eidx = (idx_in[end] + 1):size(cX,1)
        sigma_trans[sidx,:] = sigma_points[sidx,:]
        sigma_trans[eidx + ddims,:] = sigma_points[eidx,:]
    end

    # and convert
    for i = 1:size(sigma_points,2)
        sigma_trans[idx_out,i] = oT(iT(sigma_points[idx_in,i]))  # type it as the input type before converting to the output type
    end

    if xbar == nothing
        xbar_out = nothing
    elseif isa(xbar, FixedSizeArray.FixedArray) 
        xbar_out = (ddims == 0) ? setindex(xbar, oT(iT(xbar[idx_in])), idx_out) : Vec(xbar(sidx)..., oT(iT(xbar[idx_in]))..., xbar(eidx)...)
    elseif (ddims == 0)
        xbar_out = copy(xbar)
        xbar_out[idx_ind] = oT(iT(xbar_out[idx_in]))
    else
        xbar_out = [xbar(sidx); oT(iT(xbar[idx_in])); xbar(eidx)]
    end    

    return covfromsigmapoints(sigma_points, weights, xbar_out), xbar_out  # return the mean as well
end



# generate sigma points for a specifified covaraince
function getsigmapoints(cX, xbar=nothing)

    if !isa(cX, AbstractMatrix)
        cX = Matrix(cX)  # chol works with Matrix not other types
    end

    # save some time and use and existing ukf approach instead of actual maths
    n = size(cX,1)

    W0 = -0.2
    R = full(chol((n/(1-W0)) * cX, Val{:L}))::Matrix{Float64} 

    # weights
    weights = zeros(2*n+1)
    weights[1] = W0
    weights[2:end] = (1.0-W0)/(2.0*n)

    # sigma points
    if xbar != nothing
        sigma_points = [xbar  (xbar .+  R)   xbar .- R]
    else
        sigma_points = [zeros(size(cX,1))  R   -R]
    end

    #TODO: Need to put in special catch if angle sigma points get too big so they don't wrap around
    return sigma_points, weights

end


# generate covariance from sigma points
function covfromsigmapoints(sigma_points, weights, xbar=nothing)

    # recover the covariance
    #cX_out = (sigmaPointWeightsObs'.* sigmaPoints) * sigmaPoints'
    if xbar != nothing
        sigma_points = sigma_points .- xbar
    end
    cov = sigma_points * (sigma_points' .* weights)  # maybe better than the above (less copy but mult across rows)

    #This is poop!! Force covariance matrix symmetric. TODO: implement sqrt
    #kf or other symmetric update scheme
    cov = 0.5*(cov + cov')
    return cov

end
